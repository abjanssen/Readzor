"""Run orchestration: task generation, the worker pool, output files and summary."""

import csv
import itertools
import logging
import multiprocessing as mp
import os
import signal
import sys
import threading
import time
import traceback

from isal import igzip

from readzor import constants
from readzor import stats as S
from readzor.fastq_io import cleanup_stdin_temp_files, basename_file, common_name_parts, detect_phred_offset, find_paired_files, is_gz_file, lazy_fastq_blobs, open_fastq_writer
from readzor.progress import ProgressTracker, SharedCounterTracker, count_reads_estimated
from readzor.worker import unified_worker, worker_initializer

logger = logging.getLogger(__name__)

def input_handler(unspecified_files, unpaired_files, paired_files, interleaved_files, output_dir, threads, chunk_size, show_progress, stdout, interleaved_out, write_rejected, discard_singles, parameters, progress_counter=None, pool=None):
    """
    Top-level orchestrator: sort the inputs, open the output files, and run
    all chunks of all inputs through a single multiprocessing pool.

    Steps:
      1. Detect gzip compression of every input (fills GZIP_DETECTION).
      2. Auto-sort ``unspecified_files`` into pairs, interleaved and
         unpaired files (find_paired_files()) and merge them with the
         explicitly given groups.
      3. If show_progress is set, estimate the read count of every input
         (count_reads_estimated()) and start a ProgressTracker.
      4. Unless streaming to stdout, open the output files in output_dir.
         Each paired or interleaved input gets a unique name prefix: the
         common name of its two files (common_name_parts()), or the
         interleaved file's name, with "_2", "_3", ... appended if the
         prefix is already taken.
      5. Stream chunks of all unpaired, then interleaved, then paired
         inputs to the pool and write the results as they come back. At
         most threads * 5 chunks are in flight at once, to bound memory
         use. With parameters["ordered_output"], results are written in
         input order. In test-run mode only the first ``threads`` chunks
         of each input are processed.
      6. Close all outputs; gzip outputs that received no data get a
         valid, empty gzip member.

    With parameters["stats"] (--stats), the read statistics every worker
    returns with its chunk are merged per input into a StatsCollector,
    keyed like the returned file statistics.

    In stdout mode, only kept reads are written to stdout (for paired and
    interleaved inputs, the interleaved pairs); no FASTQ files are opened.

    Args:
        unspecified_files (list[str] | None): Files with unknown pairing
            status, to auto-detect.
        unpaired_files (list[str] | None): Explicitly provided unpaired
            FASTQ files.
        paired_files (list[tuple[str, str]] | None): Explicitly provided
            (R1, R2) file pairs.
        interleaved_files (list[str] | None): Explicitly provided
            interleaved FASTQ files.
        output_dir (str): Directory to write output files to.
        threads (int): Number of worker processes in the pool.
        chunk_size (int): Number of reads (read pairs for paired and
            interleaved inputs) per processing chunk.
        show_progress (bool): If True, renders a live progress bar to stderr.
        stdout (bool): If True, streams surviving reads to stdout instead
            of writing output files.
        interleaved_out (bool): If True, writes surviving pairs of each
            paired or interleaved input interleaved into a single file,
            instead of separate R1/R2 files.
        write_rejected (bool): If True, also writes rejected reads to file.
        discard_singles (bool): If True, no singleton files are opened for
            paired and interleaved inputs.
        parameters (dict): Run parameters.
        progress_counter (multiprocessing.Value | None): When the run is
            split over coordinator processes (readzor.split), the counter
            that processed reads are added to; the parent process draws the
            progress bar. Overrides show_progress. Defaults to None.
        pool (multiprocessing.pool.Pool | None): Worker pool to use, created
            with worker_initializer and the same ``parameters``. A coordinator
            (readzor.split) passes its own pool, so it is created once and
            reused for all of its inputs instead of once per input; it is
            left running. If None, a pool is created here and shut down at
            the end. Defaults to None.

    Returns:
        tuple[dict, stats.StatsCollector | None]:
            - Summary statistics per input. Unpaired files are keyed by
              their filepath, with "kept" and "rejected" counts. Paired and
              interleaved inputs are keyed by their name prefix, with
              "kept_pairs", "kept_R1_singletons", "kept_R2_singletons",
              "rejected_R1" and "rejected_R2" counts, since each mate is
              trimmed and filtered independently before reconciliation.
              Every entry also has "type": "unpaired", "interleaved" or
              "paired", and "input"/"output" file name strings.
            - The merged read statistics (--stats), keyed the same way, or
              None without --stats.
    """
    for file in (unspecified_files or []):
        constants.GZIP_DETECTION[file] = is_gz_file(file)
    for file in (unpaired_files or []):
        constants.GZIP_DETECTION[file] = is_gz_file(file)
    for file in (interleaved_files or []):
        constants.GZIP_DETECTION[file] = is_gz_file(file)
    for f1, f2 in (paired_files or []):
        constants.GZIP_DETECTION[f1] = is_gz_file(f1)
        constants.GZIP_DETECTION[f2] = is_gz_file(f2)
    auto_paired, auto_interleaved, auto_unpaired = find_paired_files(unspecified_files)
    unpaired = auto_unpaired + (unpaired_files or [])
    paired = auto_paired + (paired_files or [])
    interleaved = auto_interleaved + (interleaved_files or [])
    if progress_counter is not None:
        tracker = SharedCounterTracker(progress_counter)
    elif show_progress:
        for file in unpaired:
            constants.ESTIMATED_READ_COUNTS[file] = count_reads_estimated(file)
        for file in interleaved:
            constants.ESTIMATED_READ_COUNTS[file] = count_reads_estimated(file)
        for f1, f2 in paired:
            constants.ESTIMATED_READ_COUNTS[f1] = count_reads_estimated(f1)
            constants.ESTIMATED_READ_COUNTS[f2] = count_reads_estimated(f2)
        total_reads = sum(constants.ESTIMATED_READ_COUNTS.values())
        tracker = ProgressTracker(total_reads = total_reads)
    else:
        class NullTracker:
            def update(self, n): pass
            def close(self, interrupted=False): pass
            def clear_line(self): pass
        tracker = NullTracker()
    constants.ACTIVE_PROGRESS_TRACKER = tracker
    file_stats = {}
    file_writing_handles = {}
    
    out_ext = ".fastq.gz" if parameters["gzip_output"] else ".fastq"
    
    for file in unpaired:
        input_filename = os.path.basename(file)
        output_filename = basename_file(file) + ("_filtered.fastq.gz" if parameters["gzip_output"] else "_filtered.fastq")
        file_stats[file] = {
            "type": "unpaired",
            "input": input_filename,
            "output": output_filename,
            "kept": 0, 
            "rejected": 0
        }
        if not stdout:
            file_writing_handles[file] = open_fastq_writer(file, output_dir, gzip_output=parameters["gzip_output"])
        if write_rejected:
            suffix = basename_file(file) + "_rejected"
            file_writing_handles[f"{file}_rejected"] = open_fastq_writer(f"{suffix}", output_dir, gzip_output=parameters["gzip_output"])
    pair_keys = {}
    used_prefixes = {basename_file(f) for f in unpaired}

    def unique_prefix(base_prefix):
        """
        Return base_prefix, or base_prefix with the first free numeric
        suffix (_2, _3, ...) if it is already used by another input, and
        mark the returned prefix as used.
        """
        prefix, n = base_prefix, 2
        while prefix in used_prefixes:
            prefix = f"{base_prefix}_{n}"
            n += 1
        used_prefixes.add(prefix)
        return prefix

    def pair_output_keys(prefix):
        """
        Map each output role ("paired", "paired_1", "paired_2", "singles_1",
        "singles_2", "rejected", "rejected_1", "rejected_2") to its file
        handle key for one paired/interleaved input, depending on the
        interleaved-output, discard-singles and write-rejected settings.
        Returns an empty dict in stdout mode.
        """
        if stdout:
            return {}
        if interleaved_out:
            keys = {"paired": f"{prefix}_interleaved"}
        else:
            keys = {"paired_1": f"{prefix}_R1_paired", "paired_2": f"{prefix}_R2_paired"}
        if not discard_singles:
            keys["singles_1"] = f"{prefix}_R1_unpaired"
            keys["singles_2"] = f"{prefix}_R2_unpaired"
        if write_rejected:
            if interleaved_out:
                keys["rejected"] = f"{prefix}_rejected"
            else:
                keys["rejected_1"] = f"{prefix}_R1_rejected"
                keys["rejected_2"] = f"{prefix}_R2_rejected"
        return keys
    
    def get_output_files_list(prefix):
        """Build a comma-separated string of all generated output filenames for paired/interleaved inputs."""
        outputs = []
        if interleaved_out:
            outputs.append(f"{prefix}_interleaved_filtered{out_ext}")
        else:
            outputs.append(f"{prefix}_R1_paired_filtered{out_ext}")
            outputs.append(f"{prefix}_R2_paired_filtered{out_ext}")
        if not discard_singles:
            outputs.append(f"{prefix}_R1_unpaired_filtered{out_ext}")
            outputs.append(f"{prefix}_R2_unpaired_filtered{out_ext}")
        if write_rejected:
            if interleaved_out:
                outputs.append(f"{prefix}_rejected{out_ext}")
            else:
                outputs.append(f"{prefix}_R1_rejected{out_ext}")
                outputs.append(f"{prefix}_R2_rejected{out_ext}")
        return ", ".join(outputs)

    def open_pair_outputs(prefix):
        """Open every output file one paired/interleaved input needs (see pair_output_keys)."""
        for key in pair_output_keys(prefix).values():
            file_writing_handles[key] = open_fastq_writer(key, output_dir, gzip_output=parameters["gzip_output"])

    for file in interleaved:
        prefix = unique_prefix(basename_file(file))
        pair_keys[(file, file)] = prefix
        file_stats[prefix] = {
            "type": "interleaved",
            "input": os.path.basename(file),
            "output": get_output_files_list(prefix),
            "kept_pairs": 0, 
            "kept_R1_singletons": 0, 
            "kept_R2_singletons": 0, 
            "rejected_R1": 0, 
            "rejected_R2": 0
        }
        open_pair_outputs(prefix)
    for file1, file2 in paired:
        prefix = unique_prefix(common_name_parts([os.path.basename(file1), os.path.basename(file2)]))
        pair_keys[(file1, file2)] = prefix
        file_stats[prefix] = {
            "type": "paired",
            "input": f"{os.path.basename(file1)}, {os.path.basename(file2)}",
            "output": get_output_files_list(prefix),
            "kept_pairs": 0, 
            "kept_R1_singletons": 0, 
            "kept_R2_singletons": 0, 
            "rejected_R1": 0, 
            "rejected_R2": 0
        }
        open_pair_outputs(prefix)

    def unified_chunk_streamer():
        """
        Yield the tasks of all unpaired, then interleaved, then paired
        inputs. In test-run mode, each input is limited to ``threads`` tasks.
        """
        chunks_per_file = threads if parameters["testrun"] else None
        if unpaired:
            logger.info("Started processing unpaired files.")
            for file in unpaired:
                gen = generate_unpaired_tasks(
                    filepaths=[file],
                    chunk_size=chunk_size,
                    parameters=parameters,
                    filetype="unpaired",
                )
                if chunks_per_file is not None:
                    gen = itertools.islice(gen, chunks_per_file)
                yield from gen
            logger.info("Finished processing unpaired files.")
        if interleaved:
            logger.info("Started processing interleaved files.")
            for file in interleaved:
                gen = generate_unpaired_tasks(
                    filepaths=[file],
                    chunk_size=chunk_size,
                    parameters=parameters,
                    filetype="interleaved",
                )
                if chunks_per_file is not None:
                    gen = itertools.islice(gen, chunks_per_file)
                yield from gen
            logger.info("Finished processing interleaved files.")
        if paired:
            logger.info("Started processing paired files.")
            for pair in paired:
                gen = generate_paired_tasks(
                    files=[pair], chunk_size=chunk_size, parameters=parameters
                )
                if chunks_per_file is not None:
                    gen = itertools.islice(gen, chunks_per_file)
                yield from gen
            logger.info("Finished processing paired files.")
    chunk_stream = unified_chunk_streamer()

    max_in_flight = threads * 5
    backpressure = threading.Semaphore(max_in_flight)
    stop_feeding = threading.Event()
    def bounded_chunk_stream():
        """
        Yield tasks from chunk_stream, blocking while threads * 5 tasks are
        still unprocessed (the semaphore is released for each result).
        Stops early once stop_feeding is set (see stop_and_drain()).
        """
        for item in chunk_stream:
            backpressure.acquire()
            if stop_feeding.is_set():
                return
            yield item

    def stop_and_drain(results, timeout=10):
        """
        Called when the result loop is interrupted (Ctrl+C or an error).

        Stops feeding new tasks and waits for the tasks already sent to the
        workers (at most threads * 5) to come back, discarding their results.
        Afterwards every worker is idle, so the pool can be shut down
        safely. Shutting down the pool while it is still busy can hang
        forever: the task-feeding thread waits for a backpressure permit
        that is never released, or a worker is stopped in the middle of
        sending a large result through the result pipe.

        Returns:
            bool: True if all outstanding results were received; False on
                timeout or a second Ctrl+C (the pool is then shut down while
                still busy, as before).
        """
        stop_feeding.set()
        backpressure.release(max_in_flight + 1)
        logger.info("Shutdown signal received, preparing to shut down.")
        deadline = time.monotonic() + timeout
        while True:
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                logger.warning("Waited %s seconds, processing seems stuck; forcing shutdown instead.", timeout)
                return False
            try:
                results.next(timeout=remaining)
            except StopIteration:
                return True
            except mp.TimeoutError:
                logger.warning("Waited %s seconds, processing seems stuck; forcing shutdown instead.", timeout)
                return False
            except KeyboardInterrupt:
                logger.warning("Second shutdown signal received; forcing shutdown.")
                return False
            except Exception:
                continue

    def force_exit(error):
        """
        Exit immediately when stop_and_drain() could not wait for the
        workers (second Ctrl+C or timeout). The pool is still busy, so its
        normal shutdown could hang; instead the workers are killed, the
        progress bar, output files and stdin temporary files are cleaned
        up, and the process exits without running the pool's shutdown code.
        Exit status: 130 after Ctrl+C, 1 after an error.
        """
        signal.signal(signal.SIGINT, signal.SIG_IGN)
        for child in mp.active_children():
            child.kill()
        tracker.close(interrupted=True)
        constants.ACTIVE_PROGRESS_TRACKER = None
        for handle in file_writing_handles.values():
            try:
                handle.close()
            except OSError:
                pass
        cleanup_stdin_temp_files()
        if isinstance(error, KeyboardInterrupt):
            exit_code = 130
            logger.info("Readzor interrupted by user — shutting down.")
            if not parameters["verbose"]:
                print("Readzor interrupted by user — shutting down.", file=sys.stderr)
        else:
            exit_code = 1
            logger.error("Readzor stopped by an error.", exc_info=error)
            traceback.print_exception(error)
        logging.shutdown()
        sys.stdout.flush()
        sys.stderr.flush()
        os._exit(exit_code)

    collector = S.StatsCollector() if parameters["stats"] else None
    finished = False
    sys.setswitchinterval(0.0005)
    try:
        own_pool = pool is None
        if own_pool:
            pool = mp.Pool(threads, initializer=worker_initializer, initargs=(parameters,))
        try:
            submit = pool.imap if parameters["ordered_output"] else pool.imap_unordered
            results = submit(unified_worker, bounded_chunk_stream(), chunksize=1)
            try:
                for result in results:
                    backpressure.release()
                    report, result = result[-1], result[:-1]
                    if result[0] == "unpaired":
                        _, filepath, chunk_results, kept, rejected, rejected_reads = result
                        if parameters["stdout"]:
                            if chunk_results:
                                sys.stdout.buffer.write(chunk_results)
                        else:
                            if chunk_results:
                                file_writing_handles[filepath].write(chunk_results)
                            if write_rejected and rejected_reads:
                                file_writing_handles[f"{filepath}_rejected"].write(rejected_reads)
                
                        file_stats[filepath]["kept"] += kept
                        file_stats[filepath]["rejected"] += rejected
                        if collector is not None:
                            collector.add(filepath, report)
                        tracker.update(kept + rejected)
                    elif result[0] == "paired":
                        _, file1, file2, paired_out_1, paired_out_2, R1_singles_out, R2_singles_out, num_paired, num_R1_singles, num_R2_singles, rejected_1, rejected_2, rejected_R1, rejected_R2 = result
                        common_prefix = pair_keys[(file1, file2)]
                        if parameters["stdout"]:
                            if paired_out_1:
                                sys.stdout.buffer.write(paired_out_1)
                        else:
                            payloads = {
                                "paired": paired_out_1,
                                "paired_1": paired_out_1,
                                "paired_2": paired_out_2,
                                "singles_1": R1_singles_out,
                                "singles_2": R2_singles_out,
                                "rejected": rejected_R1,
                                "rejected_1": rejected_R1,
                                "rejected_2": rejected_R2,
                            }
                            for role, handle_key in pair_output_keys(common_prefix).items():
                                records = payloads[role]
                                if records:
                                    file_writing_handles[handle_key].write(records)
                        file_stats[common_prefix]["kept_pairs"] += num_paired
                        file_stats[common_prefix]["kept_R1_singletons"] += num_R1_singles
                        file_stats[common_prefix]["kept_R2_singletons"] += num_R2_singles
                        file_stats[common_prefix]["rejected_R1"] += rejected_1
                        file_stats[common_prefix]["rejected_R2"] += rejected_2
                        if collector is not None:
                            collector.add(common_prefix, report)
                        tracker.update(num_paired * 2 + num_R1_singles + num_R2_singles + rejected_1 + rejected_2)
            except BaseException as error:
                try:
                    drained = stop_and_drain(results)
                except KeyboardInterrupt:
                    drained = False
                if not drained:
                    force_exit(error)
                signal.signal(signal.SIGINT, lambda signum, frame: os._exit(130))
                raise
            finished = True
        finally:
            if own_pool:
                pool.terminate()
    finally:
        tracker.close(interrupted = not finished) 
        constants.ACTIVE_PROGRESS_TRACKER = None
        for handle in file_writing_handles.values():
            if not handle.closed:
                if parameters["gzip_output"] and handle.tell() == 0:
                    handle.write(igzip.compress(b""))
                handle.close()
    return file_stats, collector

def generate_unpaired_tasks(filepaths, chunk_size, parameters, filetype = None):
    """
    Lazily split unpaired or interleaved FASTQ files into worker tasks.

    For each file, logs basic file information, detects the Phred offset
    once, and then yields one task per chunk, so that a single global pool
    can process chunks from many files. Chunks are raw FASTQ text from
    lazy_fastq_blobs(); the workers parse them. Logs the processing time
    and rate once a file is exhausted.

      - filetype "unpaired": each task holds a blob of ``chunk_size``
        records and has type "unpaired".
      - filetype "interleaved": each task holds a blob of 2 * ``chunk_size``
        records. The task has type "paired", with file1 == file2 ==
        filepath and "interleaved" set; the worker splits the records
        alternately into R1 (records 1, 3, 5, ...) and R2 (records 2, 4,
        6, ...), so it is processed exactly like a two-file pair.

    Args:
        filepaths (list[str]): Paths to unpaired or interleaved FASTQ files.
        chunk_size (int): Number of reads per chunk (read pairs for
            interleaved files).
        parameters (dict): Run parameters (Phred detection settings, and
            the gzip and discard_singles settings copied into paired tasks).
        filetype (str | None): "unpaired" or "interleaved".

    Yields:
        dict: A task dictionary with the task type, filepath(s), the chunk
            as a raw FASTQ blob, and precomputed metadata (Phred offset(s),
            and for paired tasks the gzip and discard_singles settings).

    Raises:
        TypeError: If filetype is neither "unpaired" nor "interleaved".
    """
    for filepath in filepaths:
        start_time = time.monotonic()
        logger.info("%s: started processing.", os.path.basename(filepath))
        if constants.GZIP_DETECTION[filepath]:
            logger.info("%s: gzip format detected.", os.path.basename(filepath))
            if constants.ESTIMATED_ZIP_RATIO.get(filepath) is not None:
                logger.info("%s: estimated gzip compression ratio: %s.", os.path.basename(filepath), constants.ESTIMATED_ZIP_RATIO.get(filepath))
        else:
            logger.info("%s: text format detected.", os.path.basename(filepath))
        logger.info("%s: file size: %s bytes.", os.path.basename(filepath), os.path.getsize(filepath))
        phred_offset = detect_phred_offset(
            filepath=filepath,
            reads_for_phred_offset=parameters["reads_for_phred_offset"],
            phred_offset=parameters["phred_offset"]
        )
        logger.info("%s: Phred offset of %s detected.", os.path.basename(filepath), phred_offset)
        if constants.ESTIMATED_READ_COUNTS.get(filepath) is not None:
            logger.info("%s: estimated bytes per read: %s.", os.path.basename(filepath), constants.ESTIMATED_BYTE_PER_READ[filepath])
            logger.info("%s: (estimated) read count: %s.", os.path.basename(filepath), constants.ESTIMATED_READ_COUNTS[filepath])
        if filetype == "unpaired":
            total_reads = 0
            for blob, n_records in lazy_fastq_blobs(filepath, chunk_size):
                total_reads += n_records
                yield {
                    "type": "unpaired",
                    "filepath": filepath,
                    "chunk": blob,
                    "phred_offset": phred_offset
                }
            elapsed = time.monotonic() - start_time
            if total_reads > 0 and elapsed > 0:
                logger.info(
                    "%s: finished processing %d reads in %.2fs (%.0f reads/sec).",
                    os.path.basename(filepath), total_reads, elapsed, total_reads / elapsed
                )
            else:
                logger.info("%s: finished processing in %.2fs.", os.path.basename(filepath), elapsed)
        elif filetype == "interleaved":
            total_reads = 0
            for blob, n_records in lazy_fastq_blobs(filepath, 2 * chunk_size):
                total_reads += (n_records + 1) // 2
                yield {
                    "type": "paired",
                    "interleaved": True,
                    "file1": filepath,
                    "file2": filepath,
                    "chunk1": blob,
                    "chunk2": None,
                    "phred_offset_1": phred_offset,
                    "phred_offset_2": phred_offset,
                    "gzip_output": parameters["gzip_output"],
                    "gzip_level": parameters["gzip_level"],
                    "discard_singles": parameters["discard_singles"]
                }
            elapsed = time.monotonic() - start_time
            if total_reads > 0 and elapsed > 0:
                logger.info(
                    "%s: finished processing %d read pairs in %.2fs (%.0f pairs/sec).",
                    os.path.basename(filepath), total_reads, elapsed, total_reads / elapsed
                )
            else:
                logger.info("%s: finished processing in %.2fs.", os.path.basename(filepath), elapsed)
        else:
            raise TypeError(f"file {filepath} returned filetype {filetype}, which is not recognized.")

def generate_paired_tasks(files, chunk_size, parameters):
    """
    Lazily split paired-end file pairs into worker tasks.

    For each (R1, R2) pair, logs basic file information, detects the Phred
    offset of each file once, and then yields one task per chunk of
    ``chunk_size`` read pairs, as one raw FASTQ blob per file
    (lazy_fastq_blobs()), so that a single global pool can process chunks
    from many files. Mates are matched by position: record k of the R1
    blob is the mate of record k of the R2 blob, also in
    process_paired_chunk(). Read IDs are not compared.
    Logs the processing time and rate once both files are exhausted.

    Args:
        files (list[tuple[str, str]]): (file1, file2) path tuples for
            paired FASTQ files.
        chunk_size (int): Number of read pairs per chunk.
        parameters (dict): Run parameters (Phred detection settings, and
            the gzip and discard_singles settings copied into each task).

    Yields:
        dict: A "paired" task dictionary with the filepaths, the R1 and R2
            blobs and precomputed metadata.

    Raises:
        ValueError: If one file of a pair runs out of reads before the
            other (different read counts, e.g. a truncated file).
    """
    for pair in files:
        file1, file2 = pair
        start_time = time.monotonic()
        logger.info("%s and %s: started processing.", os.path.basename(file1), os.path.basename(file2))
        for f in (file1, file2):
            if constants.GZIP_DETECTION[f]:
                logger.info("%s: gzip format detected.", os.path.basename(f))
                if constants.ESTIMATED_ZIP_RATIO.get(f) is not None:
                    logger.info("%s: estimated gzip compression ratio: %s.", os.path.basename(f), constants.ESTIMATED_ZIP_RATIO.get(f))
            else:
                logger.info("%s: text format detected.", os.path.basename(f))
            logger.info("%s: file size: %s bytes.", os.path.basename(f), os.path.getsize(f))
        
        phred_offset_1 = detect_phred_offset(filepath=file1, reads_for_phred_offset=parameters["reads_for_phred_offset"], phred_offset=parameters["phred_offset"])
        phred_offset_2 = detect_phred_offset(filepath=file2, reads_for_phred_offset=parameters["reads_for_phred_offset"], phred_offset=parameters["phred_offset"])
        logger.info("%s: Phred offset of %s detected.", os.path.basename(file1), phred_offset_1)
        logger.info("%s: Phred offset of %s detected.", os.path.basename(file2), phred_offset_2)
        logger.info("%s: (estimated) read count: %s.", os.path.basename(file1), constants.ESTIMATED_READ_COUNTS.get(file1, "unknown"))
        logger.info("%s: (estimated) read count: %s.", os.path.basename(file2), constants.ESTIMATED_READ_COUNTS.get(file2, "unknown"))
        blobs_1 = lazy_fastq_blobs(file1, chunk_size)
        blobs_2 = lazy_fastq_blobs(file2, chunk_size)
        total_pairs = 0
        for (blob1, n1), (blob2, n2) in itertools.zip_longest(blobs_1, blobs_2, fillvalue=(None, 0)):
            if n1 != n2:
                raise ValueError(
                    f"Mismatched read counts in paired files {file1} and {file2}. "
                    f"Files must have identical read counts (possible file corruption or truncation)."
                )
            total_pairs += n1
            yield {
                "type": "paired",
                "file1": file1,
                "file2": file2,
                "chunk1": blob1,
                "chunk2": blob2,
                "phred_offset_1": phred_offset_1,
                "phred_offset_2": phred_offset_2,
                "gzip_output": parameters["gzip_output"],
                "gzip_level": parameters["gzip_level"],
                "discard_singles": parameters["discard_singles"]
            }
        elapsed = time.monotonic() - start_time
        if total_pairs > 0 and elapsed > 0:
            logger.info(
                "%s and %s: finished processing %d read pairs in %.2fs (%.0f pairs/sec).",
                os.path.basename(file1), os.path.basename(file2), total_pairs, elapsed, total_pairs / elapsed,
            )
        else:
            logger.info(
                "%s and %s: finished processing in %.2fs.",
                os.path.basename(file1), os.path.basename(file2), elapsed,
            )

def write_summary_and_statistics(summary_results, parameters, output_dir, collector=None):
    """
    Write results_summary.txt to output_dir.

    Without --stats (collector is None), the tab-separated file contains a
    "[Paired Reads]" and an "[Unpaired Reads]" table, each detailing input
    files, output files, and read counts.

    With --stats, it is one tab-separated table with one row per input and
    mate (unpaired inputs: mate "all"; paired and interleaved inputs: "R1"
    and "R2"): input/output files, kept and rejected counts, length and
    quality before and after trimming, one group of columns per trimming
    module and one column per adapter sequence. Cells that do not apply
    are "NA". See readzor.stats.build_summary_grid().

    Args:
        summary_results (dict): File statistics from input_handler().
        parameters (dict): Run parameters.
        output_dir (str): Directory in which to create the file.
        collector (stats.StatsCollector | None): Read statistics from
            input_handler() (--stats), or None.

    Returns:
        None
    """
    summary_path = os.path.join(output_dir, "results_summary.txt")
    if collector is not None:
        S.write_summary(summary_path, summary_results, collector, display_name=basename_file)
        return

    paired_data = []
    unpaired_data = []
    
    for key, stats in summary_results.items():
        if "kept" in stats:  # Unpaired entry
            input_str = stats.get("input", os.path.basename(key))
            output_str = stats.get("output", basename_file(key) + ("_filtered.fastq.gz" if parameters.get("gzip_output") else "_filtered.fastq"))
            unpaired_data.append((input_str, output_str, stats["kept"], stats["rejected"]))
        else:  # Paired or interleaved entry
            input_str = stats.get("input", key)
            output_str = stats.get("output", key)
            paired_data.append((
                input_str,
                output_str,
                stats["kept_pairs"],
                stats["kept_R1_singletons"],
                stats["kept_R2_singletons"],
                stats["rejected_R1"],
                stats["rejected_R2"],
            ))
            
    with open(summary_path, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        
        if paired_data:
            f.write("[Paired Reads]\n")
            writer.writerow(["Input File(s)", "Output File(s)", "Kept_Pairs", "Kept_R1_Singletons", "Kept_R2_Singletons", "Rejected_R1", "Rejected_R2"])
            writer.writerows(paired_data)
            f.write("\n")
            
        if unpaired_data:
            f.write("[Unpaired Reads]\n")
            writer.writerow(["Input File", "Output File", "Kept", "Rejected"])
            writer.writerows(unpaired_data)

def log_parameters(parameters):
    """
    Log all run parameters, sorted by name, one per line and tagged
    [PARAMETER] so they are easy to grep out of the log file. Bytes values
    are decoded, and lists, tuples and sets are joined with commas.

    Args:
        parameters (dict): Dictionary of parameter names and their values.

    Returns:
        None
    """
    
    def _format_value(value):
        """Convert a parameter value into a readable string."""
        if isinstance(value, bytes):
            return value.decode('utf-8', errors='replace')
        if isinstance(value, (list, tuple, set)):
            return ", ".join(_format_value(v) for v in value)
        return str(value)

    logger.info("--- Start of run parameters ---")
    for key, value in sorted(parameters.items()):
        logger.info("[PARAMETER] %s: %s", key, _format_value(value))
    logger.info("---- End of run parameters ----")