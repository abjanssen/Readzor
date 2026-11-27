"""Split one run over several coordinator processes.

A single main process can only feed a limited number of workers: it reads,
frames, pickles, receives and writes every byte itself. When more workers are
requested than one main process can feed, the run is shared by several
coordinator processes, each with its own worker pool.

The inputs (an unpaired file, an interleaved file or an R1/R2 pair) are put in
a shared queue, largest first. Each coordinator takes the next input as soon
as it has finished the previous one, so all coordinators stay busy until the
queue is empty, whatever the file sizes. All coordinators write into the same
results folder and log file; their summaries are merged into one
results_summary.txt.

A single input is never divided, and a run with one input is never split.
"""

import logging
import math
import multiprocessing as mp
import os
import queue as queue_module
import traceback

from readzor import constants
from readzor import stats as S
from readzor.fastq_io import basename_file, common_name_parts, find_paired_files, is_gz_file
from readzor.progress import ProgressTracker, count_reads_estimated

logger = logging.getLogger(__name__)

DEFAULT_MAX_WORKERS_PER_COORDINATOR = 6
GZIP_SIZE_FACTOR = 4  # rough uncompressed/compressed ratio, only used to order the inputs

def resolve_inputs(parameters):
    """
    Detect gzip compression of every input and sort the unspecified files
    into pairs, interleaved and unpaired files, as input_handler() would.
    The result can be passed to input_handler() as explicit groups (with
    unspecified_files=None).

    Returns:
        tuple[list[str], list[tuple[str, str]], list[str]]: unpaired, paired
            and interleaved inputs.
    """
    singles = ((parameters["unspecified_files"] or []) + (parameters["unpaired_files"] or [])
               + (parameters["interleaved_files"] or []))
    for file in singles:
        constants.GZIP_DETECTION[file] = is_gz_file(file)
    for file_1, file_2 in (parameters["paired_files"] or []):
        constants.GZIP_DETECTION[file_1] = is_gz_file(file_1)
        constants.GZIP_DETECTION[file_2] = is_gz_file(file_2)
    auto_paired, auto_interleaved, auto_unpaired = find_paired_files(parameters["unspecified_files"])
    unpaired = auto_unpaired + (parameters["unpaired_files"] or [])
    paired = auto_paired + (parameters["paired_files"] or [])
    interleaved = auto_interleaved + (parameters["interleaved_files"] or [])
    return unpaired, paired, interleaved

def input_weight(files):
    """Estimated uncompressed size of an input, used to process large inputs first."""
    weight = 0
    for file in files:
        size = os.path.getsize(file)
        weight += size * GZIP_SIZE_FACTOR if constants.GZIP_DETECTION[file] else size
    return weight

def make_bundles(unpaired, paired, interleaved):
    """
    Group the inputs into bundles; a bundle is always processed as a whole
    by one coordinator, in one input_handler() call.

    Normally every input is its own bundle. Inputs whose output files would
    get the same name prefix (e.g. equally named files from different
    folders) share a bundle, because output-name deduplication (_2, _3, ...)
    only works within one input_handler() call; otherwise the second would
    fail with FileExistsError.

    Returns:
        list[dict]: Bundles with "unpaired", "paired" and "interleaved" lists
            and a "weight" (estimated uncompressed bytes), largest first.
    """
    bundles = {}
    def bundle_for(prefix):
        return bundles.setdefault(prefix, {"unpaired": [], "paired": [], "interleaved": [], "weight": 0})
    for file in unpaired:
        bundle = bundle_for(basename_file(file))
        bundle["unpaired"].append(file)
        bundle["weight"] += input_weight([file])
    for file in interleaved:
        bundle = bundle_for(basename_file(file))
        bundle["interleaved"].append(file)
        bundle["weight"] += input_weight([file])
    for file_1, file_2 in paired:
        bundle = bundle_for(common_name_parts([os.path.basename(file_1), os.path.basename(file_2)]))
        bundle["paired"].append((file_1, file_2))
        bundle["weight"] += input_weight([file_1, file_2])
    return sorted(bundles.values(), key=lambda b: b["weight"], reverse=True)

def decide_n_coordinators(parameters, n_bundles):
    """
    Number of coordinator processes to use: enough that none has more than
    parameters["max_workers_per_process"] workers, but never more than the
    number of bundles (a coordinator without input would only hold idle
    workers). Returns 1 (no split) for --stdout (two processes cannot share
    one FASTQ stream) and --testrun.
    """
    if parameters["stdout"] or parameters["testrun"]:
        return 1
    max_workers = parameters.get("max_workers_per_process") or DEFAULT_MAX_WORKERS_PER_COORDINATOR
    wanted = math.ceil(parameters["threads"] / max_workers)
    return max(1, min(wanted, n_bundles, parameters["threads"]))

def split_threads(total_threads, n_coordinators):
    """
    Divide all workers as evenly as possible over the coordinators. No
    workers are given up for the coordinators' own main processes: a
    coordinator feeding a few workers only needs a fraction of a core.
    """
    base, extra = divmod(total_threads, n_coordinators)
    return [base + (1 if i < extra else 0) for i in range(n_coordinators)]

def coordinator_main(index, threads, parameters, output_dir, work_queue, result_queue, abort, progress_counter=None):
    """
    Entry point of one coordinator process: take bundles from work_queue
    until it gets None (or abort is set), run input_handler() on each, and
    send (index, status, file_stats, collector) for all of them back.
    """
    # The parent owns the stdin temporary file; a coordinator must not delete
    # it (input_handler's forced-exit path cleans up STDIN_TEMP_FILES).
    constants.STDIN_TEMP_FILES.clear()
    file_stats = {}
    merged = S.StatsCollector() if parameters["stats"] else None
    try:
        from readzor.runner import input_handler
        coordinator_parameters = dict(parameters, threads=threads, show_progress=False)
        logger.info("Coordinator %d (pid %d) started with %d worker(s).", index + 1, os.getpid(), threads)
        n_done = 0
        while not abort.is_set():
            bundle = work_queue.get()
            if bundle is None:
                break
            stats, collector = input_handler(
                unspecified_files=None,
                unpaired_files=bundle["unpaired"],
                paired_files=bundle["paired"],
                interleaved_files=bundle["interleaved"],
                output_dir=output_dir,
                threads=threads,
                chunk_size=coordinator_parameters["chunk_size"],
                show_progress=False,
                stdout=False,
                interleaved_out=coordinator_parameters["interleaved_out"],
                write_rejected=coordinator_parameters["write_rejected"],
                discard_singles=coordinator_parameters["discard_singles"],
                parameters=coordinator_parameters,
                progress_counter=progress_counter,
            )
            file_stats.update(stats)
            if merged is not None and collector is not None:
                merged.inputs.update(collector.inputs)
            n_done += 1
        logger.info("Coordinator %d finished (%d input(s)).", index + 1, n_done)
        result_queue.put((index, "ok", file_stats, merged))
    except KeyboardInterrupt:
        abort.set()
        result_queue.put((index, "interrupted", None, None))
    except BaseException:
        abort.set()
        result_queue.put((index, "error", traceback.format_exc(), None))

def run_split(bundles, n_coordinators, parameters, output_dir):
    """
    Start n_coordinators coordinator processes that share the bundles
    through a queue, wait for all of them and merge their results.

    Returns:
        tuple[dict, stats.StatsCollector | None]: As input_handler().

    Raises:
        KeyboardInterrupt: If any coordinator was interrupted.
        RuntimeError: If any coordinator failed or died.
    """
    shares = split_threads(parameters["threads"], n_coordinators)
    logger.info(
        "Splitting the run over %d coordinator processes (%s workers; at most %d per process), "
        "sharing %d input(s) through a queue.",
        n_coordinators, " + ".join(map(str, shares)),
        parameters.get("max_workers_per_process") or DEFAULT_MAX_WORKERS_PER_COORDINATOR, len(bundles))
    tracker = None
    progress_counter = None
    if parameters["show_progress"]:
        # Estimate before the coordinators start, so they inherit the
        # estimates (used in their log messages) instead of redoing them.
        for bundle in bundles:
            files = bundle["unpaired"] + bundle["interleaved"] + [f for pair in bundle["paired"] for f in pair]
            for file in files:
                constants.ESTIMATED_READ_COUNTS[file] = count_reads_estimated(file)
        progress_counter = mp.Value("q", 0)
        tracker = ProgressTracker(total_reads=sum(constants.ESTIMATED_READ_COUNTS.values()))
        constants.ACTIVE_PROGRESS_TRACKER = tracker
    reported = 0

    def report_progress():
        nonlocal reported
        if tracker is not None:
            done = progress_counter.value
            if done != reported:
                tracker.update(done - reported)
                reported = done

    work_queue = mp.Queue()
    for bundle in bundles:  # largest first, so the last inputs to finish are small
        work_queue.put(bundle)
    for _ in range(n_coordinators):
        work_queue.put(None)  # one stop signal per coordinator
    result_queue = mp.Queue()
    abort = mp.Event()

    processes = []
    for index, threads in enumerate(shares):
        process = mp.Process(
            target=coordinator_main,
            args=(index, threads, parameters, output_dir, work_queue, result_queue, abort, progress_counter),
            name=f"readzor-coordinator-{index + 1}",
        )
        process.start()
        processes.append(process)

    results = {}
    finished = False
    try:
        while len(results) < len(processes):
            report_progress()
            try:
                index, status, payload, collector = result_queue.get(timeout=0.2)
                results[index] = (status, payload, collector)
                continue
            except queue_module.Empty:
                pass
            dead = [i for i, p in enumerate(processes) if i not in results and not p.is_alive()]
            if not dead:
                continue
            # A result may have arrived between the timeout and the liveness check.
            while True:
                try:
                    index, status, payload, collector = result_queue.get_nowait()
                    results[index] = (status, payload, collector)
                except queue_module.Empty:
                    break
            for i in dead:
                if i not in results:
                    exit_code = processes[i].exitcode
                    results[i] = ("interrupted" if exit_code == 130 else "died", exit_code, None)
                    abort.set()  # stop the others from starting new inputs
    except KeyboardInterrupt:
        # SIGINT reaches the coordinators too; give them time to shut down cleanly.
        abort.set()
        for process in processes:
            process.join(timeout=15)
        for process in processes:
            if process.is_alive():
                process.kill()
        raise
    finally:
        if tracker is not None:
            report_progress()
            finished = all(r[0] == "ok" for r in results.values()) and len(results) == len(processes)
            tracker.close(interrupted=not finished)
            constants.ACTIVE_PROGRESS_TRACKER = None
    for process in processes:
        process.join()

    file_stats = {}
    merged = S.StatsCollector() if parameters["stats"] else None
    failures = []
    interrupted = False
    for index in sorted(results):
        status, payload, collector = results[index]
        if status == "ok":
            file_stats.update(payload)
            if merged is not None and collector is not None:
                merged.inputs.update(collector.inputs)
        elif status == "interrupted":
            interrupted = True
        elif status == "error":
            failures.append(f"Coordinator {index + 1} failed:\n{payload}")
        else:
            failures.append(f"Coordinator {index + 1} exited unexpectedly (exit code {payload}).")
    if interrupted:
        raise KeyboardInterrupt
    if failures:
        raise RuntimeError("\n".join(failures))
    # Same order as a normal run: unpaired, then interleaved, then paired inputs.
    type_order = {"unpaired": 0, "interleaved": 1, "paired": 2}
    file_stats = dict(sorted(file_stats.items(), key=lambda item: (type_order.get(item[1].get("type"), 3), item[0])))
    return file_stats, merged
