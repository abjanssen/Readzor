#!/usr/bin/env python3

##### Import packages #####
import argparse
from collections import defaultdict
import datetime
import itertools
import logging
import multiprocessing as mp
import os
import random
import re
import shlex
import shutil
import signal
import sys
import tempfile
import time
import threading

from numpy.lib.stride_tricks import sliding_window_view
from isal import igzip as gzip
import numpy as np

##### Definition of constant values #####
WORKER_PARAMETERS = None
ESTIMATED_ZIP_RATIO = {}
ESTIMATED_READ_COUNTS = {}
STDIN_TEMP_FILES = []
ESTIMATED_BYTE_PER_READ = {}
GZIP_DETECTION = {}
VERSION = "0.2.1"
PHRED_ALLOWED = bytes(range(33, 127))
DEFAULT_ADAPTERS = [
    ["TruSeq3", [
        ["TruSeq3_R1_short", "AGATCGGAAGAGCACA"],  # first 16 of full seq
        ["TruSeq3_R2_short", "AGATCGGAAGAGCGTC"],  # first 16 of full seq
    ]],
    ["TruSeq2", [
        ["TruSeq2", "AGATCGGAAGAGCGGTTCAG"],
    ]],
    ["TruSeq_small_RNA", [
        ["TruSeq_small_RNA", "TGGAATTCTCGGGTGCCAAGG"],
    ]],
    ["Nextera", [
        ["Nextera", "CTGTCTCTTATACACATCT"],
    ]],
    ["Illumina_RNA", [
        ["Illumina_RNA", "ACTGTCTCTTATACACATCT"],
    ]],
]
FULL_AUTO_PRESERVED_DESTS = {"input_files", "input_paired", "input_unpaired", "full_auto"}
FIELD_SEP = b"\x1f"
FULL_AUTO_OVERRIDES = {
    "endqual_filter_flag": True,
    "adapter_filter_flag": True,
    "n_filter": True,
    "min_length_output_perc": 33,
    "gzip": True,
    "progress": True
}
NUCL_ATCG = b"ATCG"
NUCL_ATCGN = b"ATCGN"
ACTIVE_PROGRESS_TRACKER = None
PHRED64_TO_33 = bytes.maketrans(
    bytes(range(59, 127)),
    bytes(max(33, b - 31) for b in range(59, 127))
)
PHRED33_TO_64 = bytes.maketrans(
    bytes(range(33, 127)),
    bytes(min(126, b + 31) for b in range(33, 127))
)

##### Logging #####
logger = logging.getLogger("readzor")
file_only_logger = logging.getLogger("readzor.file_only")

class ProgressAwareStreamHandler(logging.StreamHandler):
    """
    A StreamHandler that clears any in-progress progress-bar line before
    emitting a log record, so the two don't get interleaved on the same
    terminal line.
    """
    def emit(self, record):
        if ACTIVE_PROGRESS_TRACKER is not None:
            ACTIVE_PROGRESS_TRACKER.clear_line()
        super().emit(record)

def setup_logging(output_dir = None, verbose = False, parameters = None):
    """
    Configure Readzor's logger.
    Called from parse_args() twice — once with no arguments, before parsing,
    to attach a console handler so any parser.error()/early messages are
    visible; and again immediately after args are parsed, to apply the
    user's --verbose choice to that same console handler. Called a third
    time from main() once the timestamped output directory exists, to attach
    a file handler so the full run is captured on disk
    (Readzor_results_.../readzor.log) regardless of console verbosity.
    Args:
        output_dir (str | None): Timestamped results directory to write
            readzor.log into. If None, only the console handler is
            (re)configured.
        verbose (bool): If True, console output includes DEBUG-level
            messages. The log file always captures DEBUG and above.
        parameters (dict | None): Run parameters. Required (and only
            consulted) when output_dir is given, to log whether full-auto
            and/or test-run mode is active.
    """
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter(
        "%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )
    console_level = logging.DEBUG if verbose else (logging.CRITICAL + 1)
    if not logger.handlers:
        console = ProgressAwareStreamHandler(sys.stderr)
        console.setLevel(console_level)
        console.setFormatter(formatter)
        console.name = "console"
        logger.addHandler(console)
    else:
        for handler in logger.handlers:
            if getattr(handler, "name", None) == "console":
                handler.setLevel(console_level)
    if output_dir is not None and not any(getattr(h, "name", None) == "file" for h in logger.handlers):
        log_path = os.path.join(output_dir, "Readzor_log.txt")
        file_handler = logging.FileHandler(log_path)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        file_handler.name = "file"
        logger.addHandler(file_handler)
        logger.info("Readzor %s starting.", VERSION)
        used_command = " ".join(map(shlex.quote, [sys.executable] + sys.argv))
        logger.info("Command used: %s", used_command)
        if parameters["full_auto"]:
            logger.warning("--full-auto/-GO specified; ignoring all other input parameters (except input file parameters).")
        if parameters["testrun"]:
            logger.info("Test run output directory temporarily created at %s", output_dir)
            logger.info("Test run log file temporarily initialized at %s", log_path)
        else:
            logger.info("Output directory created at %s", output_dir)
            logger.info("Log file initialized at %s", log_path)
        file_only_logger.setLevel(logging.DEBUG)
        file_only_logger.propagate = False
        file_only_logger.addHandler(file_handler)
        
##### Progress tracker #####
def estimate_bytes_per_read(filepath, sample_size=10):
    """
    Estimates the average on-disk (uncompressed) bytes of a file's reads 
    by sampling for header, sequence, plus-line, and
    quality lengths and adding back the newline stripped by `lazy_fastq`.

    Args:
        filepath (str): Path to the FASTQ file (.fastq, .fq, or gzip-compressed).
        sample_size (int): Maximum number of records to sample. Defaults to 10.

    Returns:
        int: Estimated number of bytes one record occupies on disk,
            including line-ending characters.

    Raises:
        ValueError: If the file contains no readable FASTQ records.
    """
    total_bytes = 0
    n_records = 0
    for record in lazy_fastq(filepath):
        if n_records >= sample_size:
            break
        header, sequence, plus, quality = record.split(FIELD_SEP)
        header_bytes = len(header) + 1
        seq_bytes = len(sequence) + 1
        plus_bytes = len(plus) + 1
        qual_bytes = len(quality) + 1
        total_bytes += header_bytes + seq_bytes + plus_bytes + qual_bytes
        n_records += 1
    if n_records == 0:
        raise ValueError(f"No FASTQ records found in '{filepath}'; cannot estimate bytes per read.")
    estimate = int(total_bytes / n_records)
    ESTIMATED_BYTE_PER_READ[filepath] = estimate
    return estimate

def estimate_gzip_ratio(filepath, sample_bytes=50 * 1024 * 1024):
    """
    Estimate a gzip file's compression ratio (uncompressed / compressed) by
    decompressing a leading sample, rather than the whole file.

    Args:
        filepath (str): Path to the gzip-compressed file.
        sample_bytes (int): Target number of uncompressed bytes to sample
            before stopping. Defaults to 50 MiB.

    Returns:
        float | None: The uncompressed/compressed ratio, or None if the
        sample is too small to give a stable ratio (e.g. the whole file is
        tiny) -- caller should fall back to a fixed default ratio in that case.
    """
    compressed_read = 0
    uncompressed_read = 0
    with open(filepath, 'rb') as raw:
        decompressor = gzip.GzipFile(fileobj=raw)
        while uncompressed_read < sample_bytes:
            chunk = decompressor.read(1024 * 1024)
            if not chunk:
                break
            uncompressed_read += len(chunk)
        compressed_read = raw.tell()
    if compressed_read == 0 or uncompressed_read < 1024 * 1024:
        return None
    return uncompressed_read / compressed_read

def chunk_size_setter(chunk_size):
    """
    Resolves the chunk size to use for parallel processing.
    If a chunk size is explicitly given, it's returned unchanged.
    Otherwise, detects whether a supported job scheduler's tools are
    available on the system (SLURM, PBS/Torque, SGE, or LSF) and picks
    a larger default chunk size for scheduler-managed (typically
    higher-resource) systems, or a smaller default otherwise.
    Args:
        chunk_size (int | None): User-specified chunk size. Defaults to
            None (auto-detect).
    Returns:
        int: The resolved chunk size — 20000 on scheduler-managed systems,
            1000 otherwise, unless overridden.
    """
    if chunk_size is not None:
        return chunk_size
    scheduler_tools = ('sinfo','sbatch','squeue','qsub','qstat','bsub','bjobs')
    if any(shutil.which(tool) is not None for tool in scheduler_tools):
        return 20000
    else:
        return 1000

def count_reads_estimated(filepath, default_gzip_ratio=4):
    """
    Estimate the number of reads in a FASTQ file from file size and a
    per-read byte estimate, without a full parse.
 
    For gzip files, estimates the uncompressed size via a sampled
    compression ratio (falls back to `default_gzip_ratio` if the file is
    too small to sample reliably).

    Args:
        filepath (str): Path to the FASTQ file (.fastq, .fq, or gzip-compressed).
        default_gzip_ratio (float): Fallback compression ratio to use when
            the file is too small to sample a stable ratio. Defaults to 4.

    Returns:
        int: Estimated number of reads in the file, always >= 1.
    """
    file_size = os.path.getsize(filepath)
    bytes_per_read = estimate_bytes_per_read(filepath)
    if GZIP_DETECTION[filepath]:
        ratio = estimate_gzip_ratio(filepath)
        if ratio is None:
            ratio = default_gzip_ratio
        ESTIMATED_ZIP_RATIO[filepath] = ratio
        estimated_uncompressed_size = file_size * ratio
    else:
        estimated_uncompressed_size = file_size
    return max(1, round(estimated_uncompressed_size / bytes_per_read))

class ProgressTracker:
    """
    Tracks number of analysed reads against estimated total and renders a single-line
    progress bar to stderr. Dynamically scales bar width to fit screen size.
    """
    def __init__(self, total_reads, bar_width=50, min_interval=0.2):
        self.total = max(total_reads, 1)
        self.done = 0
        self.bar_width = bar_width
        self.min_interval = min_interval
        self._last_render = 0.0
        self._start_time = time.time()

    @staticmethod
    def _format_duration(seconds_val, concise=False):
        total_seconds = int(round(seconds_val))
        hours, remainder = divmod(total_seconds, 3600)
        minutes, seconds = divmod(remainder, 60)
        if concise:
            if hours > 0:
                return f"{hours}h {minutes:02d}m {seconds:02d}s"
            if minutes > 0:
                return f"{minutes}m {seconds:02d}s"
            return f"{seconds}s"
        def plural(value, unit):
            return f"{value} {unit}{'s' if value != 1 else ''}"
        parts = []
        if hours > 0:
            parts.append(plural(hours, 'hour'))
        if minutes > 0:
            parts.append(plural(minutes, 'minute'))
        if seconds > 0 or total_seconds < 60:
            parts.append(plural(seconds, 'second'))
        return ", ".join(parts)

    def update(self, n):
        self.done += n
        now = time.time()
        if now - self._last_render >= self.min_interval:
            self._render()
            self._last_render = now

    def _render(self):
        if not sys.stderr.isatty():
            return
        term_width = shutil.get_terminal_size(fallback=(80, 24)).columns
        frac = min(self.done / self.total, 0.999)
        elapsed = time.time() - self._start_time
        rate = self.done / elapsed if elapsed > 0 else 0
        eta = (self.total - self.done) / rate if rate > 0 else float('inf')
        eta_str = self._format_duration(eta, concise=True)
        stats = (
            f" {frac*100:5.1f}% "
            f"({self.done:,}/{self.total:,} reads). "
            f"Rate: {rate:,.0f} reads/s. Estimated time remaining: {eta_str}"
        )
        max_bar_len = term_width - len(stats) - 3
        if max_bar_len >= 5:
            effective_bar_width = min(self.bar_width, max_bar_len)
            filled = int(effective_bar_width * frac)
            progressbar = "#" * filled + "-" * (effective_bar_width - filled)
            line = f"\r[{progressbar}]{stats}"
        else:
            line = f"\r{stats.strip()}"
        line = line[: term_width - 1].ljust(term_width - 1)
        sys.stderr.write(line)
        sys.stderr.flush()

    def clear_line(self):
        """
        Blank out the current progress-bar line so other output (e.g. log
        messages) can be written to stderr without visually clashing with
        the bar. The bar redraws itself on the next call to `update()`.
        """
        if not sys.stderr.isatty():
            return
        term_width = shutil.get_terminal_size(fallback=(80, 24)).columns
        sys.stderr.write("\r" + " " * (term_width - 1) + "\r")
        sys.stderr.flush()

    def close(self):
        elapsed = time.time() - self._start_time
        rate = self.done / elapsed if elapsed > 0 else 0
        time_str = self._format_duration(elapsed, concise=False)
        stats = (
            f" 100% ({self.done:,} reads analyzed). "
            f"Average rate: {rate:,.0f} reads/s. Total time: {time_str}."
        )
        if not sys.stderr.isatty():
            sys.stderr.write(stats.strip() + "\n")
            sys.stderr.flush()
            return
        term_width = shutil.get_terminal_size(fallback=(80, 24)).columns
        max_bar_len = term_width - len(stats) - 3
        if max_bar_len >= 5:
            effective_bar_width = min(self.bar_width, max_bar_len)
            progressbar = "#" * effective_bar_width
            line = f"\r[{progressbar}]{stats}"
        else:
            line = f"\r{stats.strip()}"
        line = line[: term_width - 1].ljust(term_width - 1)
        sys.stderr.write(line + "\n")
        sys.stderr.flush()

##### Helper functions #####

def resolve_stdin_input():
    """
    Auto-detects piped stdin data and spools it to a temporary file.

    Returns:
        str: A real filesystem path safe to pass through the rest of the pipeline.

    Raises:
        SystemExit: If no data is being piped into stdin.
    """
    if sys.stdin.isatty():
        raise SystemExit("Error: No input file specified and no piped data detected on stdin.")

    fd, tmp_path = tempfile.mkstemp(suffix=".fastq", prefix="readzor_stdin_")
    with os.fdopen(fd, "wb") as tmp_file:
        shutil.copyfileobj(sys.stdin.buffer, tmp_file, length=10 * 1024 * 1024)
    STDIN_TEMP_FILES.append(tmp_path)
    return tmp_path

def cleanup_stdin_temp_files():
    """Remove any temp files created by resolve_stdin_input(), ignoring ones already gone."""
    for tmp_path in STDIN_TEMP_FILES:
        try:
            os.remove(tmp_path)
        except OSError:
            pass
    STDIN_TEMP_FILES.clear()

def group_paired_input_into_pairs(files, parser):
    """
    Group a flat list of paired-end FASTQ files into (R1, R2) tuples.

    Assumes files are supplied in consecutive R1/R2 order, e.g.
    [sample1_R1, sample1_R2, sample2_R1, sample2_R2] -> [(sample1_R1, sample1_R2), (sample2_R1, sample2_R2)].
    
    Args:
        files (list[str] | None): Flat list of input file paths, or None/empty
            if no paired input was provided.
        parser (argparse.ArgumentParser): Parser used to report a usage error
            (via parser.error) if the file count is invalid.
    
    Returns:
        list[tuple[str, str]] | None: List of (R1, R2) file path tuples, or
        None if `files` is empty/None.
    
    Raises:
        SystemExit: Raised indirectly via parser.error() if `files` contains
            an odd number of entries, since paired input requires an even count.
    """
    if not files:
        return None
    if len(files) % 2 != 0:
        parser.error(f"argument --paired: expected an even number of files (R1 R2 pairs), got {len(files)}.")
    pairs = [tuple(files[i:i+2]) for i in range(0, len(files), 2)]
    return pairs

def create_folder_structure(output_dir):
    """
    Create a timestamped results folder inside a specified output directory.
    
    Generates a new subfolder named "Readzor_results_<YYYY-MM-DD_HH-MM-SS>"
    based on the current date and time. This ensures each run's outputs are
    isolated and prevents accidental overwrites of results from previous runs.

    Args:
        output_dir (str): Path to the parent directory where the timestamped
            results folder will be created. Created if it does not exist.

    Returns:
        str: The absolute path to the newly created timestamped results folder.
    
    Raises:
        OSError: If the folder cannot be created due to permission issues or
            invalid path.
    """
    created_output_dir = os.path.join(output_dir, "Readzor_results_" + datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S"))
    os.makedirs(created_output_dir, exist_ok=True)
    return created_output_dir

def worker_determination(threads=None):
    """
    Determine the number of worker processes for parallel task execution.

    The worker count is calculated in the following order of precedence:
    1. Explicit Request: If `threads` is a valid integer, it is used.
    2. Scheduler Environment: If a supported job-scheduler CPU variable
       (SLURM, PBS, SGE, LSF, or OMP_NUM_THREADS) is set and valid, it is used.
    3. Local Default: Uses (available CPUs - 1).

    All return values are clamped to a minimum of 1 worker.

    Args:
        threads (int | None): User-requested number of base threads/CPUs.
            Ignored if not a valid int (e.g., None, float, bool). Defaults to None.

    Returns:
        int: The number of worker processes to spawn, always >= 1.
    """
    if isinstance(threads, int):
        return max(1, threads)
    for var in ("SLURM_CPUS_PER_TASK", "PBS_NCPUS", "NCPUS", "NSLOTS", "LSB_DJOB_NUMPROC", "OMP_NUM_THREADS"):
        val = os.environ.get(var)
        if val:
            try:
                return max(1, int(val))
            except ValueError:
                pass
    return max(1, (os.cpu_count() or 1) - 1)

def common_name_parts(filenames):
    """
    Extract the common underscore-separated stem shared across FASTQ filenames.
    
    Strips file extensions (.fastq, .fq, optionally .gz/.gzip) and identifies
    tokens that are identical across all input filenames at each position.
    Stops at the first position where tokens differ, allowing paired-end files
    (e.g., sample1_R1.fastq.gz and sample1_R2.fastq.gz) to be reduced to their
    shared prefix (sample1).
    
    If no common tokens are found, returns the first filename's stem as a fallback
    to ensure a non-empty identifier is always available for downstream processing.

    Args:
        filenames (list[str]): List of FASTQ filenames (basenames, not paths).
            Can be single or multiple files. Extensions are case-insensitive.
            Must not be empty.

    Returns:
        str: The common name stem, with tokens rejoined by underscores. If no
        common tokens exist, returns the first filename's stem (without extension).
        Never returns an empty string.

    Examples:
        - ["sample1_R1.fastq.gz", "sample1_R2.fastq.gz"] -> "sample1"
        - ["file.fq.gz"] -> "file"
        - ["a_b_c.fastq", "x_y_z.fastq"] -> "a_b_c"  (first file's stem as fallback)
    """
    if not filenames:
        return "unknown"
    stems = [re.sub(r'\.(fastq|fq)(\.gzip|\.gz)?$', '', file, flags=re.IGNORECASE) for file in filenames]
    token_lists = [stem.split('_') for stem in stems]
    min_len = min(len(tokens) for tokens in token_lists)
    common_tokens = []
    for i in range(min_len):
        values_at_position = {tokens[i].lower() for tokens in token_lists}
        if len(values_at_position) != 1:
            break
        common_tokens.append(token_lists[0][i])
    output = '_'.join(common_tokens) if common_tokens else stems[0]
    return output

def basename_file(filepath):
    """
    Extract the sample name from a FASTQ filepath by removing directory and extension.
    
    Strips the directory path and FASTQ file extension (.fastq, .fq, optionally
    .gz/.gzip) to isolate the sample identifier. Extension matching is case-insensitive.

    Args:
        filepath (str): Path to the FASTQ file (absolute or relative).

    Returns:
        str: The filename stem without directory or FASTQ extension.

    Examples:
        - "/data/reads/sample1_R1.fastq.gz" -> "sample1_R1"
        - "sample2_R2.FQ.GZ" -> "sample2_R2"
    """
    filename = os.path.basename(filepath)
    reduced_filename = re.sub(r'\.(fastq|fq)(\.gzip|\.gz)?$', '', filename, flags=re.IGNORECASE)
    return reduced_filename

def is_gz_file(filepath):
    """
    Determine if a file is gzip-compressed by inspecting its magic bytes.
    
    Reads the first two bytes of the file and compares them to the gzip magic
    number (0x1f 0x8b), providing reliable detection independent of filename
    or extension. Gracefully handles missing or inaccessible files.

    Args:
        filepath (str): Path to the file to check (absolute or relative).

    Returns:
        bool: True if the file is gzip-compressed; False if not compressed,
        the file does not exist, or cannot be read due to permissions.
    """
    try:
        with open(filepath, 'rb') as file:
            return file.read(2) == b'\x1f\x8b'
    except (IOError, OSError):
        return False

def lazy_fastq(filepath):
    """
    Lazily yield FASTQ records one at a time without loading the entire file.

    Detects gzip compression via magic bytes (independent of file extension).
    Gzip files are read via isal in large chunks with manual line-stitching
    across chunk boundaries (fastest path for isal's decompression throughput).
    Plain-text files are read via Python's line iterator with a 10 MB read
    buffer (fastest path for uncompressed I/O). Yields each 4-line FASTQ
    record as (header, sequence, plus, quality) tuples with line endings
    stripped.

    Args:
        filepath (str): Path to the FASTQ file (.fastq, .fq, or gzip-compressed).
    Yields:
        bytes: One record's (header, sequence, plus, quality) lines joined
            by FIELD_SEP into a single bytes object, line endings stripped.
    """
    buffer_size = 10 * 1024 * 1024
    if GZIP_DETECTION[filepath]:
        fastq_file = gzip.open(filepath, 'rb')
        with fastq_file:
            leftover = b""
            while True:
                chunk = fastq_file.read(buffer_size)
                if not chunk:
                    if leftover:
                        lines = leftover.rstrip(b"\r").split(b"\n")
                        usable = len(lines) - (len(lines) % 4)
                        it = iter(lines[:usable])
                        for group in zip(it, it, it, it):
                            yield FIELD_SEP.join(group)
                    break
                data = leftover + chunk
                last_newline = data.rfind(b"\n")
                if last_newline == -1:
                    leftover = data
                    continue
                lines = data[:last_newline].split(b"\n")
                usable = len(lines) - (len(lines) % 4)
                it = iter(lines[:usable])
                for group in zip(it, it, it, it):
                    yield FIELD_SEP.join(group)
                leftover_lines = lines[usable:]
                if leftover_lines:
                    leftover = b"\n".join(leftover_lines) + b"\n" + data[last_newline + 1:]
                else:
                    leftover = data[last_newline + 1:]
    else:
        fastq_file = open(filepath, 'rb', buffering=buffer_size)
        with fastq_file:
            lines = iter(fastq_file)
            for header in lines:
                header = header.strip()
                try:
                    sequence = next(lines).strip()
                    plus = next(lines).strip()
                    quality = next(lines).strip()
                except StopIteration:
                    break
                yield FIELD_SEP.join((header, sequence, plus, quality))

def find_paired_files(filepaths):
    """
    Groups a list of FASTQ filepaths into paired-end pairs and leftover unpaired files.
    For each file, reads the first line (its header) to determine a
    "base ID" and read number (e.g., R1/R2) via helper functions,
    without parsing the whole file. Files sharing a base ID with valid
    mate designations are paired together (ensuring correct R1, R2 order);
    leftover files, unmatchable files, or base IDs with more than two files
    are returned as unpaired.
    Args:
        filepaths (list[str]): Paths to FASTQ files to inspect and pair.
    Returns:
        tuple[list[tuple[str, str]], list[str], list[str]]:
            - A list of (file_1, file_2) tuples representing matched pairs (ordered R1, R2).
            - A list of filepaths whose first two records share a base ID
              (i.e. are interleaved R1/R2 reads within a single file).
            - A list of remaining filepaths that could not be matched with
              a valid pair and aren't interleaved.
    """
    base_ids = {}
    interleaved_flag = {}
    if filepaths is None:
        filepaths = []
    logger.info("Inspecting %s candidate file(s) for file pairing.", len(filepaths))

    for filepath in filepaths:
        try:
            if GZIP_DETECTION[filepath]:
                raw = open(filepath, 'rb', buffering=1024)
                try:
                    file = gzip.GzipFile(fileobj=raw)
                    with file:
                        lines = [file.readline() for _ in range(8)]
                finally:
                    raw.close()
            else:
                with open(filepath, 'rb', buffering=1024) as file:
                    lines = [file.readline() for _ in range(8)]
        except (IOError, OSError) as e:
            logger.warning("Could not read '%s', skipping: %s", filepath, e)
            continue
        headers = [lines[i] for i in (0, 4) if lines[i]]
        if not headers:
            logger.warning("File '%s' has an empty or missing header, skipping.", filepath)
            continue

        header_1 = headers[0].strip().lstrip(b'@')
        base_id, read_num = read_info_from_header(header_1)
        base_ids[filepath] = (base_id, read_num)
        if len(headers) == 2:
            header_2 = headers[1].strip().lstrip(b'@')
            base_id_2, read_num_2 = read_info_from_header(header_2)
            interleaved_flag[filepath] = (
                base_id is not None
                and base_id == base_id_2
            )
        else:
            interleaved_flag[filepath] = False

    pairs_dict = defaultdict(list)
    for filepath, (base_id, read_num) in base_ids.items():
        pairs_dict[base_id].append((filepath, read_num))

    pairs = []
    dropped = set()
    for base_id, file_list in pairs_dict.items():
        if len(file_list) == 2:
            (file_1, read_1), (file_2, read_2) = file_list
            if read_1 == 1 and read_2 == 2:
                logger.info("Paired '%s' (R1) with '%s' (R2).", os.path.basename(file_1), os.path.basename(file_2))
                pairs.append((file_1, file_2))
                dropped.add(file_1)
                dropped.add(file_2)
            elif read_1 == 2 and read_2 == 1:
                logger.info("Paired '%s' (R1) with '%s' (R2).", os.path.basename(file_1), os.path.basename(file_2))
                pairs.append((file_2, file_1))
                dropped.add(file_1)
                dropped.add(file_2)
            elif read_1 == read_2:
                logger.warning(
                    "Files '%s' and '%s' share base ID '%s', but also the same read number (%s). "
                    "These files might be duplicate copies. Treating as unpaired files instead.",
                    file_1, file_2, base_id.decode('utf-8', errors='replace'), read_1,
                )
            else:
                logger.warning(
                    "Files '%s' and '%s' share base ID '%s' but do not form a "
                    "valid R1/R2 pair (read numbers: %s, %s). Treating as unpaired files instead.",
                    file_1, file_2, base_id.decode('utf-8', errors='replace'), read_1, read_2
                )
        elif len(file_list) > 2:
            logger.warning(
                "Found %s files matching base ID '%s' (expected max 2 for "
                "paired-end data). Treating as unpaired files instead.",
                len(file_list), base_id.decode('utf-8', errors='replace')
            )

    remaining = [f for f in base_ids if f not in dropped]
    interleaved = [f for f in remaining if interleaved_flag.get(f)]
    unpaired = [f for f in remaining if not interleaved_flag.get(f)]

    logger.info(
        "Matched %s pair(s), %s interleaved file(s), %s file(s) left unpaired.",
        len(pairs), len(interleaved), len(unpaired),
    )
    if pairs:
        logger.info("Paired files:\n%s", "\n".join(f"  {os.path.basename(f1)} + {os.path.basename(f2)}" for f1, f2 in pairs))
    if interleaved:
        logger.info("Interleaved files:\n%s", "\n".join(f"  {os.path.basename(f)}" for f in interleaved))
    if unpaired:
        logger.info("Unpaired files:\n%s", "\n".join(f"  {os.path.basename(f)}" for f in unpaired))

    return pairs, interleaved, unpaired

def read_info_from_header(header):
    """
    Extracts the base read identifier and read number from a FASTQ header.

    Handles two common Illumina header conventions:
      - Modern Illumina: "@<id> 1:N:0:..." or "@<id> 2:N:0:..."
      - Legacy Illumina, and MGI: "@<id>/1" or "@<id>/2"
      
    If neither convention matches, the entire header is used as the base ID 
    and the read number is set to None.

    Args:
        header (str): The FASTQ header line, with the leading '@' and
            trailing newline already stripped.

    Returns:
        tuple[str, int | None]:
            - base_id: The read identifier shared by R1/R2 mates.
            - read_num: 1, 2, or None if the read number cannot be determined.
    """
    read_num = None
    header = header.strip()
    if b' ' in header:
        parts = header.split(b' ')
        base_id = parts[0]
        m = re.match(rb'([12]):?', parts[1])
        if m:
            read_num = int(m.group(1))
    else:
        m = re.search(rb'/([12])$', header)
        if m:
            read_num = int(m.group(1))
            base_id = re.sub(rb'/[12]$', b'', header)
        else:
            base_id = header
    return base_id, read_num

def detect_phred_offset(filepath, reads_for_phred_offset, phred_offset):
    """
    Auto-detects the Phred quality score encoding offset of a FASTQ file.

    If a pre-determined `phred_offset` is provided, it returns it immediately. 
    Otherwise, it samples up to `reads_for_phred_offset` reads and inspects the 
    ASCII range of their quality strings to distinguish Phred+33 (Sanger/modern Illumina) 
    from Phred+64 (older Illumina) encoding.

    Args:
        filepath (str): Path to the FASTQ file to inspect.
        reads_for_phred_offset (int): Maximum number of reads to sample.
        phred_offset (int | None): User-specified Phred offset, if already known.

    Returns:
        int: 33 or 64, the detected or provided Phred offset.

    Raises:
        ValueError: If the FASTQ file cannot be read, or if the observed 
            ASCII range is ambiguous and doesn't clearly match either encoding.
    """
    if phred_offset is not None:
        return phred_offset
    min_ascii = 127
    max_ascii = 0
    count = 0
    reader = lazy_fastq(filepath)
    try:
        for record in reader:
            if count >= reads_for_phred_offset:
                break
            _, _, _, quality = record.split(FIELD_SEP)
            q_bytes = np.frombuffer(quality, dtype=np.uint8)
            min_ascii = min(min_ascii, q_bytes.min())
            max_ascii = max(max_ascii, q_bytes.max())
            count += 1
    except (FileNotFoundError, IOError) as error:
        raise ValueError(f"Cannot read FASTQ file '{filepath}': {error}") from error
    finally:
        reader.close()
    if min_ascii < 64:
        return 33
    if max_ascii <= 104:
        return 64
    else:
        min_char = chr(min_ascii) if 32 <= min_ascii <= 126 else '?'
        max_char = chr(max_ascii) if 32 <= max_ascii <= 126 else '?'
        raise ValueError(
            f"Ambiguous Phred encoding detected in {filepath} (ASCII range {min_ascii}-{max_ascii} ['{min_char}' - '{max_char}'])."
            f"Please specify the Phred offset (33/64) manually using the --phred-offset option."
        )

def validate_fastq(header, sequence, plus, quality, n_filter, min_length_input, max_length_input):
    """
    Validates that a single FASTQ record is well-formed.

    Checks that the header starts with '@' and the plus-line starts with
    '+', that the sequence and quality lines are the same length and
    within the given length bounds, that the sequence contains only
    allowed nucleotide bytes (ACTG, or ACTGN if `n_filter` is off), and
    that the quality string contains only printable Phred-range bytes.

    Args:
        header (bytes): The header line.
        sequence (bytes): The nucleotide sequence line.
        plus (bytes): The plus-separator line.
        quality (bytes): The quality-score line.
        n_filter (bool): If True, reject reads containing any 'N' base;
            if False, 'N' is allowed alongside A/C/T/G.
        min_length_input (int | None): Minimum acceptable sequence length,
            or None to skip this check.
        max_length_input (int | None): Maximum acceptable sequence length,
            or None to skip this check.

    Returns:
        bool: True if the record is well-formed and passes all checks,
            False otherwise.
    """
    if len(header) <= 1 or header[0] != 64:
        return False
    if not plus or plus[0] != 43:
        return False
    seq_len = len(sequence)
    if seq_len != len(quality):
        return False
    if (min_length_input is not None and seq_len < min_length_input) or (max_length_input is not None and seq_len > max_length_input):
        return False
    if n_filter:
        if sequence.translate(None, NUCL_ATCG):
            return False
    else:
        if sequence.translate(None, NUCL_ATCGN):
            return False
    if quality.translate(None, PHRED_ALLOWED):
        return False
    return True

def load_adapters_from_fasta(fasta_file):
    """
    Parses adapter sequences from a FASTA file.

    Reads a FASTA-formatted file containing adapter sequences, extracting
    each adapter's name and sequence. Handles multi-line sequences (where
    a single sequence may span multiple lines in the FASTA file) by
    concatenating all lines between headers.

    Args:
        fasta_file (str): Path to the FASTA file containing adapter sequences.
            Each adapter entry must have a header line starting with ">" followed
            by one or more sequence lines. Sequence lines may be split across
            multiple lines.

    Returns:
        list[tuple[str, str]]: A list of (name, sequence) tuples, where:
            - name (str): The adapter name, extracted from the header line
              (with leading ">" and surrounding whitespace removed).
            - sequence (str): The full nucleotide sequence, with all lines
              between headers concatenated into a single string.

    Raises:
        ValueError: If a header is encountered without any following sequence lines.
    """
    result = []
    with open(fasta_file) as fastafile:
        header = fastafile.readline().rstrip()
        while header:
            if not header.strip():
                header = fastafile.readline().rstrip()
                continue
            if not header.startswith(">"):
                raise ValueError("Malformed adapter Fasta file detected.")
            output_sequence = []
            sequence = fastafile.readline().rstrip()
            while sequence and not sequence.startswith(">"):
                output_sequence.append(sequence.upper())
                sequence = fastafile.readline().rstrip()
            if len(output_sequence) == 0:
                raise ValueError(f"No sequence detected for header {header}")
            name = header.lstrip(">").strip()
            joined_sequence = "".join(output_sequence)
            result.append((name, joined_sequence))
            header = sequence
        return result

def qual_to_array(quality_list, phred_offset):
    """
    Converts a list of Phred quality strings into a 2D numeric numpy array.

    Concatenates all quality strings, reinterprets the raw bytes as an
    array of Phred-shifted integer quality scores, subtracts the Phred offset,
    and reshapes into a (n_reads, read_length) matrix for vectorized downstream processing.

    Args:
        quality_list (list[str]): Quality strings, all of equal length.
        phred_offset (int): The Phred encoding offset (33 or 64) to subtract.

    Returns:
        numpy.ndarray: Signed 8-bit integer array of shape
            (len(quality_list), read_length) with true quality scores.
    """
    quality_arr_lengths = np.array([len(q) for q in quality_list])
    n_reads = len(quality_list)
    max_len = quality_arr_lengths.max()
    min_len = quality_arr_lengths.min()
    if max_len == min_len:
        joined = b''.join(quality_list)
        array = np.frombuffer(joined, dtype=np.int8).reshape(n_reads, max_len)
        array = array - phred_offset
        chunk_padding_bool = False
        row_tilde_count = 0
        padding_mask_bool = None
        return array, chunk_padding_bool, row_tilde_count, padding_mask_bool, quality_arr_lengths, max_len, n_reads
    else:
        padded = [q.ljust(max_len, b'\x00') for q in quality_list]
        joined = b''.join(padded)
        array = np.frombuffer(joined, dtype=np.int8).reshape(n_reads, max_len) 
        padding_mask_bool = array == ord('\x00')
        row_tilde_count = np.sum(padding_mask_bool, axis=1)
        chunk_padding_bool = bool(np.any(row_tilde_count > 0))
        array = np.where(padding_mask_bool, array, array - phred_offset)
        return array, chunk_padding_bool, row_tilde_count, padding_mask_bool, quality_arr_lengths, max_len, n_reads

def seq_to_array(sequence_list, chunk_padding_bool, max_len, n_reads):
    """
    Converts a list of nucleotide sequence strings into a 2D numpy array.

    Concatenates all sequences and reinterprets the raw bytes as an
    array of ASCII codes, reshaping into a (n_reads, read_length) matrix 
    for vectorized downstream processing.

    Args:
        sequence_list (list[bytes]): Sequence byte-strings; equal length
            unless chunk_padding_bool is True.
        chunk_padding_bool (bool): If True, sequences differ in length and
            are right-padded with null bytes to max_len before conversion.
        max_len (int): Length to pad/reshape each sequence to.
        n_reads (int): Number of sequences (rows) in sequence_list.

    Returns:
        numpy.ndarray: Signed 8-bit integer array of shape
            (n_reads, max_len) of ASCII character codes.
    """
    if not chunk_padding_bool:
        joined = b''.join(sequence_list)
        array = np.frombuffer(joined, dtype=np.int8).reshape(n_reads, max_len)
        return array
    else:
        padded = [s.ljust(max_len, b'\x00') for s in sequence_list]
        joined = b''.join(padded)
        array = np.frombuffer(joined, dtype=np.int8).reshape(n_reads, max_len)
        return array

def header_mgi_to_illumina(mgi_header, barcode5, barcode7, instrument, run):
    """
    Converts an MGI-style FASTQ header into an Illumina-style header format.

    Takes strings for elements needed in an Illumina header that are missing 
    from the MGI header (instrument, run, barcodes). Parses the MGI coordinate 
    components and maps them into the standard Illumina identifier structure.

    Args:
        mgi_header (str): MGI header, with or without leading '@', in the
            form "<flowcell>L<lane>C<column>R<row><tile>/<read_num>".
        barcode5 (str): The i5 index sequence or barcode string.
        barcode7 (str): The i7 index sequence or barcode string.
        instrument (str): The sequencing instrument identifier.
        run (str): The run identifier or run number.

    Returns:
        str: Illumina-style header (including leading '@').

    Raises:
        ValueError: If the header doesn't match the expected MGI format.
    """
    mgi_header = mgi_header.lstrip(b"@").strip()
    if re.search(rb"^[\w-]+:\d+:[\w-]+:\d+:\d+:\d+:\d+\s+[12]:[YN]:\d+:", mgi_header):
        return b"@" + mgi_header
    strings = re.search(rb"^(\w+)L(\d+)C(\d+)R(\d{3})(\d+)\/([12])", mgi_header)
    if strings is None:
        raise ValueError(f"Header does not match expected MGI format: {mgi_header!r}")
    illumina_header = b"@%b:%b:%b:%d:%d:%d:%d %d:N:0:%b+%b" % (
        instrument.encode('ascii'),
        run.encode('ascii'),
        strings.group(1),
        int(strings.group(2)),
        int(strings.group(5)),
        int(strings.group(3)),
        int(strings.group(4)),
        int(strings.group(6)),
        barcode5.encode('ascii'),
        barcode7.encode('ascii')
    )
    return illumina_header

def build_pipeline(parameters):
    """Constructs a list of active trimming operations based on configuration flags.

    Filters incoming parameter flags and returns a list of executable lambda functions 
    representing active trimming strategies. Each function expects `(sequence_arr, quality_arr)` 
    and returns a tuple of `(left, right)` NumPy boundary arrays.

    Args:
        parameters (dict): Configuration dictionary containing boolean flags 
            (e.g., 'endqual_filter_flag', 'adapter_filter_flag') and strategy-specific 
            trimming arguments.

    Returns:
        list[callable]: A list of functions matching active modules.
    """
    pipeline = []
    #sequence_based
    if parameters.get("kmer_filter_flag"):
        pipeline.append(lambda seq, qual, chunk_padding_bool, row_tilde_count, padding_mask_bool: kmer_complexity_scan(seq, chunk_padding_bool, padding_mask_bool, kmer=parameters["kmer_size"], low_complex_cutoff=parameters["kmer_cutoff"], allow_n=parameters["allow_n_kmer"]))
    if parameters.get("n_trimming_flag"):
        pipeline.append(lambda seq, qual, chunk_padding_bool, row_tilde_count, padding_mask_bool: n_end_trimming(seq, padding_mask_bool, chunk_padding_bool))
    if parameters.get("poly_filter_flag"):
        pipeline.append(lambda seq, qual, chunk_padding_bool, row_tilde_count, padding_mask_bool: homopolymer_nucleotide_trimming(seq, padding_mask_bool, chunk_padding_bool, poly_length_both=parameters["poly_length_both"], poly_length_start=parameters["poly_length_start"], poly_length_end=parameters["poly_length_end"], poly_bases_both=parameters["poly_bases_both"], poly_bases_start=parameters["poly_bases_start"], poly_bases_end=parameters["poly_bases_end"]))
    if parameters.get("adapter_filter_flag"):
        pipeline.append(lambda seq, qual, chunk_padding_bool, row_tilde_count, padding_mask_bool: adapter_trimming(seq, chunk_padding_bool, row_tilde_count, adapter_sequences=parameters["adapter_sequences"], mismatches=parameters["adapter_mismatch"]))

    #quality_based
    if parameters.get("endqual_filter_flag"):
        pipeline.append(lambda seq, qual, chunk_padding_bool, row_tilde_count, padding_mask_bool: trim_ends_quality(qual, chunk_padding_bool, padding_mask_bool, min_quality_both=parameters["min_quality_both"], endqual_min_start=parameters["endqual_min_start"], endqual_min_end=parameters["endqual_min_end"]))
    if parameters.get("minimum_average_qual_pre") > 0:
        pipeline.append(lambda seq, qual, chunk_padding_bool, row_tilde_count, padding_mask_bool: average_quality_filter_wrapper(qual, chunk_padding_bool, row_tilde_count, min_avg_qual=parameters["minimum_average_qual_pre"]))
    if parameters.get("slider_filter_flag"):
        pipeline.append(lambda seq, qual, chunk_padding_bool, row_tilde_count, padding_mask_bool: sliding_window_quality(qual, chunk_padding_bool, padding_mask_bool, row_tilde_count, slider_quality=parameters["slider_quality"], slider_window=parameters["slider_window"], slider_step=parameters["slider_step"]))

    #length_based
    if parameters.get("cut_flag"):
        pipeline.append(lambda seq, qual, chunk_padding_bool, row_tilde_count, padding_mask_bool: cut_set_ends(seq, chunk_padding_bool, row_tilde_count, cut_both=parameters["cut_both"], cut_start=parameters["cut_start"], cut_end=parameters["cut_end"]))

    return pipeline

##### Writing files functions #####

def open_fastq_writer(filepath, output_dir, gzip_output):
    """
    Opens an output file handle for writing filtered FASTQ reads.

    Derives the output filename from the input filepath's basename plus a
    "_filtered" suffix (with a ".gz" extension if gzip_output is True) and 
    opens the file in exclusive binary write mode (`"xb"`).

    Args:
        filepath (str): Path to the original input FASTQ file, used to
            derive the output filename.
        output_dir (str): Directory in which to create the output file.
        gzip_output (bool): If True, names the file with a ".fastq.gz" 
            extension for gzip compression.

    Returns:
        io.BufferedIOBase: An opened binary file handle for the output file.

    Raises:
        FileExistsError: If the output file already exists (due to exclusive create mode).
    """
    basename_for_write = basename_file(filepath)
    extension = "_filtered.fastq.gz" if gzip_output else "_filtered.fastq"
    filename = basename_for_write + extension
    if "rejected" in filename:
        filename = filename.replace("_filtered", "")
    out_filepath = os.path.join(output_dir, filename)
    return open(out_filepath, "xb")

##### Processing reads functions #####
def average_quality_filter_wrapper(quality_arr, chunk_padding_bool, row_tilde_count, min_avg_qual):
    """
    Wraps average_quality_batch to conform to the pipeline's (left, right)
    boundary-array contract, since average_quality_batch is a whole-read
    pass/fail filter rather than a trimmer.
    Reads with average quality >= min_avg_qual are left untouched
    (left=0, right=read_length). Reads that fail the threshold are collapsed
    to a zero-length boundary (left=0, right=0), signaling downstream stages
    to treat them as fully discarded.
    Args:
        quality_arr (numpy.ndarray): (n_reads, read_length) array of
            per-base quality scores.
        chunk_padding_bool (bool): Whether quality_arr contains padded
            (unequal-length) reads.
        row_tilde_count (numpy.ndarray | None): (n_reads,) count of padding
            bytes per read, used to recover each read's real length when
            chunk_padding_bool is True.
        min_avg_qual (float): Minimum average quality required to keep a read.
    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: (left_cutoffs, right_cutoffs),
            each of shape (n_reads,) and dtype int16.
    """
    n_reads, length = quality_arr.shape
    if not chunk_padding_bool:
        real_lengths = np.full(n_reads, length, dtype=np.int64)
    else:
        real_lengths = length - row_tilde_count
    avg_quals = average_quality_batch(quality_arr, lefts=0, rights=real_lengths)
    passed = avg_quals >= min_avg_qual
    right_cutoffs = np.where(passed, real_lengths, 0).astype(np.int64)
    return np.zeros(n_reads, dtype=np.int8), right_cutoffs

def trim_ends_quality(quality_arr, chunk_padding_bool, padding_mask_bool, min_quality_both, endqual_min_start, endqual_min_end):
    """
    Determines per-read trim boundaries based on quality thresholds at each end.

    For each read, finds the position meeting `endqual_min_start` at the start,
    and the position meeting `endqual_min_end` at the end. Vectorized
    across all reads in a chunk. Reads with no position meeting the threshold
    are assigned zero-length cutoffs. If `endqual_filter_flag` is False,
    returns full-length arrays without trimming.

    Args:
        quality_arr (numpy.ndarray): (n_reads, read_length) array of per-base quality scores.
        chunk_padding_bool (bool): Whether quality_arr contains padded
            (unequal-length) reads.
        padding_mask_bool (numpy.ndarray | None): (n_reads, read_length)
            boolean mask marking right-padding positions, or None if the
            chunk isn't padded. Padded positions are excluded from the scan.
        min_quality_both (int | None): Minimum quality fallback to keep bases from both ends.
        endqual_min_start (int | None): Minimum quality to keep bases from the start of the read.
        endqual_min_end (int | None): Minimum quality to keep bases from the end of the read.

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: A tuple of `(start_cutoffs, end_cutoffs)` 
            arrays of shape `(n_reads,)` and dtype `int16`, giving the left and right trim 
            boundaries per read.
    """
    n_reads, length = quality_arr.shape
    if endqual_min_start is None:
        endqual_min_start = min_quality_both if min_quality_both is not None else 0
    if endqual_min_end is None:
        endqual_min_end = min_quality_both if min_quality_both is not None else 0

    if not chunk_padding_bool:
        qual_mask = quality_arr >= endqual_min_start
        start_cutoffs = qual_mask.argmax(axis=1)
        zero_rows = start_cutoffs == 0
        if zero_rows.any():
            start_good_pos = qual_mask[:, 0] | ~zero_rows
            start_cutoffs = np.where(start_good_pos, start_cutoffs, length)

        quality_arr_rev = quality_arr[:, ::-1]
        qual_mask = quality_arr_rev >= endqual_min_end
        end_cutoffs = length - qual_mask.argmax(axis=1)
        zero_end_rows = end_cutoffs == length
        if zero_end_rows.any():
            end_good_pos = qual_mask[:, 0] | ~zero_end_rows
            end_cutoffs = np.where(end_good_pos, end_cutoffs, 0)
    else:
        qual_mask = (quality_arr >= endqual_min_start) & ~padding_mask_bool
        start_cutoffs = qual_mask.argmax(axis=1)
        zero_rows = start_cutoffs == 0
        if zero_rows.any():
            start_good_pos = qual_mask[:, 0] | ~zero_rows
            start_cutoffs = np.where(start_good_pos, start_cutoffs, length)

        quality_arr_rev = quality_arr[:, ::-1]
        pad_mask_rev = padding_mask_bool[:, ::-1]
        qual_mask = (quality_arr_rev >= endqual_min_end) & ~pad_mask_rev
        first_good_pos = qual_mask.argmax(axis=1)

        tilde_count_per_row = np.sum(padding_mask_bool, axis=1)
        real_lengths = length - tilde_count_per_row
        run_length = first_good_pos - tilde_count_per_row
        end_cutoffs = real_lengths - run_length

        zero_end_rows = first_good_pos == 0
        if zero_end_rows.any():
            end_good_pos = qual_mask[:, 0] | ~zero_end_rows
            end_cutoffs = np.where(end_good_pos, end_cutoffs, 0)

    return start_cutoffs.astype(np.int64), end_cutoffs.astype(np.int64)

def homopolymer_nucleotide_trimming(sequence_arr, padding_mask_bool, chunk_padding_bool, poly_length_both, poly_length_start, poly_length_end, poly_bases_both, poly_bases_start, poly_bases_end):
    """
    Determines per-read trim boundaries to remove homopolymer runs 
    from the start and/or end of each read independently.

    Scans for specified base runs (e.g., poly-G tails or custom nucleotides) 
    that meet or exceed defined length thresholds separately for both ends 
    of the sequence array. Specified bases are checked independently 
    rather than on a heteropolymer basis.

    Args:
        sequence_arr (numpy.ndarray): (n_reads, read_length) array of
            per-base ASCII sequence codes.
        padding_mask_bool (numpy.ndarray | None): (n_reads, read_length)
            boolean mask marking right-padding positions, or None if the
            chunk isn't padded.
        chunk_padding_bool (bool): Whether sequence_arr contains padded
            (unequal-length) reads.
        poly_length_both (int): Fallback minimum run length to trigger 
            trimming on either end.
        poly_length_start (int): Minimum run length to trigger trimming 
            at the start of the read.
        poly_length_end (int): Minimum run length to trigger trimming 
            at the end of the read.
        poly_bases_both (str | None): Comma-separated bases to check independently 
            at both ends.
        poly_bases_start (str | None): Comma-separated bases to check independently 
            at the start.
        poly_bases_end (str | None): Comma-separated bases to check independently 
            at the end.

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: (left_cutoffs, right_cutoffs),
            each of shape (n_reads,) and dtype int16, giving the left and right 
            trim boundaries per read.
    """
    n_reads, length = sequence_arr.shape
    if not poly_bases_start and not poly_bases_end and not poly_bases_both:
        return np.zeros(n_reads, dtype=np.int8), np.full(n_reads, length, dtype=np.int64)
    if poly_length_both == poly_length_start == poly_length_end == 0:
        return np.zeros(n_reads, dtype=np.int8), np.full(n_reads, length, dtype=np.int64)
    
    start_bases = []
    end_bases = []
    
    if poly_bases_both is not None:
        bases = [b.strip().upper() for b in poly_bases_both.split(",") if b.strip()]
        if any(len(b) != 1 for b in bases):
            raise ValueError(f"Invalid base entry in '{poly_bases_both}': All bases must be single characters.")
        start_bases = bases
        end_bases = bases
    
    if poly_bases_start is not None:
        start_bases = [b.strip().upper() for b in poly_bases_start.split(",") if b.strip()]
        if any(len(b) != 1 for b in start_bases):
            raise ValueError(f"Invalid base entry in '{poly_bases_start}': All bases must be single characters.")
    
    if poly_bases_end is not None:
        end_bases = [b.strip().upper() for b in poly_bases_end.split(",") if b.strip()]
        if any(len(b) != 1 for b in end_bases):
            raise ValueError(f"Invalid base entry in '{poly_bases_end}': All bases must be single characters.")
    
    poly_length_start = poly_length_start if poly_length_start != 0 else poly_length_both
    poly_length_end = poly_length_end if poly_length_end != 0 else poly_length_both
    
    right_cutoffs = np.full(n_reads, length, dtype=np.int64)
    left_cutoffs = np.zeros(n_reads, dtype=np.int8)
    
    for base in start_bases:
        base_code = ord(base)
        non_base_mask = sequence_arr != base_code
        padded_mask = np.column_stack([non_base_mask, np.ones(n_reads, dtype=bool)])
        first_non_pos = padded_mask.argmax(axis=1)
        trim_amount = np.where(first_non_pos >= poly_length_start, first_non_pos, 0)
        left_cutoffs = np.maximum(left_cutoffs, trim_amount)
    
    if end_bases:
        rev_seq = np.ascontiguousarray(sequence_arr[:, ::-1])    
        if not chunk_padding_bool:
            for base in end_bases:
                base_code = ord(base)
                non_base_mask = rev_seq != base_code
                padded_mask = np.column_stack([non_base_mask, np.ones(n_reads, dtype=bool)])
                first_non_pos = padded_mask.argmax(axis=1)
                trim_amount = np.where(first_non_pos >= poly_length_end, first_non_pos, 0)
                base_right_cutoffs = length - trim_amount
                right_cutoffs = np.minimum(right_cutoffs, base_right_cutoffs)
        else:
            tilde_count_per_row = np.sum(padding_mask_bool, axis=1)
            real_length = length - tilde_count_per_row
            for base in end_bases:
                non_base_mask = (rev_seq != ord(base)) & (rev_seq != ord('\x00'))
                padded_mask = np.column_stack([non_base_mask, np.ones(n_reads, dtype=bool)])
                first_non_pos = padded_mask.argmax(axis=1)
                run_length = first_non_pos - tilde_count_per_row
                trim_amount = np.where(run_length >= poly_length_end, run_length, 0)
                base_right_cutoffs = real_length - trim_amount
                right_cutoffs = np.minimum(right_cutoffs, base_right_cutoffs)
    
    return left_cutoffs, right_cutoffs

def n_end_trimming(sequence_arr, padding_mask_bool, chunk_padding_bool):
    """
    Removes leading and trailing N-bases from each read.

    This function detects and strips runs of N's from both the ends
    of each read (any length >= 1), leaving internal N's untouched.
    Internally reuses `homopolymer_nucleotide_trimming` under the hood
    configured for N bases.

    Args:
        sequence_arr (numpy.ndarray): (n_reads, read_length) array of
            per-base ASCII sequence codes.
        padding_mask_bool (numpy.ndarray | None): (n_reads, read_length)
            boolean mask marking right-padding positions, or None if the
            chunk isn't padded.
        chunk_padding_bool (bool): Whether sequence_arr contains padded
            (unequal-length) reads.

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: (left_cutoffs, right_cutoffs),
            each of shape (n_reads,), giving the trim boundaries per read.
    """
    lefts, rights = homopolymer_nucleotide_trimming(sequence_arr, padding_mask_bool = padding_mask_bool, chunk_padding_bool = chunk_padding_bool, poly_length_both = 1, poly_length_start = 0, poly_length_end = 0, poly_bases_both = "N", poly_bases_start = None, poly_bases_end = None)
    return lefts, rights

def cut_set_ends(sequence_arr, chunk_padding_bool, row_tilde_count, cut_both, cut_start, cut_end):
    """
    Produces fixed, user-specified trim boundaries applied uniformly to
    every read (e.g. hard-trimming a known number of adapter/primer bases
    off each end, regardless of quality).

    `cut_both` sets a symmetric default trim for both ends, but is overridden
    on either side individually if `cut_start` and/or `cut_end` are also
    given (nonzero) — so a user can specify a general trim amount while
    still customizing one end specifically. A value of 0 for `cut_start` or
    `cut_end` is treated as "not given" and falls back to `cut_both`, not
    as an explicit request to trim nothing off that end.

    Resulting boundaries are clamped to a minimum of 0 — if the requested
    trim would push a boundary negative (e.g. `cut_end` larger than the
    read, or a left boundary computed past the right boundary), it is
    floored at 0 instead.

    Args:
        sequence_arr (numpy.ndarray): (n_reads, read_length) array of
            per-base ASCII sequence codes, used only for its shape.
        chunk_padding_bool (bool): Whether sequence_arr contains padded
            (unequal-length) reads.
        row_tilde_count (numpy.ndarray): (n_reads,) count of padding bytes
            per read, used to recover each read's real length when
            chunk_padding_bool is True.
        cut_both (int): Number of bases to trim off both ends. Used as a
            fallback for any side left at 0 (i.e. not given explicitly)
            via cut_start/cut_end.
        cut_start (int): Number of bases to trim from the 5' end. Takes
            priority over cut_both if nonzero.
        cut_end (int): Number of bases to trim from the 3' end. Takes
            priority over cut_both if nonzero.

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: (left_cutoffs, right_cutoffs),
            each of shape (n_reads,), giving the left and right trim
            boundaries per read, clamped to [0, read_length].
    """
    n_reads, length = sequence_arr.shape
    cut_start = cut_start if cut_start != 0 else cut_both
    cut_end = cut_end if cut_end != 0 else cut_both
    if not chunk_padding_bool:
        cut_end_pos = length - cut_end
        cut_end_pos = max(cut_end_pos, 0)
        cut_start_pos = min(cut_start, cut_end_pos)
        cut_start_pos = max(cut_start_pos, 0)
        return np.full(n_reads, cut_start_pos, dtype=np.int64), np.full(n_reads, cut_end_pos, dtype=np.int64)
    else:
        real_length = length - row_tilde_count
        end_positions = real_length - cut_end
        cut_start_pos = np.where(cut_start > end_positions, end_positions, cut_start)
        return np.full(n_reads, cut_start_pos, dtype=np.int64), end_positions

def sliding_window_quality(quality_arr, chunk_padding_bool, padding_mask_bool, row_tilde_count, slider_quality, slider_window, slider_step):
    """
    Determines per-read trim boundaries using a sliding-window quality scan,
    finding the longest (and highest-average-quality, as tiebreaker) stretch
    of the read where every window of `slider_window` bases has a mean
    quality above `slider_quality`.

    Rather than trimming only from the ends inward, it identifies the best
    surviving internal stretch of acceptable quality and reports its
    boundaries.

    Args:
        quality_arr (numpy.ndarray): (n_reads, read_length) array of
            per-base quality scores.
        chunk_padding_bool (bool): Whether quality_arr contains padded
            (unequal-length) reads.
        padding_mask_bool (numpy.ndarray | None): (n_reads, read_length)
            boolean mask marking right-padding positions, or None if the
            chunk isn't padded. Padded positions are always treated as
            failing windows.
        row_tilde_count (numpy.ndarray | None): (n_reads,) count of padding
            bytes per read, used to recover each read's real length when
            chunk_padding_bool is True.
        slider_quality (int): Minimum acceptable mean quality within a window.
        slider_window (int): Number of bases per sliding window.
        slider_step (int): Step size between successive window start positions.

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: (left_cutoffs, right_cutoffs),
            each of shape (n_reads,), giving the best surviving
            [left, right) region per read. Reads with no failing windows keep
            their full length; reads that fail everywhere get a zero-length region.
    """
    n_reads, length = quality_arr.shape

    if not chunk_padding_bool:
        if length < slider_window:
            return np.zeros(n_reads, dtype=np.int16), np.full(n_reads, length, dtype=np.int16)

        last_possible_start = length - slider_window
        window_starts = np.arange(0, last_possible_start + 1, slider_step)
        if window_starts[-1] != last_possible_start:
            window_starts = np.append(window_starts, last_possible_start)

        cumsum = np.cumsum(quality_arr, axis=1, dtype=np.int64)
        cumsum = np.concatenate([np.zeros((n_reads, 1), dtype=np.int64), cumsum], axis=1)
        window_sums = cumsum[:, window_starts + slider_window] - cumsum[:, window_starts]
        failed_mask = window_sums < (slider_quality * slider_window)

        bad_positions = np.zeros((n_reads, length), dtype=bool)
        n_windows = len(window_starts)

        if slider_step == 1:
            for offset in range(slider_window):
                bad_positions[:, offset:offset + n_windows] |= failed_mask
        else:
            for j, start in enumerate(window_starts):
                bad_positions[:, start:start + slider_window] |= failed_mask[:, j:j + 1]

        real_lengths = np.full(n_reads, length, dtype=np.int64)
    else:
        real_lengths = length - row_tilde_count
        too_short = real_lengths < slider_window

        if length < slider_window:
            left_cutoffs = np.zeros(n_reads, dtype=np.int16)
            right_cutoffs = real_lengths.astype(np.int64)
            return left_cutoffs, right_cutoffs

        cumsum = np.cumsum(quality_arr, axis=1, dtype=np.int64)
        cumsum = np.concatenate([np.zeros((n_reads, 1), dtype=np.int64), cumsum], axis=1)
        pad_cumsum = np.cumsum(padding_mask_bool.astype(np.int64), axis=1)
        pad_cumsum = np.concatenate([np.zeros((n_reads, 1), dtype=np.int64), pad_cumsum], axis=1)

        last_possible_start = length - slider_window
        window_starts = np.arange(0, last_possible_start + 1, slider_step)
        if window_starts[-1] != last_possible_start:
            window_starts = np.append(window_starts, last_possible_start)

        window_sums = cumsum[:, window_starts + slider_window] - cumsum[:, window_starts]
        window_pad_counts = pad_cumsum[:, window_starts + slider_window] - pad_cumsum[:, window_starts]
        window_has_pad = window_pad_counts > 0

        failed_mask = (window_sums < (slider_quality * slider_window)) | window_has_pad

        bad_positions = np.zeros((n_reads, length), dtype=bool)
        n_windows = len(window_starts)
        if slider_step == 1:
            for offset in range(slider_window):
                bad_positions[:, offset:offset + n_windows] |= failed_mask
        else:
            for j, start in enumerate(window_starts):
                bad_positions[:, start:start + slider_window] |= failed_mask[:, j:j + 1]

        bad_positions |= padding_mask_bool

    good_positions = ~bad_positions
    no_bad = ~bad_positions.any(axis=1)
    all_bad = bad_positions.all(axis=1)
    left_cutoffs = np.zeros(n_reads, dtype=np.int64)
    right_cutoffs = np.zeros(n_reads, dtype=np.int64)
    right_cutoffs[no_bad] = real_lengths[no_bad].astype(np.int64)
    needs_stretch_search = ~no_bad & ~all_bad
    if not needs_stretch_search.any():
        if chunk_padding_bool:
            left_cutoffs[too_short] = 0
            right_cutoffs[too_short] = real_lengths[too_short].astype(np.int64)
        return left_cutoffs, right_cutoffs

    padded = np.zeros((needs_stretch_search.sum(), length + 2), dtype=bool)
    padded[:, 1:-1] = good_positions[needs_stretch_search]
    diffs = np.diff(padded.view(np.int8), axis=1)
    local_rows, start_cols = np.where(diffs == 1)
    _, end_cols = np.where(diffs == -1)
    global_rows = np.where(needs_stretch_search)[0][local_rows]
    run_lengths = end_cols - start_cols
    
    run_sums = cumsum[global_rows, end_cols] - cumsum[global_rows, start_cols]
    run_means = run_sums / run_lengths
    order = np.lexsort((-run_means, -run_lengths, global_rows))
    sorted_rows = global_rows[order]
    first_in_group = np.concatenate([[0], np.flatnonzero(sorted_rows[1:] != sorted_rows[:-1]) + 1])
    best_idx = order[first_in_group]

    left_cutoffs[global_rows[best_idx]] = start_cols[best_idx]
    right_cutoffs[global_rows[best_idx]] = end_cols[best_idx]

    if chunk_padding_bool:
        left_cutoffs[too_short] = 0
        right_cutoffs[too_short] = real_lengths[too_short].astype(np.int16)

    return left_cutoffs, right_cutoffs

def adapter_trimming(sequence_arr, chunk_padding_bool, row_tilde_count, adapter_sequences, mismatches):
    """
    Determines per-read trim boundaries to remove specific adapter sequences.
    For each read, finds the earliest position where any provided adapter
    sequence appears (as an exact substring if mismatches <= 0, or as an
    approximate match allowing up to `mismatches` substitutions if
    mismatches > 0, via vectorized Hamming-distance search), and trims the
    read at that position. 
    Args:
        sequence_arr (numpy.ndarray): (n_reads, read_length) array of
            per-base ASCII sequence codes (uint8).
        chunk_padding_bool (bool): Whether sequence_arr contains padded
            (unequal-length) reads.
        row_tilde_count (numpy.ndarray | None): (n_reads,) count of padding
            bytes per read, used to keep matches within each read's real
            length when chunk_padding_bool is True.
        adapter_sequences (list[bytes]): List of adapter byte-sequences to
            search for. 
        mismatches (int): Number of allowed mismatches (substitutions) when
            searching for adapters. If = 0, uses exact substring matching.
            If > 0, uses vectorized Hamming-distance fuzzy matching
            (substitutions only, no indels).
    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: (left_cutoffs, right_cutoffs),
            each of shape (n_reads,) and dtype int16, giving the left and right
            trim boundaries per read. Left cutoffs are always 0 (3'-end trimming only).
    """
    n_reads, length = sequence_arr.shape
    all_bytes = sequence_arr.tobytes()
    right_cutoffs = np.zeros(n_reads, dtype=np.int64)
    if mismatches == 0:
        if not chunk_padding_bool:
            for i in range(n_reads):
                row_bytes = all_bytes[i*length:(i+1)*length]
                best = length
                for adapter_bytes in adapter_sequences:
                    pos = row_bytes.find(adapter_bytes)
                    if pos != -1 and pos < best:
                        best = pos
                right_cutoffs[i] = best
        else:
            real_lengths = length - row_tilde_count
            for i in range(n_reads):
                row_bytes = all_bytes[i*length:(i+1)*length]
                real_len = real_lengths[i]
                row_bytes_real = row_bytes[:real_len]
                best = real_len
                for adapter_bytes in adapter_sequences:
                    pos = row_bytes_real.find(adapter_bytes)
                    if pos != -1 and pos < best:
                        best = pos
                right_cutoffs[i] = best
    else:
        right_cutoffs[:] = length if not chunk_padding_bool else (length - row_tilde_count)

        if not chunk_padding_bool:
            for adapter_bytes in adapter_sequences:
                adapter_arr = np.frombuffer(adapter_bytes, dtype=np.uint8)
                windows = sliding_window_view(sequence_arr, window_shape=len(adapter_arr), axis=1)
                mismatch_matrix = (windows != adapter_arr).sum(axis=2)
                valid_mask = mismatch_matrix <= mismatches
                has_match = valid_mask.any(axis=1)

                if has_match.any():
                    first_match_col = valid_mask.argmax(axis=1)
                    right_cutoffs[has_match] = np.minimum(
                        right_cutoffs[has_match],
                        first_match_col[has_match].astype(np.int16)
                    )
        else:
            real_lengths = length - row_tilde_count
            for adapter_bytes in adapter_sequences:
                adapter_arr = np.frombuffer(adapter_bytes, dtype=np.uint8)
                length_adapter = len(adapter_arr)
                windows = sliding_window_view(sequence_arr, window_shape=length_adapter, axis=1)
                mismatch_matrix = (windows != adapter_arr).sum(axis=2)
                window_ends = np.arange(length_adapter, length + 1)
                within_bounds = window_ends <= real_lengths[:, None]
                valid_mask = (mismatch_matrix <= mismatches) & within_bounds
                has_match = valid_mask.any(axis=1)

                if has_match.any():
                    first_match_col = valid_mask.argmax(axis=1)
                    right_cutoffs[has_match] = np.minimum(
                        right_cutoffs[has_match],
                        first_match_col[has_match].astype(np.int16)
                    )
    return np.zeros(n_reads, dtype=np.int16), right_cutoffs
    
def average_quality_batch(quality_arr, lefts, rights):
    """
    Computes the mean quality score within a per-read [left, right) window,
    vectorized across all reads at once.

    Args:
        quality_arr (numpy.ndarray): (n_reads, read_length) array of
            per-base quality scores.
        lefts (numpy.ndarray): Per-read left boundary (inclusive), shape (n_reads,).
        rights (numpy.ndarray): Per-read right boundary (exclusive), shape (n_reads,).

    Returns:
        numpy.ndarray: Per-read mean quality within [left, right), shape
            (n_reads,). Reads with an empty window (right <= left) get 0.0.
    """
    n_reads, length = quality_arr.shape
    cumsum = np.empty((n_reads, length + 1), dtype=np.int32)
    cumsum[:, 0] = 0
    np.cumsum(quality_arr, axis=1, dtype=np.int32, out=cumsum[:, 1:])
    row_indices = np.arange(n_reads)
    sums = cumsum[row_indices, rights] - cumsum[row_indices, lefts]
    counts = rights - lefts
    avg_qualities = np.divide(sums, counts, out=np.zeros_like(sums, dtype=np.float64), where=counts > 0)
    return avg_qualities

def kmer_complexity_scan(sequence_arr, chunk_padding_bool, padding_mask_bool, kmer, low_complex_cutoff, allow_n):
    """
    Counts k-mers per read and flags reads falling below the specified complexity
    cutoff for removal.

    Args:
        sequence_arr (numpy.ndarray): A 2D array of ASCII sequence codes
            of shape (n_reads, length).
        chunk_padding_bool (bool): Whether sequence_arr contains padded
            (unequal-length) reads.
        padding_mask_bool (numpy.ndarray | None): (n_reads, length) boolean
            mask marking right-padding positions, or None if the chunk
            isn't padded. K-mer windows overlapping padding are excluded.
        kmer (int | str | iterable): The length or lengths of the k-mers to evaluate.
        low_complex_cutoff (float): The percentage threshold of unique k-mers 
            relative to the maximum possible windows. If the ratio falls below 
            this value, the read is flagged as low complexity.
        allow_n (bool): If True, encodes 'N' using a 3-bit system (base-8); 
            if False, uses a 2-bit system strictly for 'A', 'C', 'G', and 'T' (base-4).

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: A tuple containing two 1D NumPy arrays of dtype `np.int32`:
            - First array: An array of zeros of shape (n_reads,) acting as primary status flags.
            - Second array: An array of shape (n_reads,) where complex reads retain their 
              original length value and low-complexity reads are set to 0.

    Raises:
        ValueError: If any chosen `kmer` length is greater than the sequence `length`.
    """
    n_reads, length = sequence_arr.shape
    if isinstance(kmer, str):
        kmer_list = [int(k.strip()) for k in kmer.split(',')]
    elif hasattr(kmer, '__iter__'):
        kmer_list = [int(k) for k in kmer]
    else:
        kmer_list = [int(kmer)]

    for k in kmer_list:
        if k > length:
            raise ValueError(f"k-mer length {k} is greater than sequence length {length}")
        if k > 21:
            raise ValueError(
                f"Chosen k-mer length {k} is too big, maximum safe k-mer is 21. "
            )


    global_passed = np.ones(n_reads, dtype=bool)
    mapping = np.zeros(256, dtype=np.int8)
    if allow_n:
        mapping[ord('A')] = 0
        mapping[ord('C')] = 1
        mapping[ord('G')] = 2
        mapping[ord('T')] = 3
        mapping[ord('N')] = 4
        mapping[ord('\x00')] = 5
        bits_per_base = 3
    else:
        mapping[ord('A')] = 0
        mapping[ord('C')] = 1
        mapping[ord('G')] = 2
        mapping[ord('T')] = 3
        mapping[ord('\x00')] = 4
        bits_per_base = 3
    int_matrix_full = mapping[sequence_arr]

    for k in kmer_list:
        if not np.any(global_passed):
            break
        if k > length:
            raise ValueError(f"k-mer length {k} is greater than sequence length {length}")
            
        max_kmers = length - k + 1
        alphabet_size = 5 if allow_n else 4
        max_possible_unique = min(max_kmers, alphabet_size**k)
        
        kmer_ints = np.zeros((n_reads, max_kmers), dtype=np.int64)
        
        if not chunk_padding_bool:
            for i in range(k):
                kmer_ints = (kmer_ints << bits_per_base) | int_matrix_full[:, i:i+max_kmers]

            sorted_kmers = np.sort(kmer_ints, axis=1)
            is_new = np.empty_like(sorted_kmers, dtype=bool)
            is_new[:, 0] = True
            np.not_equal(sorted_kmers[:, 1:], sorted_kmers[:, :-1], out=is_new[:, 1:])
            unique_counts = is_new.sum(axis=1)

            ratio = unique_counts / max_possible_unique
            global_passed &= (ratio >= (low_complex_cutoff / 100))
        else:
            window_has_pad = np.zeros((n_reads, max_kmers), dtype=bool)
            for i in range(k):
                kmer_ints = (kmer_ints << bits_per_base) | int_matrix_full[:, i:i+max_kmers]
                window_has_pad |= padding_mask_bool[:, i:i+max_kmers]
            kmer_ints = np.where(window_has_pad, -1, kmer_ints)

            sorted_kmers = np.sort(kmer_ints, axis=1)
            is_new = np.empty_like(sorted_kmers, dtype=bool)
            is_new[:, 0] = True
            np.not_equal(sorted_kmers[:, 1:], sorted_kmers[:, :-1], out=is_new[:, 1:])
            unique_counts = is_new.sum(axis=1)

            valid_kmer_counts = (~window_has_pad).sum(axis=1)
            any_pad_in_row = window_has_pad.any(axis=1)
            unique_counts = unique_counts - any_pad_in_row.astype(np.int64)

            # Cap the valid k-mer denominator per read by alphabet_size**k
            denom = np.minimum(valid_kmer_counts, alphabet_size**k)

            ratio = np.where(denom > 0, unique_counts / denom, 0.0)
            global_passed &= (ratio >= (low_complex_cutoff / 100))

    second_array = np.where(global_passed, length, 0).astype(np.int64)
    return np.zeros(n_reads, dtype=np.int16), second_array

##### Unpaired reads workflow functions #####
def process_unpaired_chunk(chunk, phred_offset, minimum_average_qual_post, gzip_output, gzip_level, min_length_output, max_length_output, min_length_output_perc, max_length_output_perc, write_rejected, parameters):
    """
    Validates, quality-trims, and length/quality-filters a chunk of unpaired
    FASTQ reads.

    Combines multiple independent trimming strategies (quality-threshold end
    trimming, sliding-window quality trimming, poly-G/poly-X tail trimming,
    adapter trimming, fixed-position end trimming, N-end trimming, and k-mer
    complexity scanning) by taking the most conservative (innermost) boundary
    from each, then discards reads that fall outside the acceptable length or
    average-quality range after trimming.

    Args:
        chunk (list[bytes]): FIELD_SEP-joined (header, sequence, plus,
            quality) records, as yielded by `lazy_fastq`.
        phred_offset (int): Phred encoding offset (33 or 64).
        minimum_average_qual_post (float): Minimum acceptable mean quality
            after trimming. 0 disables this filter.
        gzip_output (bool): If True, compresses output records with Isal gzip.
        gzip_level (int): Isal gzip compression level (0-3).
        min_length_output (int | None): Minimum acceptable read length
            after trimming, as an absolute count.
        max_length_output (int | None): Maximum acceptable read length
            after trimming, as an absolute count.
        min_length_output_perc (int | None): Minimum acceptable read length
            after trimming, as a percentage of the read's input length.
            Ignored if min_length_output is given.
        max_length_output_perc (int | None): Maximum acceptable read length
            after trimming, as a percentage of the read's input length.
            Ignored if max_length_output is given.
        write_rejected (bool): If True, also collect rejected/discarded
            records (invalid on input or filtered out after trimming).
        parameters (dict): Dictionary of configuration parameters, used to
            build the trimming pipeline and to resolve N-filtering,
            MGI-to-Illumina header conversion, and Phred re-encoding.

    Returns:
        tuple[bytes, int, int, bytes | list]: A tuple containing:
            - Formatted and optionally gzipped FASTQ records for surviving reads.
            - Count of kept reads.
            - Count of rejected reads (invalid input records plus records
              filtered out after trimming).
            - Rejected records: optionally gzipped bytes if write_rejected
              is True, otherwise an empty list.
    """
    if not chunk or len(chunk) == 0:
        return b""
    valid_headers = []
    valid_sequences = []
    valid_pluses = []
    valid_qualities = []
    rejected_reads = []
    rejected = 0
    for r in chunk:
        header, sequence, plus, quality = r.split(FIELD_SEP)
        if validate_fastq(header, sequence, plus, quality, n_filter = parameters["n_filter"], min_length_input = parameters["min_length_input"], max_length_input = parameters["max_length_input"]):
            valid_headers.append(header)
            valid_sequences.append(sequence)
            valid_pluses.append(plus)
            valid_qualities.append(quality)
        else:
            rejected += 1
            if write_rejected:
                rejected_reads.append(b"\n".join((header, sequence, plus, quality)) + b"\n")
    if not valid_headers:
        return [], 0, rejected
    if parameters["mgi_convert_flag"]:
        valid_pluses = [plus + b"_OriginalHeader:" + header for plus, header in zip(valid_pluses, valid_headers)]
        valid_headers = [header_mgi_to_illumina(header, parameters["mgi_bc5"], parameters["mgi_bc7"], parameters["mgi_instrument"], parameters["mgi_run"]) for header in valid_headers]
    quality_arr, chunk_padding_bool, row_tilde_count, padding_mask_bool, raw_lengths, max_len, n_reads = qual_to_array(quality_list = valid_qualities, phred_offset = phred_offset)
    sequence_arr = seq_to_array(sequence_list = valid_sequences, chunk_padding_bool = chunk_padding_bool, max_len = max_len, n_reads = n_reads)
    n_reads, length = quality_arr.shape
    left_list = [np.zeros(n_reads, dtype=np.int64)]
    right_list = [np.full(n_reads, length, dtype=np.int64) - row_tilde_count]
    for step in build_pipeline(parameters):
        left, right = step(sequence_arr, quality_arr, chunk_padding_bool, row_tilde_count, padding_mask_bool)
        left_list.append(left)
        right_list.append(right)
    lefts = np.maximum.reduce(left_list)
    rights = np.minimum.reduce(right_list)
    if (min_length_output is not None or max_length_output is not None
            or min_length_output_perc is not None or max_length_output_perc is not None):
        lengths_out = rights - lefts
        if max_length_output is not None:
            effective_max = max_length_output
        elif max_length_output_perc is not None:
            effective_max = (raw_lengths * max_length_output_perc / 100).astype(int)
        else:
            effective_max = np.inf
        if min_length_output is not None:
            effective_min = min_length_output
        elif min_length_output_perc is not None:
            effective_min = (raw_lengths * min_length_output_perc / 100).astype(int)
        else:
            effective_min = 0
        length_mask = (lengths_out <= effective_max) & (lengths_out >= effective_min) & (lengths_out > 0)
    else:
        lengths_out = rights - lefts
        length_mask = lengths_out > 0
    if minimum_average_qual_post > 0:
        average_quals = average_quality_batch(quality_arr, lefts, rights)
        qual_mask = average_quals >= minimum_average_qual_post
    else:
        qual_mask = np.ones(n_reads, dtype=bool)
    keep_mask = length_mask & qual_mask
    if parameters["phred_out"] == 33 and phred_offset == 64:
        qual_table = PHRED64_TO_33
    elif parameters["phred_out"] == 64 and phred_offset == 33:
        qual_table = PHRED33_TO_64
    else:
        qual_table = None
    results = []
    for header, sequence, plus_line, quality, left, right, keep in zip(
            valid_headers, valid_sequences, valid_pluses, valid_qualities, lefts, rights, keep_mask):
        if not keep:
            rejected += 1
            if write_rejected:
                rejected_reads.append(b"\n".join((header, sequence, plus_line, quality)) + b"\n")
            continue
        qual_out = quality[left:right]
        if qual_table is not None:
            qual_out = qual_out.translate(qual_table)
        results.append(b"\n".join((header, sequence[left:right], plus_line, qual_out)) + b"\n")
    len_results = len(results)
    results = b"".join(results)
    if write_rejected:
        rejected_reads = b"".join(rejected_reads)
    if gzip_output:
        try:
            results = gzip.compress(results, compresslevel=gzip_level)
            if write_rejected:
                rejected_reads = gzip.compress(rejected_reads, compresslevel=gzip_level)
        except Exception as e:
            raise RuntimeError(f"Compression failed inside worker: {str(e)}") from None
    return results, len_results, rejected, rejected_reads

def generate_unpaired_tasks(filepaths, chunk_size, parameters, filetype = None):
    """
    Lazily yields individual chunks alongside their filepaths and precomputed parameters,
    allowing a single global pool to process chunks from multiple files concurrently.
    Args:
        filepaths (list[str]): List of paths to unpaired or interleaved FASTQ files.
        chunk_size (int): Number of reads per chunk (record pairs for
            interleaved files, so 2*chunk_size records are read at once).
        parameters (dict): Dictionary of configuration parameters.
        filetype (str | None): Either "unpaired" or "interleaved", selecting
            how records are chunked and tagged in the yielded tasks.
    Yields:
        dict: A task dictionary containing the task type, filepath, chunk data,
            and precomputed metadata.
    """
    for filepath in filepaths:
        start_time = time.monotonic()
        logger.info("%s: started processing.", os.path.basename(filepath))
        if GZIP_DETECTION[filepath]:
            logger.info("%s: gzip format detected.", os.path.basename(filepath))
        else:
            logger.info("%s: text format detected.", os.path.basename(filepath))
        logger.info("%s: file size: %s bytes.", os.path.basename(filepath), os.path.getsize(filepath))
        phred_offset = detect_phred_offset(
            filepath=filepath,
            reads_for_phred_offset=parameters["reads_for_phred_offset"],
            phred_offset=parameters["phred_offset"]
        )
        if ESTIMATED_ZIP_RATIO.get(filepath) is not None:
            logger.info("%s: gzip format detected.", os.path.basename(filepath))
            logger.info("%s: estimated gzip compression ratio: %s.", os.path.basename(filepath), ESTIMATED_ZIP_RATIO.get(filepath))
        logger.info("%s: Phred offset of %s detected.", os.path.basename(filepath), phred_offset)
        logger.info("%s: estimated bytes per read: %s.", os.path.basename(filepath), ESTIMATED_BYTE_PER_READ[filepath])
        logger.info("%s: (estimated) read count: %s.", os.path.basename(filepath), ESTIMATED_READ_COUNTS.get(filepath, "unknown"))
        reads_iter = lazy_fastq(filepath)
        chunk_number = 0
        if filetype == "unpaired":
            while True:
                chunk_number += 1
                chunk = list(itertools.islice(reads_iter, chunk_size))
                if not chunk:
                    elapsed = time.monotonic() - start_time
                    estimated_reads = ESTIMATED_READ_COUNTS.get(filepath)
                    if estimated_reads is not None and elapsed > 0:
                        reads_per_sec = estimated_reads / elapsed
                        logger.info(
                            "%s: finished processing in %.2fs (%.0f reads/sec).",
                            os.path.basename(filepath), elapsed, reads_per_sec)
                    else:
                        logger.info("%s: finished processing in %.2fs.", os.path.basename(filepath), elapsed)
                    break
                yield {
                    "type": "unpaired",
                    "filepath": filepath,
                    "chunk": chunk,
                    "phred_offset": phred_offset
                }
        elif filetype == "interleaved":
            while True:
                chunk_number += 1
                chunk = list(itertools.islice(reads_iter, 2*chunk_size))
                if not chunk:
                    elapsed = time.monotonic() - start_time
                    estimated_reads = ESTIMATED_READ_COUNTS.get(filepath)
                    if estimated_reads is not None and elapsed > 0:
                        reads_per_sec = estimated_reads / elapsed
                        logger.info(
                            "%s: finished processing in %.2fs (%.0f reads/sec).",
                            os.path.basename(filepath), elapsed, reads_per_sec)
                    else:
                        logger.info("%s: finished processing in %.2fs.", os.path.basename(filepath), elapsed)
                    break
                chunk_1 = chunk[0::2]
                chunk_2 = chunk[1::2]
                yield {
                    "type": "paired",
                    "file1": filepath,
                    "file2": filepath,
                    "chunk1": chunk_1,
                    "chunk2": chunk_2,
                    "phred_offset_1": phred_offset,
                    "phred_offset_2": phred_offset,
                    "gzip_output": parameters["gzip_output"], 
                    "gzip_level": parameters["gzip_level"],
                    "discard_singletons": parameters["discard_singletons"]
                }
        else:
            raise TypeError(f"file {filepath} returned filetype {filetype}, which is not recognized.")                
                
def process_unpaired_task_flat(task, parameters):
    """
    Worker wrapper for flat task queue execution on unpaired files.

    Args:
        task (dict): A task dictionary containing chunk data, filepaths, and metadata.
        parameters (dict): Dictionary of configuration parameters.

    Returns:
        tuple: A tuple containing (type, filepath, chunk_results, kept, rejected).
    """
    chunk_results, kept, rejected, rejected_reads = process_unpaired_chunk(
        chunk=task["chunk"],
        phred_offset=task["phred_offset"],
        minimum_average_qual_post=parameters["minimum_average_qual_post"],
        gzip_output = parameters["gzip_output"],
        gzip_level = parameters["gzip_level"],
        min_length_output = parameters["min_length_output"],
        max_length_output = parameters["max_length_output"],
        min_length_output_perc = parameters["min_length_output_perc"],
        max_length_output_perc =parameters["max_length_output_perc"],
        write_rejected= parameters["write_rejected"],
        parameters=parameters
    )
    return task["type"], task["filepath"], chunk_results, kept, rejected, rejected_reads

##### Paired reads workflow funtions #####
def process_paired_task_flat(task, parameters):
    """
    Worker wrapper for flat task queue execution on paired-end files.

    Args:
        task (dict): A task dictionary containing chunks for R1 and R2, filepaths, and metadata.
        parameters (dict): Dictionary of configuration parameters.

    Returns:
        tuple: A tuple containing (type, file1, file2, paired_out_1, paired_out_2,
            R1_singles_out, R2_singles_out, num_paired, num_R1_singles,
            num_R2_singles, rejected_1, rejected_2). R1 and R2 singleton
            records and rejected counts are kept separate since each mate
            is trimmed and filtered independently before reconciliation.
    """
    file1 = task["file1"]
    file2 = task["file2"]
    paired_out_1, paired_out_2, R1_singles_out, R2_singles_out, num_paired, num_R1_singles, num_R2_singles, rejected_1, rejected_2, rejected_R1, rejected_R2 = process_paired_chunk(
        chunks = (task["chunk1"],task["chunk2"] ),
        phred_offset_1=task["phred_offset_1"],
        phred_offset_2=task["phred_offset_2"],
        gzip_output = task["gzip_output"],
        gzip_level = task["gzip_level"],
        discard_singletons=task["discard_singletons"],
        parameters=parameters
    )
    return task["type"], file1, file2, paired_out_1, paired_out_2, R1_singles_out, R2_singles_out, num_paired, num_R1_singles, num_R2_singles, rejected_1, rejected_2, rejected_R1, rejected_R2

def generate_paired_tasks(files, chunk_size, parameters):
    """
    Lazily yields individual paired-end chunks alongside their filepaths and precomputed parameters,
    allowing a single global pool to process chunks from multiple files concurrently.

    Args:
        files (list[tuple[str, str]]): List of (file1, file2) path tuples for paired FASTQ files.
        chunk_size (int): Number of read pairs per chunk.
        parameters (dict): Dictionary of configuration parameters.

    Yields:
        dict: A task dictionary containing the task type, filepaths, chunk data,
            and precomputed metadata.
    """
    for pair in files:
        file1, file2 = pair
        start_time = time.monotonic()
        logger.info("%s and %s: started processing.", os.path.basename(file1), os.path.basename(file2))
        if GZIP_DETECTION[file1]:
            logger.info("%s: gzip format detected.", os.path.basename(file1))
        else:
            logger.info("%s: text format detected.", os.path.basename(file1))
        if GZIP_DETECTION[file2]:
            logger.info("%s: gzip format detected.", os.path.basename(file2))
        else:
            logger.info("%s: text format detected.", os.path.basename(file2))
        logger.info("%s: file size: %s bytes.", os.path.basename(file1), os.path.getsize(file1))
        logger.info("%s: file size: %s bytes.", os.path.basename(file2), os.path.getsize(file2))
        phred_offset_1 = detect_phred_offset(filepath = file1, reads_for_phred_offset = parameters["reads_for_phred_offset"], phred_offset = parameters["phred_offset"])
        phred_offset_2 = detect_phred_offset(filepath = file2, reads_for_phred_offset = parameters["reads_for_phred_offset"], phred_offset = parameters["phred_offset"])
        if ESTIMATED_ZIP_RATIO.get(file1) is not None:
            logger.info("%s: gzip format detected.", os.path.basename(file1))
            logger.info("%s: estimated gzip compression ratio: %s.", os.path.basename(file1), ESTIMATED_ZIP_RATIO.get(file1))
        if ESTIMATED_ZIP_RATIO.get(file2) is not None:
            logger.info("%s: gzip format detected.", os.path.basename(file2))
            logger.info("%s: estimated gzip compression ratio: %s.", os.path.basename(file2), ESTIMATED_ZIP_RATIO.get(file2))
        logger.info("%s: Phred offset of %s detected.", os.path.basename(file1), phred_offset_1)
        logger.info("%s: Phred offset of %s detected.", os.path.basename(file2), phred_offset_2)
        logger.info("%s: (estimated) read count: %s.", os.path.basename(file1), ESTIMATED_READ_COUNTS.get(file1, "unknown"))
        logger.info("%s: (estimated) read count: %s.", os.path.basename(file2), ESTIMATED_READ_COUNTS.get(file2, "unknown"))
        reads_iter_1 = lazy_fastq(file1)
        reads_iter_2 = lazy_fastq(file2)
        while True:
            chunk1 = list(itertools.islice(reads_iter_1, chunk_size))
            chunk2 = list(itertools.islice(reads_iter_2, chunk_size))
            if not chunk1 and not chunk2:
                elapsed = time.monotonic() - start_time
                estimated_reads = ESTIMATED_READ_COUNTS.get(file1)
                if estimated_reads is not None and elapsed > 0:
                    reads_per_sec = estimated_reads / elapsed
                    logger.info(
                        "%s and %s: finished processing in %.2fs (%.0f read pairs/sec).",
                        os.path.basename(file1), os.path.basename(file2), elapsed, reads_per_sec)
                else:
                    logger.info(
                        "%s and %s: finished processing in %.2fs.",
                        os.path.basename(file1), os.path.basename(file2), elapsed
                    )
                return
            if len(chunk1) != len(chunk2):
                raise ValueError(
                    f"Mismatched read counts in paired files {file1} and {file2}. "
                    f"Files must have identical read counts (possible file corruption or truncation)."
                )
            yield {
                "type": "paired",
                "file1": file1,
                "file2": file2,
                "chunk1": chunk1,
                "chunk2": chunk2,
                "phred_offset_1": phred_offset_1,
                "phred_offset_2": phred_offset_2,
                "gzip_output": parameters["gzip_output"], 
                "gzip_level": parameters["gzip_level"]
            }

def trim_reads(records, phred_offset, minimum_average_qual_post, min_length_output, max_length_output, min_length_output_perc, max_length_output_perc, write_rejected, parameters):
    """
    Validates, quality-trims, and length/quality-filters a batch of FASTQ
    reads, keyed by their base (mate-independent) read ID.

    Same trimming logic as `process_unpaired_chunk`, but returns a dict
    keyed by base read ID rather than a flat list of formatted strings —
    this allows the paired workflow to later match up surviving R1/R2 mates
    by ID.

    Args:
        records (list[bytes]): FIELD_SEP-joined (header, sequence, plus,
            quality) records, as yielded by `lazy_fastq`.
        phred_offset (int): Phred encoding offset (33 or 64).
        minimum_average_qual_post (float): Minimum acceptable mean quality
            after trimming. 0 disables this filter.
        min_length_output (int | None): Minimum acceptable read length
            after trimming, as an absolute count.
        max_length_output (int | None): Maximum acceptable read length
            after trimming, as an absolute count.
        min_length_output_perc (int | None): Minimum acceptable read length
            after trimming, as a percentage of the read's input length.
            Ignored if min_length_output is given.
        max_length_output_perc (int | None): Maximum acceptable read length
            after trimming, as a percentage of the read's input length.
            Ignored if max_length_output is given.
        write_rejected (bool): If True, also collect rejected/discarded
            records (invalid on input or filtered out after trimming).
        parameters (dict): Dictionary of configuration parameters, used to
            build the trimming pipeline and to resolve N-filtering,
            MGI-to-Illumina header conversion, and Phred re-encoding.

    Returns:
        tuple[dict[bytes, bytes], int, list[bytes]]: A tuple containing:
            - A dictionary mapping each surviving read's base ID to its
              formatted FASTQ record (not gzip-compressed; compression, if
              any, happens later once R1/R2 mates are reconciled).
            - Total count of rejected reads (invalid input records plus
              records filtered out after trimming).
            - A list of rejected record strings, populated only if
              write_rejected is True.
    """
    valid_headers = []
    valid_sequences = []
    valid_pluses = []
    valid_qualities = []
    rejected_reads = []
    rejected = 0
    for r in records:
        header, sequence, plus, quality = r.split(FIELD_SEP)
        if validate_fastq(header, sequence, plus, quality, n_filter = parameters["n_filter"], min_length_input = parameters["min_length_input"], max_length_input = parameters["max_length_input"]):
            valid_headers.append(header)
            valid_sequences.append(sequence)
            valid_pluses.append(plus)
            valid_qualities.append(quality)
        else:
            rejected += 1
            if write_rejected:
                rejected_reads.append(b"\n".join((header, sequence, plus, quality)) + b"\n")
    if not valid_headers:
        return {}, len(records)
    if parameters["mgi_convert_flag"]:
        valid_pluses = [plus + b"_OriginalHeader:" + header for plus, header in zip(valid_pluses, valid_headers)]
        valid_headers = [header_mgi_to_illumina(header, parameters["mgi_bc5"], parameters["mgi_bc7"], parameters["mgi_instrument"], parameters["mgi_run"]) for header in valid_headers]
    quality_arr, chunk_padding_bool, row_tilde_count, padding_mask_bool, raw_lengths, max_len, n_reads = qual_to_array(quality_list = valid_qualities, phred_offset = phred_offset)
    sequence_arr = seq_to_array(sequence_list = valid_sequences, chunk_padding_bool = chunk_padding_bool, max_len = max_len, n_reads = n_reads)
    n_reads, length = sequence_arr.shape
    left_list = [np.zeros(n_reads, dtype=np.int64)]
    right_list = [np.full(n_reads, length, dtype=np.int64) - row_tilde_count]
    for step in build_pipeline(parameters):
        left, right = step(sequence_arr, quality_arr, chunk_padding_bool, row_tilde_count, padding_mask_bool)
        left_list.append(left)
        right_list.append(right)
    lefts = np.maximum.reduce(left_list)
    rights = np.minimum.reduce(right_list)
    if (min_length_output is not None or max_length_output is not None
            or min_length_output_perc is not None or max_length_output_perc is not None):
        lengths_out = rights - lefts
    
        if max_length_output is not None:
            effective_max = max_length_output
        elif max_length_output_perc is not None:
            effective_max = (raw_lengths * max_length_output_perc / 100).astype(int)
        else:
            effective_max = np.inf
    
        if min_length_output is not None:
            effective_min = min_length_output
        elif min_length_output_perc is not None:
            effective_min = (raw_lengths * min_length_output_perc / 100).astype(int)
        else:
            effective_min = 0
        length_mask = (lengths_out <= effective_max) & (lengths_out >= effective_min) & (lengths_out > 0)
    else:
        lengths_out = rights - lefts
        length_mask = lengths_out > 0
    if minimum_average_qual_post > 0:
        avg_quals = average_quality_batch(quality_arr, lefts, rights)
        qual_mask = avg_quals >= minimum_average_qual_post
    else:
        qual_mask = np.ones(n_reads, dtype=bool)
    keep_mask = length_mask & qual_mask
    if parameters["phred_out"] == 33 and phred_offset == 64:
        qual_table = PHRED64_TO_33
    elif parameters["phred_out"] == 64 and phred_offset == 33:
        qual_table = PHRED33_TO_64
    else:
        qual_table = None
    survivors = {}
    for i, keep in enumerate(keep_mask):
        if not keep:
            rejected += 1
            if write_rejected:
                rejected_reads.append(b"\n".join((valid_headers[i], valid_sequences[i], valid_pluses[i], valid_qualities[i])) + b"\n")
            continue
        left, right = int(lefts[i]), int(rights[i])
        seq_out = valid_sequences[i][left:right]
        qual_out = valid_qualities[i][left:right]
        if qual_table is not None:
            qual_out = qual_out.translate(qual_table)
        base_id, _ = read_info_from_header(valid_headers[i])
        survivors[base_id] = b"\n".join((valid_headers[i], seq_out, valid_pluses[i], qual_out)) + b"\n"
    return survivors, rejected, rejected_reads
        
def process_paired_chunk(chunks, phred_offset_1, phred_offset_2, gzip_output, gzip_level, discard_singletons, parameters):
    """
    Trims and filters one paired chunk of R1/R2 reads, then reconciles the
    two mates by base read ID to determine which reads survive as intact
    pairs versus as orphaned singletons.

    Args:
        chunks (tuple[list, list]): (chunk1, chunk2) — record lists for R1
            and R2 respectively, covering the same reads in the same order.
        phred_offset_1 (int): Phred encoding offset (33 or 64) for R1.
        phred_offset_2 (int): Phred encoding offset (33 or 64) for R2.
        gzip_output (bool): If True, compresses output records with Isal gzip.
        gzip_level (int): Isal gzip compression level (0-3).
        parameters (dict): Dictionary of configuration parameters, including
            `stdout`/`interleaved_out` (which route both mates' surviving
            reads into the R1 output as interleaved records) and
            `write_rejected`.

    Returns:
        tuple[bytes, bytes, bytes, bytes, int, int, int, int, int]: A tuple containing:
            - Surviving R1 records whose R2 mate also survived (optionally gzipped).
            - Surviving R2 records whose R1 mate also survived (optionally gzipped).
            - Surviving R1 records whose mate did not survive, treated as unpaired singletons (optionally gzipped).
            - Surviving R2 records whose mate did not survive, treated as unpaired singletons (optionally gzipped).
            - Count of surviving read pairs.
            - Count of surviving R1 singleton reads.
            - Count of surviving R2 singleton reads.
            - Count of R1 reads rejected during trimming.
            - Count of R2 reads rejected during trimming.
    """
    chunk1, chunk2 = chunks
    survivors_1, rejected_1, rejected_R1 = trim_reads(chunk1, phred_offset_1, minimum_average_qual_post = parameters["minimum_average_qual_post"], min_length_output = parameters["min_length_output"], max_length_output = parameters["max_length_output"], min_length_output_perc = parameters["min_length_output_perc"], max_length_output_perc = parameters["max_length_output_perc"], write_rejected = parameters["write_rejected"], parameters = parameters)
    survivors_2, rejected_2, rejected_R2 = trim_reads(chunk2, phred_offset_2, minimum_average_qual_post = parameters["minimum_average_qual_post"], min_length_output = parameters["min_length_output"], max_length_output = parameters["max_length_output"], min_length_output_perc = parameters["min_length_output_perc"], max_length_output_perc = parameters["max_length_output_perc"], write_rejected = parameters["write_rejected"], parameters = parameters)
    paired_out_1 = []
    paired_out_2 = []
    singles_out_1 = []
    singles_out_2 = []
    for bid, record in survivors_1.items():
        mate = survivors_2.pop(bid, None)
        if mate is not None:
            paired_out_1.append(record)
            paired_out_2.append(mate)
        else:
            singles_out_1.append(record)
    singles_out_2 = list(survivors_2.values())
    num_R1_singles = len(singles_out_1)
    num_R2_singles = len(singles_out_2)
    num_paired = len(paired_out_1)
    
    if discard_singletons:
        singles_out_1 = b""
        singles_out_2 = b""
        
    interleave = parameters["stdout"] or parameters["interleaved_out"]
    if interleave:
        paired_out_1 = b"".join(r1 + r2 for r1, r2 in zip(paired_out_1, paired_out_2))
        paired_out_1 += b"".join(singles_out_1) + b"".join(singles_out_2)
        rejected_R1 = b"".join(rejected_R1) + b"".join(rejected_R2)
        rejected_R2 = b""
        paired_out_2 = b""
        singles_out_1 = b""
        singles_out_2 = b""
    else:
        paired_out_1 = b"".join(paired_out_1)
        paired_out_2 = b"".join(paired_out_2)
        rejected_R1 = b"".join(x.encode("utf-8") if isinstance(x, str) else x for x in rejected_R1)
        rejected_R2 = b"".join(x.encode("utf-8") if isinstance(x, str) else x for x in rejected_R2)
        singles_out_1 = b"".join(singles_out_1)
        singles_out_2 = b"".join(singles_out_2)

    if gzip_output:
        try:
            if paired_out_1:
                paired_out_1 = gzip.compress(paired_out_1, compresslevel=gzip_level)
            if paired_out_2:
                paired_out_2 = gzip.compress(paired_out_2, compresslevel=gzip_level)
            if rejected_R1:
                rejected_R1 = gzip.compress(rejected_R1, compresslevel=gzip_level)
            if rejected_R2:
                rejected_R2 = gzip.compress(rejected_R2, compresslevel=gzip_level)
            if singles_out_1:
                singles_out_1 = gzip.compress(singles_out_1, compresslevel=gzip_level)
            if singles_out_2:
                singles_out_2 = gzip.compress(singles_out_2, compresslevel=gzip_level)
        except Exception as e:
            raise RuntimeError(f"Compression failed inside worker: {str(e)}") from None
    return paired_out_1, paired_out_2, singles_out_1, singles_out_2, num_paired, num_R1_singles, num_R2_singles, rejected_1, rejected_2, rejected_R1, rejected_R2

##### Input handler functions #####
def worker_initilizer(parameters):
    """
    Pool initializer: runs once per worker process at startup, storing
    `parameters` in a worker-global so it doesn't need to be pickled and
    sent again with every individual task.

    Args:
        parameters (dict): Configuration parameters, pickled and sent to
            each worker exactly once when the pool spins it up.
    """
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    global WORKER_PARAMETERS
    WORKER_PARAMETERS = parameters

def unified_worker(task):
    """
    Routes a processing task to the appropriate handler based on its type.
    Serves as a dispatcher for the multiprocessing pool, examining the task's
    type field and delegating to either paired-end or unpaired read processing.

    Args:
        task (dict): A task dictionary containing at minimum:
            - type (str): Either "paired" or "unpaired", indicating which
              processing workflow to apply.
            - parameters (dict): Configuration parameters passed through to the handler.

    Returns:
        tuple: The return value from the appropriate handler function.

    Raises:
        ValueError: If the task type is unknown or missing.
    """
    task_type = task.get("type")
    parameters = WORKER_PARAMETERS

    if task_type == "paired":
        return process_paired_task_flat(task, parameters=parameters)
    elif task_type == "unpaired":
        return process_unpaired_task_flat(task, parameters=parameters)
    else:
        raise ValueError(f"Unknown or missing task type: {task_type}")

def input_handler(unspecified_files, unpaired_files, paired_files, interleaved_files, output_dir, threads, chunk_size, show_progress, stdout, interleaved_out, write_rejected, discard_singletons, parameters):
    """
    Top-level orchestrator that separates input files into paired and
    unpaired groups, sets up file writers, and runs the multiprocessing pool
    workflow across all inputs.

    Args:
        unspecified_files (list[str]): Files with unknown pairing status to auto-detect.
        unpaired_files (list[str]): Explicitly provided unpaired input FASTQ files.
        paired_files (list[tuple[str, str]]): Explicitly provided paired FASTQ file pairs.
        interleaved_files (list[str]): Explicitly provided interleaved FASTQ files.
        output_dir (str): Directory to write output files to.
        threads (int | None): Number of worker threads/processes for the pool.
        chunk_size (int): Number of reads per processing chunk.
        show_progress (bool): If True, renders a live progress bar to stderr.
        stdout (bool): If True, streams surviving reads to stdout instead
            of writing output files.
        interleaved_out (bool): If True, interleaves surviving paired reads
            into a single output file per pair instead of separate R1/R2 files.
        write_rejected (bool): If True, also writes rejected/discarded reads
            to file (ignored when stdout is True).
        parameters (dict): Dictionary of configuration parameters.

    Returns:
        dict: Mapping of file identifiers to their summary statistics.
            For each unpaired file (keyed by its filepath), the value is
            a dict with ``kept`` and ``rejected`` counts. For each paired
            group (keyed by the pair's common filename prefix), the value
            is a dict with ``kept_pairs``, ``kept_R1_singletons``,
            ``kept_R2_singletons``, ``rejected_R1``, and ``rejected_R2``
            counts, since each mate is trimmed and filtered independently
            before reconciliation.
    """
    for file in (unspecified_files or []):
        GZIP_DETECTION[file] = is_gz_file(file)
    for file in (unpaired_files or []):
        GZIP_DETECTION[file] = is_gz_file(file)
    for file in (interleaved_files or []):
        GZIP_DETECTION[file] = is_gz_file(file)
    for f1, f2 in (paired_files or []):
        GZIP_DETECTION[f1] = is_gz_file(f1)
        GZIP_DETECTION[f2] = is_gz_file(f2)
    auto_paired, auto_interleaved, auto_unpaired = find_paired_files(unspecified_files)
    unpaired = auto_unpaired + (unpaired_files or [])
    paired = auto_paired + (paired_files or [])
    interleaved = auto_interleaved + (interleaved_files or [])
    for file in unpaired:
        ESTIMATED_READ_COUNTS[file] = count_reads_estimated(file)
    for file in interleaved:
        ESTIMATED_READ_COUNTS[file] = count_reads_estimated(file)
    for f1, f2 in paired:
        ESTIMATED_READ_COUNTS[f1] = count_reads_estimated(f1)
        ESTIMATED_READ_COUNTS[f2] = count_reads_estimated(f2)
    if show_progress:    
        total_reads = sum(ESTIMATED_READ_COUNTS.values())
        tracker = ProgressTracker(total_reads)
    else:
        class _NullTracker:
            def update(self, n): pass
            def close(self): pass
            def clear_line(self): pass
        tracker = _NullTracker()
    global ACTIVE_PROGRESS_TRACKER
    ACTIVE_PROGRESS_TRACKER = tracker
    file_stats = {}
    file_writing_handles = {}
    for file in unpaired:
        file_stats[file] = {"kept": 0, "rejected": 0}
        if not stdout:
            file_writing_handles[file] = open_fastq_writer(file, output_dir, gzip_output=parameters["gzip_output"])
        if write_rejected:
            file_writing_handles[f"{file}_rejected"] = open_fastq_writer(file, output_dir, gzip_output=parameters["gzip_output"])
    pair_keys = {}
    used_prefixes = set()
    for file in interleaved:
        base_prefix = common_name_parts([os.path.basename(file), os.path.basename(file)])
        pair_keys[(file, file)] = base_prefix
        file_stats[base_prefix] = {"kept_pairs": 0, "kept_R1_singletons": 0, "kept_R2_singletons": 0, "rejected_R1": 0, "rejected_R2": 0}
        if not stdout:
            if interleaved_out:
                key = f"{base_prefix}"
                file_writing_handles[key] = open_fastq_writer(
                    key,
                    output_dir,
                    gzip_output=parameters["gzip_output"]
                    )
                if write_rejected:
                    for suffix in ["_rejected"]:
                        key = f"{base_prefix}{suffix}"
                        file_writing_handles[key] = open_fastq_writer(
                            key,
                            output_dir,
                            gzip_output=parameters["gzip_output"]
                            )
            else:
                suffixes = ["_R1_paired", "_R2_paired"]
                if not discard_singletons:
                    suffixes += ["_R1_unpaired", "_R2_unpaired"]
                for suffix in suffixes:
                    key = f"{base_prefix}{suffix}"
                    file_writing_handles[key] = open_fastq_writer(
                        key,
                        output_dir,
                        gzip_output=parameters["gzip_output"]
                        )
                if write_rejected:
                    for suffix in ["_R1_rejected", "_R2_rejected"]:
                        key = f"{base_prefix}{suffix}"
                        file_writing_handles[key] = open_fastq_writer(
                            key,
                            output_dir,
                            gzip_output=parameters["gzip_output"]
                            )
    for pair in paired:
        file1, file2 = pair
        base_prefix = common_name_parts([os.path.basename(file1), os.path.basename(file2)])
        common_prefix = base_prefix
        suffix_n = 2
        while common_prefix in used_prefixes:
            common_prefix = f"{base_prefix}_{suffix_n}"
            suffix_n += 1
        used_prefixes.add(common_prefix)
        pair_keys[(file1, file2)] = common_prefix
        file_stats[common_prefix] = {"kept_pairs": 0, "kept_R1_singletons": 0, "kept_R2_singletons": 0, "rejected_R1": 0, "rejected_R2": 0}
        if not stdout:
            if interleaved_out:
                key = f"{common_prefix}"
                file_writing_handles[key] = open_fastq_writer(
                    key,
                    output_dir,
                    gzip_output=parameters["gzip_output"]
                    )
                if write_rejected:
                    for suffix in ["_rejected"]:
                        key = f"{common_prefix}{suffix}"
                        file_writing_handles[key] = open_fastq_writer(
                            key,
                            output_dir,
                            gzip_output=parameters["gzip_output"]
                            )
            else:
                suffixes = ["_R1_paired", "_R2_paired"]
                if not discard_singletons:
                    suffixes += ["_R1_unpaired", "_R2_unpaired"]
                for suffix in suffixes:
                    key = f"{base_prefix}{suffix}"
                    file_writing_handles[key] = open_fastq_writer(
                        key,
                        output_dir,
                        gzip_output=parameters["gzip_output"]
                        )
                if write_rejected:
                    for suffix in ["_R1_rejected", "_R2_rejected"]:
                        key = f"{base_prefix}{suffix}"
                        file_writing_handles[key] = open_fastq_writer(
                            key,
                            output_dir,
                            gzip_output=parameters["gzip_output"]
                            )

    def unified_chunk_streamer():
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

    backpressure = threading.Semaphore(threads * 5)
    def bounded_chunk_stream():
        for item in chunk_stream:
            backpressure.acquire()
            yield item

    try:
        with mp.Pool(threads, initializer=worker_initilizer, initargs=(parameters,)) as pool:
            submit = pool.imap if parameters["ordered_output"] else pool.imap_unordered
            for result in submit(unified_worker, bounded_chunk_stream(), chunksize=1):
                backpressure.release()
                if result[0] == "unpaired":
                    _, filepath, chunk_results, kept, rejected, rejected_reads = result
                    if parameters["stdout"]:
                        if chunk_results:
                            sys.stdout.write(chunk_results)
                    else:
                        if chunk_results:
                            file_writing_handles[filepath].write(chunk_results)
                        if write_rejected and rejected_reads:
                            file_writing_handles[f"{filepath}_rejected"].write(rejected_reads)
                
                    file_stats[filepath]["kept"] += kept
                    file_stats[filepath]["rejected"] += rejected
                    tracker.update(kept + rejected)
                elif result[0] == "paired":
                    _, file1, file2, paired_out_1, paired_out_2, R1_singles_out, R2_singles_out, num_paired, num_R1_singles, num_R2_singles, rejected_1, rejected_2, rejected_R1, rejected_R2 = result
                    common_prefix = pair_keys[(file1, file2)]
                    if parameters["stdout"]:
                        if paired_out_1:
                            sys.stdout.buffer.write(paired_out_1)
                    else:
                        if interleaved_out:
                            if paired_out_1:
                                file_writing_handles[common_prefix].write(paired_out_1)
                            if write_rejected:
                                if rejected_R1:
                                    file_writing_handles[f"{common_prefix}_rejected"].write(rejected_R1)
                        else:
                            writes = [
                                (f"{common_prefix}_R1_paired", paired_out_1),
                                (f"{common_prefix}_R2_paired", paired_out_2),
                            ]
                            if not discard_singletons:
                                writes.extend([
                                    (f"{common_prefix}_R1_unpaired", R1_singles_out),
                                    (f"{common_prefix}_R2_unpaired", R2_singles_out)
                                ])
                            if write_rejected:
                                writes.extend([
                                    (f"{common_prefix}_R1_rejected", rejected_R1),
                                    (f"{common_prefix}_R2_rejected", rejected_R2)
                                ])
                            for handle_key, records in writes:
                                if records:
                                    file_writing_handles[handle_key].write(records)
                    file_stats[common_prefix]["kept_pairs"] += num_paired
                    file_stats[common_prefix]["kept_R1_singletons"] += num_R1_singles
                    file_stats[common_prefix]["kept_R2_singletons"] += num_R2_singles
                    file_stats[common_prefix]["rejected_R1"] += rejected_1
                    file_stats[common_prefix]["rejected_R2"] += rejected_2
                    tracker.update(num_paired * 2 + num_R1_singles + num_R2_singles + rejected_1 + rejected_2)
    finally:
        tracker.close()
        ACTIVE_PROGRESS_TRACKER = None
        for handle in file_writing_handles.values():
                if not handle.closed:
                    handle.close()
    return file_stats

##### Input handling #####
class CleanHelpFormatter(argparse.HelpFormatter):
    """Custom argparse HelpFormatter that cleans up comma spacing, adjusts
    the starting column position of help explanations, removes the
    empty line between a description and the options that follow it,
    suppresses empty metavar placeholder artifacts (like '[ ...]'),
    and allows manual paragraph breaks (via '\\n\\n') in descriptions.
    """
    def __init__(
        self, prog, indent_increment=2, max_help_position=50, width=None
    ):
        super().__init__(
            prog,
            indent_increment=indent_increment,
            max_help_position=max_help_position,
            width=width,
        )

    def _format_args(self, action, default_metavar):
        """Suppress argument placeholder formatting (e.g., '[ ...]') when

        metavar is empty.
        """
        if action.metavar == "" or action.metavar == ("",):
            return ""
        return super()._format_args(action, default_metavar)

    def _format_action_invocation(self, action):
        """Fixes ' --cut-both , -cb' -> ' --cut-both, -cb'."""
        invocation = super()._format_action_invocation(action)
        return invocation.replace(" ,", ",")

    def _fill_text(self, text, width, indent):
        """Preserve manual '\\n\\n' paragraph breaks in description text,

        wrapping each paragraph individually rather than collapsing the
        whole description into one wrapped block.
        """
        paragraphs = text.split("\n\n")
        return "\n\n".join(
            super(CleanHelpFormatter, self)._fill_text(p, width, indent)
            if p.strip()
            else ""
            for p in paragraphs
        )

    def format_help(self):
        """Removes the blank line argparse inserts between a group's

        description and its first option, while leaving blank lines
        between groups (and manual '\\n\\n' description breaks) intact.
        """
        help_text = super().format_help()
        help_text = re.sub(r"\n\n(?=  -)", "\n", help_text)
        return help_text

def print_adapters():
    """
    Prints all built-in adapter groups and their sequences to stdout,
    formatted and aligned by group. Called for --list-adapters.
    """
    print("\nBuilt-in adapter groups and sequences:\n")
    all_entries = [entry for _, entries in DEFAULT_ADAPTERS for entry in entries]
    name_width = max(len(name) for name, _ in all_entries) + 2
    for group_name, entries in DEFAULT_ADAPTERS:
        print(f"  [{group_name}]")
        for name, sequence in entries:
            print(f"    {name:<{name_width}} {sequence}")
        print()

def print_full_auto_help(parser):
    """
    Prints the effective settings --full-auto/-GO applies, grouped by the
    parser's argument groups, using each argument's registered default
    (or its full-auto override, where one exists).
    """
    print(
        "\n"
        "When --full-auto/-GO is specified, only input parameters --input-files/-i, --input-paired/-ip, and --input-unpaired/-iu are respected."
        "\nIf none of these are provided, Readzor will auto-detect FASTQ files in the current working directory."
        "\nAll other user-provided parameters are ignored."
        "\n"
        "\nIn full automatic mode, the default settings are used, with the following specific changes:"
        "\n"
    )
    for group in parser._action_groups:
        rows = []
        for action in group._group_actions:
            dest = action.dest
            if dest not in FULL_AUTO_OVERRIDES:
                continue
            option = action.option_strings[0] if action.option_strings else dest
            rows.append((option, FULL_AUTO_OVERRIDES[dest]))
        if not rows:
            continue
        print(f"{group.title}:")
        for option, value in rows:
            print(f"    {option:<28} {value}")
        print()

def parse_args():
    """
    Parses Readzor's command-line arguments into a resolved parameters dict.
    
    Supports three main input modes: fully automatic operation (--full-auto,
    which autodetects files and ignores all other options), a flat list of
    FASTQ files (--files), or an explicit R1/R2 pair (--paired). The latter
    two are mutually exclusive.
    
    All trimming, quality, and runtime parameters default to None if not
    specified, signaling to downstream code that an adaptive/platform-
    dependent default should be resolved later (e.g. based on read length
    or Slurm detection) rather than being hardcoded here.
    
    Returns:
        dict: A parameters dictionary, keyed one-for-one with the resolved
            CLI options (input files/pairing, all trimming/filtering module
            flags and their thresholds, output/gzip/threading/chunking
            settings, and run-mode flags such as `testrun` and `full_auto`).
            See the "--- Store parameters ---" block below this docstring
            for the authoritative, up-to-date list of keys.

    Side Effects:
        If --version, --list-adapters, or --help is passed, prints the
        corresponding output and exits the program immediately (via
        `sys.exit()`). If neither --input-files, --input-paired,
        --input-unpaired, --input-interleaved, piped stdin, nor --full-auto
        provides input, calls `parser.error(...)`, which prints a usage
        message to stderr and exits with a non-zero status.
    """
    setup_logging()

    parser = argparse.ArgumentParser(
        prog="readzor",
        description="Readzor: a modular FASTQ quality trimming pipeline.\n\n"
                     "All modules are off by default. To use a module, specify a "
                     "module flag. Further specifications with module settings possible.",
        formatter_class=CleanHelpFormatter,
        add_help=False
    )

    general_group = parser.add_argument_group("General settings")
    general_group.add_argument(
        "--help", "-h", action = "store_true", default = False,
        help='[FLAG] Show this help message and exit. Combine with --full-auto/-GO for more information on fully automatic mode.'

    )
    general_group.add_argument(
        "--version", "-v", action = "store_true", default = False,
        help="[FLAG] Show Readzor version and exit."
    )
    general_group.add_argument(
        "--list-adapters", action = "store_true", default = False,
        help="[FLAG] Show all built-in adapter sequences and exit."
    )
    general_group.add_argument(
        "--full-auto", "-GO", action = "store_true", default = False,
        help="[FLAG] Run Readzor in fully automatic mode. Combine with --help/-h for more information on fully automatic mode."
    )
    general_group.add_argument(
    "--progress", action="store_true", default=False,
    help="[FLAG] Show a live progress bar and estimated time remaining during processing, based on estimated read counts. Default: off."
    )
    general_group.add_argument(
    "--verbose", action="store_true", default=False,
    help="[FLAG] Write verbose output to terminal, in addition to the log file. Default: off."
    )

    input_group = parser.add_argument_group(
        "Input options",
        "Specify input FASTQ files using any combination of --input-files, --input-paired, and --input-unpaired. "
        "Lists with any combination of regular (fastq/fq), and gzipped (fastq.gz/fq.gz) files accepted."
        "If no input flags are provided, data will be read directly from standard input (stdin)."
    )
    input_group.add_argument(
        "--input-files", "-i", nargs='+', default = None, metavar = "",
        help="FASTQ files of unspecified pairing (or data from stdin). Paired and unpaired files will be auto-detected."
    )
    input_group.add_argument(
        "--input-paired", "-ip", nargs='+', default=None, metavar = "",
        help="Paired-end FASTQ files, given as one or more R1/R2 pairs, e.g. --input-paired sample1_R1 sample1_R2 sample2_R1 sample2_R2"
    )
    input_group.add_argument(
        "--input-interleaved", "-ii", nargs='+', default = None, metavar = "",
        help="Interleaved FASTQ files. Note: these files will be split in forward and reverse reads."
    )
    input_group.add_argument(
        "--input-unpaired", "-iu", nargs='+', default = None, metavar = "",
        help="Unpaired FASTQ files."
    )

    output_group = parser.add_argument_group("Output options")
    output_group.add_argument(
        "--gzip", action="store_true", default = False,
        help="[FLAG] Compress filtered FASTQ files in gzip format. Default: off."
    )
    output_group.add_argument(
        "--gzip-level", type=int, default = 1, metavar = "", choices=range(0, 4),
        help="Set gzip compression level. Higher compression decreases processing speed. Possible values: 0-3. Default: 1."
    )
    output_group.add_argument(
        "--output", "-o", type=str, default = None, metavar = "",
        help="Path to directory in which the timestamped results folder will be created. Default: current working directory."
    )
    output_group.add_argument(
        "--stdout", action="store_true", default=False,
        help="Stream resulting FASTQ reads to stdout. Forces --interleaved-out for paired and interleaved files. Overrides --verbose and --progress. Overridden to 'off' by --gzip. Note: all files will be streamed on end, without any seperators."
    )
    output_group.add_argument(
        "--interleaved-out", action="store_true", default=False,
        help="Interleave surviving FASTQ reads of paired and interleaved input files, resulting in one output file."
    )    
    output_group.add_argument(
        "--write-rejected", action="store_true", default=False,
        help="Write rejected reads to file. Either one (for unpaired and when --interleaved-out is set), or two (for forward and reverse reads) are produced. Overridden to 'off' when --stdout is set. Default: off."
    )  
    output_group.add_argument(
        "--discard-singletons", action="store_true", default=False,
        help="Discard singletons. For paired and interleaved reads, single surviving reads will be discarded instead of written to a seperate file. No effect on unpaired read filtering. Default: off."
    )  

    general_quality_group = parser.add_argument_group("General output filter options")
    general_quality_group.add_argument(
        "--min-length-input", type=int, default = None, metavar = "",
        help="Minimum length for input read. Default: off."
    )
    general_quality_group.add_argument(
        "--max-length-input", type=int, default = None, metavar = "",
        help="Maximum length for input read. Default: off."
    )
    general_quality_group.add_argument(
        "--min-length-output", type=int, default = None, metavar = "", 
        help="Minimum length of output read. Default: off."
    )
    general_quality_group.add_argument(
        "--max-length-output", type = int, default = None, metavar = "", 
        help="Maximum length of output read. Default: off."
    )
    general_quality_group.add_argument(
        "--min-length-output-perc", type=int, default = None, metavar = "", 
        help="Minimum length of output read as percentage of its input read. Overridden by --min-length-output. Default: off."
    )
    general_quality_group.add_argument(
        "--max-length-output-perc", type = int, default = None, metavar = "", 
        help="Maximum length of output read as percentage of its input read. Overridden by --max-length-output. Default: off."
    )
    general_quality_group.add_argument(
        "--min-average-qual-pre", type=int, default = 0, metavar = "", choices=range(0, 128),
        help="Minimum average quality of input read. Default: 0."
    )
    general_quality_group.add_argument(
        "--min-average-qual-post", type=int, default = 0, metavar = "", choices=range(0, 128),
        help="Minimum average quality of output read. Default: 0."
    )
    general_quality_group.add_argument(
        "--n-filter", action="store_true", default=False,
        help="[FLAG] Reject raw reads containing N bases anywhere in read. Default: off."
    )

    trim_ends_group = parser.add_argument_group("Set-length end trimming",
                                                "Trim a set number of bases of the ends of each read, independent of sequence or quality."
                                                )
    trim_ends_group.add_argument(
        "--cut-flag", "-cf",action="store_true", default = False,
        help="[FLAG] Turn on set-length end trimming module. Default: off."
    )
    trim_ends_group.add_argument(
        "--cut-start", "-cs", type=int, default = 0, metavar="", 
        help="Number of bases to trim from the start of the read. Default: 0."
    )
    trim_ends_group.add_argument(
        "--cut-end", "-ce", type=int, default = 0, metavar="", 
        help="Number of bases to trim from the end of the read. Default: 0."
    )
    trim_ends_group.add_argument(
        "--cut-both", "-cb", type=int, default = 0, metavar="", 
        help="Number of bases to trim from both ends of the read. Overwritten by --cut-start and --cut-end. Default: 0."
    )

    quality_ends_group = parser.add_argument_group("Quality-dependent end trimming",
                                                   "Trim the ends of each read, dependent on quality. Ends of reads will be trimmed up to first position that fulfills quality requirement.")
    quality_ends_group.add_argument(
        "--endqual-filter-flag", "-ef", action="store_true", default = False,
        help="[FLAG] Turn on quality-dependent end trimming. Default: off."
    )
    quality_ends_group.add_argument(
        "--endqual-min-start", "-ems", type = int, default = None, metavar="", choices=range(0, 127),
        help="Specific phred score threshold for the start of the read. Default: 25."
    )
    quality_ends_group.add_argument(
        "--endqual-min-end", "-eme", type=int, default = None, metavar="",choices=range(0, 127),
        help="Specific phred score threshold for the end of the read.  Default: 25."
    )
    quality_ends_group.add_argument(
        "--endqual-min-both", "-emb", type=int, default = 25, metavar="", choices=range(0, 127),
        help="Phred score threshold for the quality trimming of read ends. Overwritten by --endqual-min-start and --endqual-min-end. Default: 25."
    )

    n_ends_group = parser.add_argument_group("N nucleotide end-trimming",
                                             "Trim the ends of each read for N bases. Redundant when --n-filter is set.")
    n_ends_group.add_argument(
        "--n-trimming-flag", "-ntf", action="store_true", default = False,
        help="[FLAG] Turn on N nucleotide end-trimming. Default: off."
    )

    sliding_window_group = parser.add_argument_group("Sliding window quality trimming",
                                                     "Trim the reads for quality based on a sliding window of size X, moved with stepsize Y. Longest portion survives in case of mid-read quality dropoff.")
    sliding_window_group.add_argument(
        "--slider-filter-flag", "-sf", action="store_true", default = False,
        help="[FLAG] Turn on sliding window quality trimming module. Default: off."
    )
    sliding_window_group.add_argument(
        "--slider-window", "-sw", type = int, default = 5, metavar="",
        help="Window size over which average quality is calculated. Default: 5."
    )
    sliding_window_group.add_argument(
        "--slider-quality", "-sq", type = int, default = 20, metavar="",
        help="Minimum average quality in sliding window. Default: 20."
    )
    sliding_window_group.add_argument(
        "--slider-step", "-ss", type = int, default = 1, metavar="",
        help="Sliding window step size. Default: 1."
    )

    homopolymer_nucleotide_trimming = parser.add_argument_group("Homopolymer nucleotide trimming",
                                                                "Illumina NovaSeq, NextSeq, and MiniSeq use a two-color chemistry, in which guanine bases are unlabeled. In event of short fragments, this can result in homolopolymer G calls at the end of reads.")
    homopolymer_nucleotide_trimming.add_argument(
        "--poly-filter-flag", "-pf", action="store_true", default = False,
        help="[FLAG] Turn on homopolymer read-end trimming module. Default: off."
    )
    homopolymer_nucleotide_trimming.add_argument(
        "--poly-bases-start", "-pbs", type = str, default = None, metavar = "",
        help='Base(s) to check for a homopolymer run at start of read. Comma-separated bases are checked independently. Default: none.'
    )
    homopolymer_nucleotide_trimming.add_argument(
        "--poly-bases-end", "-pbe", type = str, default = "G", metavar = "",
        help='Base(s) to check for a homopolymer run at end of read. Comma-separated bases are checked independently. Default: "G".'
    )
    homopolymer_nucleotide_trimming.add_argument(
        "--poly-bases-both", "-pbb", type = str, default = None, metavar = "",
        help='Base(s) to check for a homopolymer run at both read ends. Comma-separated bases are checked independently. Overwritten by poly_bases_start and poly_bases_end. Default: none.'
    )
    homopolymer_nucleotide_trimming.add_argument(
        "--poly-length-start", "-pls", type = int, default = 10, metavar = "",
        help="Minimum length of homopolymer run at start of read required to trigger trimming. Default: 10."
    )
    homopolymer_nucleotide_trimming.add_argument(
        "--poly-length-end", "-ple", type = int, default = 10, metavar = "",
        help="Minimum length of homopolymer run at end of read required to trigger trimming. Default: 10."
    )
    homopolymer_nucleotide_trimming.add_argument(
        "--poly-length-both", "-plb", type = int, default = 0, metavar = "",
        help="Minimum length of homopolymer run at start and end of read required to trigger trimming. Default: 0."
    )

    adapter_trimming = parser.add_argument_group("Adapter trimming",
                                                 "Trim reads for Illumina adapter sequences. Standard sequences included are TruSeq3 universal and index adapters, and Nextera adapters. Only exactly matching sequences are trimmed. Adapter trimming is performed independent of quality.")
    adapter_trimming.add_argument(
        "--adapter-filter-flag", "-af", action="store_true", default = False,
        help='[FLAG] Turn on adapter trimming module. Default: off.'
    )
    adapter_trimming.add_argument(
        "--adapter-group", "-ag", nargs = "+", default = "Nextera", choices = [name for name, _ in DEFAULT_ADAPTERS], metavar="",
        help='Specify the group(s) of adapters to be used. Ignored if --adapter-fasta-excl is set. Choices: Illumina_RNA, Nextera, TruSeq2, TruSeq3, TruSeq_small_RNA. Default: Nextera.'
    )
    adapter_trimming.add_argument(
        "--adapter-mismatch", "-am", type = int, default = 0, metavar="",
        help="Number of mismatches allowed in adapter finding. Note: setting mismatches to > 0, infers significant processing constrains. Default: 0."
    )
    adapter_trimming.add_argument(
        "--adapter-fasta-add", "-ad", type = str, default = None, metavar="",
        help="Fasta file with adapter sequences to trim for, in addition to predefined sequences."
    )
    adapter_trimming.add_argument(
        "--adapter-fasta-excl", "-ax", type = str, default = None, metavar="",
        help="Fasta file with adapter sequences to trim for, excluding predefined and additional sequences specified."
    )

    low_complexity_group = parser.add_argument_group("Low complexity filtering",
                                                     "Detect complexity of reads using kmer-based nucleotide frequencies. Low complex reads discarded entirely.")
    low_complexity_group.add_argument(
        "--kmer-filter-flag", "-kf", action="store_true", default=False,
        help="[FLAG] Turn on the kmer-based complexity filtering module. Default: off."
    )
    low_complexity_group.add_argument(
        "--kmer-size", "-ks", default = 4, metavar="",
        help="Kmer length for kmer-based complexity filtering. Comma-separated values are checked independently. Maximum value: 21. Default: 4."
    )
    low_complexity_group.add_argument(
        "--kmer-cutoff", "-kc", type = int, default = 50, metavar="",
        help="Minimum percentage of unique k-mers (relative to the maximum possible) required to pass the complexity filter. Higher values are stricter. Default: 50."
    )
    mgi_convert_group = parser.add_argument_group("MGI header conversion",
                                                  "Convert read header from MGI (BGI) format to Illumina format. Original header will be stored in the placeholder line. Conversion is necessary for downstream analysis with tools such as samtools")
    mgi_convert_group.add_argument(
        "--mgi-convert-flag", "-mf", action="store_true", default=False,
        help="[FLAG] Turn on the MGI-to-Illumina header conversion module. Default: off."
    )
    mgi_convert_group.add_argument(
        "--mgi-bc5", "-m5", type = str, default = "PLACEHOLDERi5", metavar="",
        help="Input a i5 barcode for Illumina header conversion. Default: 'PLACEHOLDERi5'."
    )
    mgi_convert_group.add_argument(
        "--mgi-bc7", "-m7", type = str, default = "PLACEHOLDERi7", metavar="",
        help="Input a i7 barcode for Illumina header conversion. Default: 'PLACEHOLDERi7'."
    )
    mgi_convert_group.add_argument(
        "--mgi-instrument", "-mi", type = str, default = "PLACEHOLDERinstrument", metavar="",
        help="Instrument name for Illumina header conversion. Default: 'PLACEHOLDERinstrument'."
    )
    mgi_convert_group.add_argument(
        "--mgi-run", "-mr", type = str, default = "PLACEHOLDERrun", metavar="",
        help="Run ID for Illumina header conversion. Default: 'PLACEHOLDERrun'."
    )

    advanced_group = parser.add_argument_group("Advanced options",
                                               "Further options that can be specified to alter the behaviour of Readzor.")
    advanced_group.add_argument(
        "--threads", "-t", type = int, default = None, metavar="",
        help="Number of threads to use. Default: platform-dependent through auto-detection: detection of assigned CPUs on HPC clusters, all-1 otherwise. Fallback: 1."
    )
    advanced_group.add_argument(
        "--reads-for-phred-offset", type = int, default = 500, metavar="",
        help="Number of reads to sample per file for detection of Phred quality encoding offset. Default: 500."
    )
    advanced_group.add_argument(
        "--chunk-size", type=int, default = None, metavar="",
        help="Number of reads per chunk sent to each worker. Default: platform-dependent (empirically set to 20,000 for HPC cluster systems, 1000 otherwise). Changing can alter processing speed."
    )
    advanced_group.add_argument(
        "--phred-offset", type = int, choices=[33, 64], default = None, metavar="",
        help="Define phred offset for all FASTQ files. When set, per-file auto-detection will not be performed. Possible values: 33, 64. Default: off (auto-detection per file)."
    )
    advanced_group.add_argument(
        "--testrun", action="store_true", default=False,
        help="[FLAG] Perform a test run according to specified settings. Implies --verbose. Default: off."
    )
    advanced_group.add_argument(
        "--ordered-output", action="store_true", default=False,
        help="[FLAG] Force writing output reads in the same order they appear in the input file. Default: off."
    )
    advanced_group.add_argument(
        "--phred-out", type=int, choices=[33, 64], default=None, metavar="",
        help="Convert Phred encoding from 33 to 64, and vice versa. Possible values: 33, 64. Default: off."
    )
    
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()
    setup_logging(verbose=args.verbose)

    if args.help and args.full_auto:
        print_full_auto_help(parser)
        sys.exit()
    if args.help:
        parser.print_help()
        sys.exit()

    if args.version:
        print(f"Readzor version: {VERSION}")
        sys.exit()

    if args.list_adapters:
        print_adapters()
        sys.exit()

    has_inputs = bool(
        args.input_files or 
        args.input_paired or 
        args.input_unpaired or 
        args.input_interleaved
    )
    
    if not has_inputs:
        if not sys.stdin.isatty():
            args.input_files = [resolve_stdin_input()]
        elif args.full_auto:
            cwd = os.getcwd()
            pattern = re.compile(r'\.(fastq|fq)(\.gz|\.gzip)?$', re.IGNORECASE)
            args.input_files = [
                os.path.join(cwd, f)
                for f in os.listdir(cwd)
                if os.path.isfile(os.path.join(cwd, f)) and pattern.search(f)
            ]
            if not args.input_files:
                parser.error(
                    f"--full-auto was set but no FASTQ files were detected in {cwd}."
                )
        else:
            parser.error(
                "You must specify input files (--input-files, --input-paired, etc.), "
                "pipe data via stdin, or use --full-auto."
            )
    
    if args.full_auto:
        for action in parser._actions:
            dest = action.dest
            if dest == "help" or dest in FULL_AUTO_PRESERVED_DESTS:
                continue
            reset_value = FULL_AUTO_OVERRIDES.get(dest, action.default)
            setattr(args, dest, reset_value)
            
        if not args.verbose:
            print("[WARNING] --full-auto/-GO specified; ignoring all other input parameters (except input file parameters).")

    # --- Store parameters ---
    parameters = {}
    parameters["testrun"] = args.testrun
    parameters["full_auto"] = args.full_auto
    parameters["unspecified_files"] = args.input_files
    parameters["unpaired_files"] = args.input_unpaired
    parameters["paired_files"] = group_paired_input_into_pairs(files = args.input_paired, parser = parser)
    parameters["interleaved_files"] = args.input_interleaved
    parameters["min_quality_both"] = args.endqual_min_both
    parameters["endqual_min_start"] = args.endqual_min_start
    parameters["endqual_min_end"] = args.endqual_min_end
    parameters["minimum_average_qual_post"] = args.min_average_qual_post
    parameters["minimum_average_qual_pre"] = args.min_average_qual_pre    
    parameters["cut_start"] = args.cut_start
    parameters["cut_end"] = args.cut_end
    parameters["cut_both"] = args.cut_both
    parameters["n_trimming_flag"] = args.n_trimming_flag
    parameters["slider_window"] = args.slider_window
    parameters["slider_step"] = args.slider_step
    parameters["slider_quality"] = args.slider_quality
    parameters["gzip_output"] = args.gzip
    parameters["gzip_level"] = args.gzip_level
    parameters["reads_for_phred_offset"] = args.reads_for_phred_offset
    parameters["adapter_filter_flag"] = args.adapter_filter_flag
    parameters["adapter_fasta_add"] = args.adapter_fasta_add
    parameters["adapter_fasta_excl"] = args.adapter_fasta_excl
    parameters["n_filter"] = args.n_filter
    parameters["phred_offset"] = args.phred_offset
    parameters["threads"] = args.threads
    parameters["chunk_size"] = args.chunk_size
    parameters["output_dir"] = args.output or os.getcwd()
    parameters["mgi_convert_flag"] = args.mgi_convert_flag
    parameters["mgi_bc5"] = args.mgi_bc5
    parameters["mgi_bc7"] = args.mgi_bc7
    parameters["mgi_instrument"] = args.mgi_instrument
    parameters["mgi_run"] = args.mgi_run
    parameters["kmer_filter_flag"] = args.kmer_filter_flag
    parameters["kmer_size"] = args.kmer_size
    parameters["kmer_cutoff"] = args.kmer_cutoff
    parameters["allow_n_kmer"] = not args.n_filter
    parameters["cut_flag"] = args.cut_flag
    parameters["endqual_filter_flag"] = args.endqual_filter_flag
    parameters["slider_filter_flag"] = args.slider_filter_flag
    parameters["poly_filter_flag"] = args.poly_filter_flag
    parameters["poly_bases_start"] = args.poly_bases_start
    parameters["poly_bases_end"] = args.poly_bases_end
    parameters["poly_bases_both"] = args.poly_bases_both
    parameters["poly_length_start"] = args.poly_length_start
    parameters["poly_length_end"] = args.poly_length_end
    parameters["poly_length_both"] = args.poly_length_both
    parameters["show_progress"] = args.progress
    parameters["verbose"] = args.verbose
    parameters["adapter_mismatch"] = args.adapter_mismatch
    parameters["ordered_output"] = args.ordered_output
    parameters["phred_out"] = args.phred_out
    parameters["min_length_input"] = args.min_length_input
    parameters["max_length_input"] = args.max_length_input
    parameters["min_length_output"] = args.min_length_output
    parameters["max_length_output"] = args.max_length_output
    parameters["min_length_output_perc"] = args.min_length_output_perc
    parameters["max_length_output_perc"] = args.max_length_output_perc
    parameters["stdout"] = args.stdout
    parameters["interleaved_out"] = args.interleaved_out
    parameters["write_rejected"] = args.write_rejected
    parameters["adapter_group"] = args.adapter_group
    parameters["discard_singles"] = args.discard_singles
    
    if parameters["gzip_output"]:
        parameters["stdout"] = False
        
    if parameters["stdout"]:
        parameters["verbose"] = False
        parameters["progress"] = False
        parameters["write_rejected"] = False
    
    if parameters["n_filter"]:
        parameters["n_trimming_flag"] = False

    if parameters.get("adapter_filter_flag"):
        if parameters.get("adapter_fasta_excl"):
            raw_adapters = load_adapters_from_fasta(parameters["adapter_fasta_excl"])
            if raw_adapters == []:
                raise ValueError(f"No sequences in file '{parameters['adapter_fasta_excl']}' detected.")
        else:
            selected_groups = parameters.get("adapter_group")
            if selected_groups:
                selected_groups = set(selected_groups)
                flat_default_adapters = [entry for name, entries in DEFAULT_ADAPTERS if name in selected_groups for entry in entries]
            else:
                flat_default_adapters = [entry for _, entries in DEFAULT_ADAPTERS for entry in entries]
            if parameters.get("adapter_fasta_add"):
                raw_adapters = flat_default_adapters + load_adapters_from_fasta(parameters["adapter_fasta_add"])
            else:
                raw_adapters = flat_default_adapters
        parameters["adapter_sequences"] = list({seq.encode('utf-8') for _, seq in raw_adapters})
    else:
        parameters["adapter_sequences"] = []

    parameters["threads"] = worker_determination(parameters["threads"])
    parameters["chunk_size"] = chunk_size_setter(parameters["chunk_size"])
    
    return parameters

##### Wrap up functions #####
def write_summary_and_statistics(summary_results, parameters, output_dir):
    """
    Write summary statistics and parameters to output text files.
 
    The function creates two files in the specified output directory:
    ``results_summary.txt`` containing per-file summary counts and
    ``parameters.txt`` containing the parameter key-value pairs used
    for the analysis.
    
    Args:
        summary_results (dict): Dictionary mapping file names to dictionaries
            of summary statistics. Entries may contain either ``kept`` and
            ``rejected`` counts, or ``kept_pairs``, ``kept_R1_singletons``,
            ``kept_R2_singletons``, ``rejected_R1``, and ``rejected_R2``
            counts for paired entries (R1/R2 rejects tracked separately
            since the mates are filtered and split independently).
        parameters (dict): Dictionary of parameter names and their values to
            write to the parameters output file.
        output_dir (str): Path to the directory where the output files will
            be created.
    
    Returns:
        None
    """
    paired_data = []
    unpaired_data = []
    for file_path, counts in summary_results.items():
        filename = os.path.basename(file_path)
        if "kept" in counts:
            unpaired_data.append((filename, counts["kept"], counts["rejected"]))
        else:
            paired_data.append((
                filename,
                counts["kept_pairs"],
                counts["kept_R1_singletons"],
                counts["kept_R2_singletons"],
                counts["rejected_R1"],
                counts["rejected_R2"],
            ))
    summary_path = os.path.join(output_dir, "results_summary.txt")
    with open(summary_path, "w", encoding="utf-8") as f:
        if paired_data:
            f.write("[Paired Reads]\n")
            f.write("Pair with common prefix\tKept_Pairs\tKept_R1_Singletons\tKept_R2_Singletons\tRejected_R1\tRejected_R2\n")
            for item in paired_data:
                f.write(f"{item[0]}\t{item[1]}\t{item[2]}\t{item[3]}\t{item[4]}\t{item[5]}\n")
            f.write("\n")
        if unpaired_data:
            f.write("[Unpaired Reads]\n")
            f.write("Filename\tKept\tRejected\n")
            for item in unpaired_data:
                f.write(f"{item[0]}\t{item[1]}\t{item[2]}\n")
  
def log_parameters(parameters):
    """
    Log all run parameters to the logger, one per line, tagged [PARAMETER]
    so they're easy to grep out of the log file. Replaces the old standalone
    parameters.txt dump.
 
    Args:
        parameters (dict): Dictionary of parameter names and their values.
 
    Returns:
        None
    """
    
    def _format_value(value):
        if isinstance(value, bytes):
            return value.decode('utf-8', errors='replace')
        if isinstance(value, (list, tuple, set)):
            return ", ".join(_format_value(v) for v in value)
        return str(value)

    logger.info("--- Start of run parameters ---")
    for key, value in sorted(parameters.items()):
        logger.info("[PARAMETER] %s: %s", key, _format_value(value))
    logger.info("---- End of run parameters ----")
        
def print_final_message(stdout):
    """
    Prints Readzor's completion message to the console and logs it to file only,
    including a citation request and a randomly selected sign-off phrase.
    Called at the end of a successful run.
    Returns:
        None
    """
    if stdout:
        return
    sign_off_messages = [
    "Please come again!",
    "Thanks for trimming with Readzor!",
    "May your reads be long and your adapters be gone!",
    "Until next time, happy analyzing!",
    "See you soon!",
    "Thanks for using Readzor!",
    "Happy analyzing!",
    "Good luck with your data!",
    "Base-ically, we're done here. See you!",
    "Readzor, signing off!"
    ]
    citation_notice = "If you find Readzor useful, please consider citing:"
    citation = (
        "Axel B. Janssen\n"
        "Readzor: A modular and user-friendly Swiss-army knife approach to short-read sequencing processing.\n"
        "2026\n"
        "https://doi.org/10.5281/zenodo.22336649"
    )
    sign_off = random.choice(sign_off_messages)

    print(f"\n{citation_notice}")
    print(f"\n{citation}")
    print(f"\n{sign_off}\n")

    file_only_logger.info(citation_notice)
    file_only_logger.info(citation)
    file_only_logger.info(sign_off)

##### Main #####
def main():
    """
    Main execution block for Readzor.

    Initializes command-line argument parsing and parameter configuration, creates the output
    directory structure, handles input stream processing and multithreaded trimming
    via the top-level orchestrator, and exports summary metrics and final parameters.
    Finishes with a sign-off message.
    """
    for method in ['fork', 'forkserver', 'spawn']:
        if method in mp.get_all_start_methods():
            try:
                mp.set_start_method(method)
                break
            except RuntimeError:
                pass
    parameters = parse_args()
    if parameters["testrun"]:
        test_run(parameters)
        return
    try:
        created_output_dir = create_folder_structure(parameters["output_dir"])
        setup_logging(output_dir = created_output_dir, verbose = parameters["verbose"], parameters = parameters)
        log_parameters(parameters)
        summary_results = input_handler(unspecified_files = parameters["unspecified_files"], unpaired_files = parameters["unpaired_files"], paired_files = parameters["paired_files"], interleaved_files = parameters["interleaved_files"], output_dir = created_output_dir, threads = parameters["threads"], chunk_size = parameters["chunk_size"], show_progress = parameters["show_progress"], stdout = parameters["stdout"], interleaved_out = parameters["interleaved_out"], discard_singletons = parameters["discard_singletons"], write_rejected = parameters["write_rejected"], parameters = parameters)
        write_summary_and_statistics(summary_results, parameters, output_dir = created_output_dir)
        logger.info("Analysis successfully completed!")
        print_final_message(stdout = parameters["stdout"])
    except KeyboardInterrupt:
        logger.info("Readzor interrupted — shutting down.")
        if not parameters["verbose"]:
            print("Interrupted by user — shutting down.", file=sys.stderr)
        sys.exit(130)
    finally:
        cleanup_stdin_temp_files()

def test_run(parameters):
    """
    Runs the full pipeline against a temporary output directory that is
    deleted on exit, to let users validate their settings without keeping
    any output. Limits each input file to one chunk per worker thread so
    the run finishes quickly. Called when --testrun is set.

    Args:
        parameters (dict): Dictionary of configuration parameters.
    """
    try:
        with tempfile.TemporaryDirectory(prefix="readzor_testrun_") as tmp_dir:
            created_output_dir = create_folder_structure(tmp_dir)
            setup_logging(output_dir = created_output_dir, verbose = True, parameters = parameters)
            log_parameters(parameters)
            summary_results = input_handler(unspecified_files = parameters["unspecified_files"], unpaired_files = parameters["unpaired_files"], paired_files = parameters["paired_files"], interleaved_files = parameters["interleaved_files"], output_dir = created_output_dir, threads = parameters["threads"], chunk_size = parameters["chunk_size"], show_progress = parameters["show_progress"], stdout = parameters["stdout"], interleaved_out = parameters["interleaved_out"],discard_singletons = parameters["discard_singletons"], parameters = parameters)
            write_summary_and_statistics(summary_results, parameters, output_dir = created_output_dir)
            logger.info("All files deleted.")
            logger.info("Testrun completed!")
    finally:
        cleanup_stdin_temp_files()
    
if __name__ == "__main__":
    main()