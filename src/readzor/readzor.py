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

from isal import igzip, igzip_threaded
import numpy as np

##### Definition of constant values #####
WORKER_PARAMETERS = None
ESTIMATED_ZIP_RATIO = {}
ESTIMATED_READ_COUNTS = {}
STDIN_TEMP_FILES = []
ESTIMATED_BYTE_PER_READ = {}
GZIP_DETECTION = {}
VERSION = "0.4.7"
PHRED_ALLOWED = bytes(range(33, 127))
DEFAULT_ADAPTERS = [
    ["TruSeq", [
        ["Read_1", "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"],  #https://support-docs.illumina.com/SHARE/AdapterSequences/Content/UDIndexes.htm
        ["Read_2", "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"]  #https://support-docs.illumina.com/SHARE/AdapterSequences/Content/UDIndexes.htm
    ]],
    ["TruSeq_small_RNA", [
        ["TruSeq_small_RNA", "TGGAATTCTCGGGTGCCAAGG"]    #https://support-docs.illumina.com/SHARE/AdapterSequences/Content/TruSeq-SmallRNA.htm
    ]],
    ["Illumina_miRNA", [
        ["TruSeq_small_RNA", "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"]     #https://support-docs.illumina.com/SHARE/AdapterSequences/Content/Illumina-miRNA-Indexes.htm
    ]],
    ["Nextera", [
        ["Nextera", "CTGTCTCTTATACACATCT"]  #https://support-docs.illumina.com/SHARE/AdapterSequences/Content/Nextera_Illumina-Sequences.htm
    ]],
    ["AmpliSeq", [
        ["AmpliSeq", "CTGTCTCTTATACACATCT"]  #https://support-docs.illumina.com/SHARE/AdapterSequences/Content/AmpliSeq-Sequences.htm
    ]],
    ["TruSeq_DNA_methylation", [
        ["Read_1", "AGATCGGAAGAGCACACGTCTGAAC"],  #https://support-docs.illumina.com/SHARE/AdapterSequences/Content/TruSeq-DNAMethyl.htm
        ["Read_2", "AGATCGGAAGAGCGTCGTGTAGGGA"],  #https://support-docs.illumina.com/SHARE/AdapterSequences/Content/TruSeq-DNAMethyl.htm
    ]],
    ["TruSeq_Ribo_profile", [
        ["TruSeq_Ribo_profile", "AGATCGGAAGAGCACACGTCT"]  #https://support-docs.illumina.com/SHARE/AdapterSequences/Content/TruSeq-RiboProfile.htm
    ]],
    ["Illumina_RNA", [
        ["Illumina_RNA", "ACTGTCTCTTATACACATCT"]    #https://support-docs.illumina.com/SHARE/AdapterSequences/Content/Nextera_Illumina-Sequences.htm
    ]]
]
FULL_AUTO_PRESERVED_DESTS = {"input_files", "input_paired", "input_unpaired", "input_interleaved", "full_auto"}
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
COMP_TABLE = bytes.maketrans(NUCL_ATCGN, b"TAGCN")
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
    printing a log record, so the two don't get interleaved on the same
    terminal line.
    """
    def emit(self, record):
        """
        Clear the active progress-bar line (if there is one), then emit
        the record as a normal StreamHandler would.
        """
        if ACTIVE_PROGRESS_TRACKER is not None:
            ACTIVE_PROGRESS_TRACKER.clear_line()
        super().emit(record)

def setup_logging(output_dir = None, verbose = False, parameters = None):
    """
    Configure the logger.

    Called from parse_args() twice: once before parsing, without arguments,
    to attach a console handler so early messages are visible, and again
    right after parsing to apply the user's --verbose choice to that same
    handler. Called a third time, from main() or test_run(), once the
    timestamped results folder exists, to attach a file handler that writes
    the full run to <output_dir>/Readzor_log.txt regardless of console
    verbosity. On that third call it also logs the Readzor version, the
    exact command line, whether full-auto and/or test-run mode is active,
    and the output folder and log file locations, and attaches the same
    file handler to ``file_only_logger`` (for messages that should only go
    to the log file). The file handler is only ever attached once.

    Args:
        output_dir (str | None): Results folder to write Readzor_log.txt
            into. If None, only the console handler is (re)configured.
        verbose (bool): If True, the console shows DEBUG-level messages and
            above. If False, the console handler is effectively silenced
            (its level is set above CRITICAL). The log file always captures
            DEBUG and above.
        parameters (dict | None): Run parameters. Required (and only
            consulted) when output_dir is given, to read the "full_auto"
            and "testrun" flags.
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

def count_reads_estimated(filepath, sample_size=50, default_gzip_ratio=4, gzip_sample_bytes=50 * 1024 * 1024):
    """
    Estimate the number of reads in a FASTQ file.

    The average number of (uncompressed) bytes per read is estimated from
    the first ``sample_size`` records: the lengths of the four lines plus
    4 for the newline characters that lazy_fastq() strips. The estimated
    uncompressed file size is then divided by this value.

    For gzip files, the uncompressed size is estimated as the file size
    times a compression ratio, measured by decompressing up to
    ``gzip_sample_bytes`` from the start of the file. If less than 1 MiB
    could be decompressed, ``default_gzip_ratio`` is used instead.

    Requires GZIP_DETECTION[filepath] to be set (see input_handler()).
    Stores the bytes-per-read estimate in ESTIMATED_BYTE_PER_READ and, for
    gzip files, the compression ratio in ESTIMATED_ZIP_RATIO, both keyed
    by ``filepath``.

    Args:
        filepath (str): Path to the FASTQ file (plain or gzip-compressed).
        sample_size (int): Maximum number of records sampled for the
            bytes-per-read estimate. Defaults to 10.
        default_gzip_ratio (float): Compression ratio used when the gzip
            sample is too small to be reliable. Defaults to 4.
        gzip_sample_bytes (int): Number of uncompressed bytes to decompress
            when measuring the gzip ratio. Defaults to 50 MiB.

    Returns:
        int: Estimated number of reads in the file, always >= 1.

    Raises:
        ValueError: If the file contains no complete FASTQ records.
    """
    blob = next(lazy_fastq_blobs(filepath, sample_size), (b"", 0))[0]
    records = list(parse_blob(blob))
    if not records:
        raise ValueError(f"No FASTQ records found in '{filepath}'; cannot estimate bytes per read.")
    bytes_per_read = len(blob) // len(records)
    ESTIMATED_BYTE_PER_READ[filepath] = bytes_per_read

    file_size = os.path.getsize(filepath)
    if GZIP_DETECTION[filepath]:
        compressed = 0
        uncompressed = 0
        with open(filepath, 'rb') as raw:
            decompressor = igzip.GzipFile(fileobj=raw)
            while uncompressed < gzip_sample_bytes:
                chunk = decompressor.read(1024 * 1024)
                if not chunk:
                    break
                uncompressed += len(chunk)
            compressed = raw.tell()
        if compressed == 0 or uncompressed < 1024 * 1024:
            ratio = default_gzip_ratio
        else:
            ratio = uncompressed / compressed
        ESTIMATED_ZIP_RATIO[filepath] = ratio
        estimated_uncompressed_size = file_size * ratio
    else:
        estimated_uncompressed_size = file_size

    return max(1, round(estimated_uncompressed_size / bytes_per_read))

class ProgressTracker:
    """
    Tracks the number of processed reads against an estimated total and
    draws a single-line progress bar, with processing rate and estimated
    time remaining, on stderr.

    The bar is only drawn when stderr is a terminal, is redrawn at most
    once every ``min_interval`` seconds, and shrinks to fit the terminal
    width (or is dropped, leaving only the statistics, on very narrow
    terminals). Because the total is an estimate, the displayed fraction
    is capped at 99.9% until close() is called.

    Args:
        total_reads (int): Estimated total number of reads (clamped to >= 1).
        bar_width (int): Maximum width of the bar in characters. Defaults to 50.
        min_interval (float): Minimum number of seconds between redraws.
            Defaults to 0.2.
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
        """
        Format a duration for display.

        Args:
            seconds_val (float): Duration in seconds. Infinite or negative
                values are shown as "?".
            concise (bool): If True, use a compact form ("1h 02m 03s",
                "2m 05s", "7s"); otherwise a spelled-out form
                ("1 hour, 2 minutes, 3 seconds").

        Returns:
            str: The formatted duration.
        """
        if seconds_val == float('inf') or seconds_val < 0 or seconds_val is None:
            return "?"
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
        """
        Add ``n`` processed reads and redraw the bar if at least
        ``min_interval`` seconds have passed since the last redraw.

        Args:
            n (int): Number of reads processed since the previous call.
        """
        self.done += n
        now = time.time()
        if now - self._last_render >= self.min_interval:
            self._render()
            self._last_render = now

    def _render(self):
        """
        Draw the current progress line on stderr, overwriting the previous
        one. Does nothing if stderr is not a terminal.
        """
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
        Does nothing if stderr is not a terminal.
        """
        if not sys.stderr.isatty():
            return
        sys.stderr.write("\r\x1b[2K")
        sys.stderr.flush()

    def close(self, interrupted=False):
        """
        Draw the final line (progress, total reads processed, average rate
        and total time) and end it with a newline. When stderr is not a
        terminal, only the statistics are written, as plain text.

        On a terminal, the rest of the line is erased after the text, so
        anything the terminal echoed onto it (such as "^C" after Ctrl+C) is
        removed.

        Args:
            interrupted (bool): If True, the run did not finish (Ctrl+C or
                an error). The line then says where it stopped ("Stopped at
                6.1%") and the bar is only filled up to that point, instead
                of claiming 100%. As the total is an estimate, the
                percentage is capped at 99.9%. Defaults to False.
        """
        elapsed = time.time() - self._start_time
        rate = self.done / elapsed if elapsed > 0 else 0
        time_str = self._format_duration(elapsed, concise=False)
        if interrupted:
            frac = min(self.done / self.total, 0.999)
            progress = f"Stopped at {frac*100:.1f}%"
        else:
            frac = 1.0
            progress = "100%"
        stats = (
            f" {progress} ({self.done:,} reads analyzed). "
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
            filled = int(effective_bar_width * frac)
            progressbar = "#" * filled + "-" * (effective_bar_width - filled)
            line = f"\r[{progressbar}]{stats}"
        else:
            line = f"\r{stats.strip()}"
        line = line[: term_width - 1].ljust(term_width - 1)
        sys.stderr.write(line + "\x1b[K\n")
        sys.stderr.flush()

##### Helper functions #####

def resolve_stdin_input(parser=None):
    """
    Spool FASTQ data piped into stdin to a temporary file.

    The rest of the pipeline needs a real, re-readable path (inputs are
    read several times: for pairing, Phred detection and processing), so
    stdin is copied to a temporary "readzor_stdin_*.fastq" file. Its path
    is recorded in STDIN_TEMP_FILES so cleanup_stdin_temp_files() can
    remove it later. Gzip-compressed data on stdin also works, since
    compression is detected from the file contents.

    Args:
        parser (argparse.ArgumentParser | None): If given, used to report
            an empty stdin through parser.error().

    Returns:
        str: Path to the temporary file holding the stdin data.

    Raises:
        SystemExit: If stdin is a terminal (nothing is piped in), or if
            stdin was empty (through parser.error() when a parser is given).
    """
    if sys.stdin.isatty():
        raise SystemExit("Error: No input file specified and no piped data detected on stdin.")

    fd, tmp_path = tempfile.mkstemp(suffix=".fastq", prefix="readzor_stdin_")
    with os.fdopen(fd, "wb") as tmp_file:
        shutil.copyfileobj(sys.stdin.buffer, tmp_file, length=10 * 1024 * 1024)
    STDIN_TEMP_FILES.append(tmp_path)

    if os.path.getsize(tmp_path) == 0:
        cleanup_stdin_temp_files()
        message = ("No input files were given and stdin was empty. Specify input with --input-files/--input-paired/--input-unpaired/--input-interleaved, pipe FASTQ data into Readzor, or use --full-auto.")
        if parser is not None:
            parser.error(message)
        raise SystemExit(f"Error: {message}")

    return tmp_path

def cleanup_stdin_temp_files():
    """
    Delete all temporary files created by resolve_stdin_input() and clear
    STDIN_TEMP_FILES. Files that no longer exist are ignored.
    """
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
    The order within each pair is taken as given; file contents are not
    inspected.
    
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
    pairs = []
    for i in range(0, len(files), 2):
        pairs.append(tuple(files[i:i+2]))
    return pairs

def create_folder_structure(output_dir):
    """
    Create a new, uniquely named results folder inside output_dir.

    The folder is named Readzor_results_<YYYY-MM-DD_HH-MM-SS>. If that name is
    already taken (e.g. several runs started in the same second with the same -o,
    as with SLURM array jobs), a numeric suffix is added: _2, _3, ...
    Each run therefore always gets its own folder, log and summary.

    Args:
        output_dir (str): Parent directory. Created if it does not exist.

    Returns:
        str: Absolute path to the newly created results folder.
    """
    os.makedirs(output_dir, exist_ok=True)
    base = os.path.join(output_dir, "Readzor_results_" + datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S"))
    created_output_dir, n = base, 2
    while True:
        try:
            os.mkdir(created_output_dir)
            return os.path.abspath(created_output_dir)
        except FileExistsError:
            created_output_dir = f"{base}_{n}"
            n += 1

def worker_determination(threads=None):
    """
    Determine the number of worker processes for the multiprocessing pool.

    The worker count is determined in the following order of precedence:
    1. Explicit request: ``threads``, if it is an int.
    2. Scheduler environment: the first of SLURM_CPUS_PER_TASK, PBS_NCPUS,
       NCPUS, NSLOTS or LSB_DJOB_NUMPROC that is set to a valid integer.
    3. Local default: the number of available CPUs minus one.

    The result is always at least 1.

    Args:
        threads (int | None): User-requested number of workers, or None to
            auto-detect. Values that are not an int are ignored.
            Defaults to None.

    Returns:
        int: The number of worker processes to spawn, always >= 1.
    """
    if isinstance(threads, int):
        return max(1, threads)
    for var in ("SLURM_CPUS_PER_TASK", "PBS_NCPUS", "NCPUS", "NSLOTS", "LSB_DJOB_NUMPROC"):
        val = os.environ.get(var)
        if val:
            try:
                return max(1, int(val))
            except ValueError:
                pass
    return max(1, (os.cpu_count() or 1) - 1)

def common_name_parts(filenames):
    """
    Find the common leading, underscore-separated name shared by FASTQ filenames.

    Strips the FASTQ extension (.fastq/.fq, optionally followed by .gz/.gzip,
    case-insensitive) from each name, splits the names on underscores, and
    keeps tokens from the start for as long as they are identical across all
    names. Tokens are compared case-insensitively; the returned tokens use
    the casing of the first filename. Used to name the outputs of a paired
    R1/R2 input, e.g. sample1_R1.fastq.gz and sample1_R2.fastq.gz share the
    prefix "sample1".

    Args:
        filenames (list[str]): FASTQ filenames (basenames, not paths).

    Returns:
        str: The common tokens, rejoined by underscores. If the very first
        token already differs, the first filename's stem (without extension)
        is returned instead. An empty input list returns "unknown".

    Examples:
        - ["sample1_R1.fastq.gz", "sample1_R2.fastq.gz"] -> "sample1"
        - ["file.fq.gz"] -> "file"
        - ["a_b_c.fastq", "x_y_z.fastq"] -> "a_b_c"  (first file's stem as fallback)
    """
    if not filenames:
        return "unknown"
    stems = []
    for file in filenames:
        stems.append(re.sub(r'\.(fastq|fq)(\.gzip|\.gz)?$', '', file, flags=re.IGNORECASE))
    token_lists = []
    for stem in stems:
        token_lists.append(stem.split('_'))
    min_len = None
    for tokens in token_lists:
        if min_len is None or len(tokens) < min_len:
            min_len = len(tokens)
    common_tokens = []
    for i in range(min_len):
        values_at_position = set()
        for tokens in token_lists:
            values_at_position.add(tokens[i].lower())
        if len(values_at_position) != 1:
            break
        common_tokens.append(token_lists[0][i])
    if common_tokens:
        output = '_'.join(common_tokens)
    else:
        output = stems[0]    
    return output

def basename_file(filepath):
    """
    Extract the sample name from a FASTQ filepath by removing directory and extension.
    
    Strips the directory path and FASTQ file extension (.fastq, .fq, optionally
    .gz/.gzip) to isolate the sample identifier. Extension matching is case-insensitive.
    Names without a FASTQ extension are returned unchanged (apart from the
    directory part).

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

def lazy_fastq_blobs(filepath, chunk_size, buffer_size = 2*1024*1024):
    """
    Lazily yield raw FASTQ text in blobs of exactly ``chunk_size`` records,
    for sending to the workers.

    The main process only frames the data: newline positions are found with
    numpy, one pass per block read, and the stream is cut after every
    4 * ``chunk_size`` lines. Each blob is sent as a single bytes object, so
    it pickles as one copy instead of one object per read. Nothing is
    split, stripped or validated here; the workers do that (parse_blob(),
    validate_fastq()).

    Newline positions already found are carried over between reads, so a
    blob larger than ``buffer_size`` is assembled without rescanning.

    At the end of the file, trailing blank lines are ignored. If the file
    ends with an incomplete record, a warning is logged and that record is
    skipped, for plain-text and gzip input alike.

    Whether the file is gzip-compressed is looked up in GZIP_DETECTION,
    which must already contain ``filepath`` (see input_handler()).

    Args:
        filepath (str): Path to the FASTQ file (plain or gzip-compressed).
        chunk_size (int): Number of FASTQ records per blob.
        buffer_size (int): Number of (decompressed) bytes requested per
            read. Defaults to 2 MiB.

    Yields:
        tuple[bytes, int]: (blob, n_records). ``blob`` is raw FASTQ text
            ending on a record boundary, including its final newline.
            ``n_records`` equals ``chunk_size`` for every blob except
            possibly the last.
    """
    lines_per_blob = 4 * chunk_size
    if GZIP_DETECTION[filepath]:
        fastq_file = igzip_threaded.open(filepath, 'rb', threads=1)
    else:
        fastq_file = open(filepath, 'rb')
    with fastq_file:
        pending = b""
        pending_newlines = np.empty(0, dtype=np.intp)
        while True:
            block = fastq_file.read(buffer_size)
            if not block:
                break
            new_newlines = np.flatnonzero(np.frombuffer(block, dtype=np.uint8) == 10)
            if pending:
                new_newlines += len(pending)
                data = pending + block
                newlines = np.concatenate((pending_newlines, new_newlines))
            else:
                data = block
                newlines = new_newlines
            n_blobs = len(newlines) // lines_per_blob
            start = 0
            for k in range(1, n_blobs + 1):
                end = int(newlines[k * lines_per_blob - 1]) + 1
                yield data[start:end], chunk_size
                start = end
            pending = data[start:]
            pending_newlines = newlines[n_blobs * lines_per_blob:] - start
        tail = pending.rstrip(b"\r\n")
        if not tail:
            return
        tail_newlines = np.flatnonzero(np.frombuffer(tail, dtype=np.uint8) == 10)
        n_lines = len(tail_newlines) + 1
        usable = n_lines - (n_lines % 4)
        if usable != n_lines:
            logger.warning(
                "%s: last record is incomplete (%d trailing line(s)); it was skipped. "
                "The file may be truncated.",
                os.path.basename(filepath), n_lines - usable)
        if usable:
            blob = tail + b"\n" if usable == n_lines else tail[:int(tail_newlines[usable - 1]) + 1]
            yield blob, usable // 4

def parse_blob(blob):
    """
    Split a raw FASTQ blob from lazy_fastq_blobs() into records. Runs in
    the worker processes.

    Carriage returns are removed (CRLF input); other whitespace is kept,
    as in the gzip branch of lazy_fastq(). The blob always ends on a record
    boundary; record contents are checked later by validate_fastq().

    Args:
        blob (bytes): Raw FASTQ text from lazy_fastq_blobs().

    Returns:
        Iterator[tuple[bytes, bytes, bytes, bytes]]: (header, sequence,
            plus, quality) tuples, line endings removed. The iterator can
            only be consumed once.
    """
    if b"\r" in blob:
        blob = blob.replace(b"\r", b"")
    line = iter(blob.split(b"\n"))
    return zip(line, line, line, line)

def find_paired_files(filepaths):
    """
    Sort FASTQ files into paired-end pairs, interleaved files and unpaired files.

    Only the first two records (first 8 lines) of each file are read. The
    base read ID and read number (1 or 2) of the first header are extracted
    with read_info_from_header(), independent of the filename. A file whose
    first two records share the same base ID is marked as interleaved.

    Files are then grouped by the base ID of their first read:
      - Exactly two files with read numbers 1 and 2 form a pair, returned in
        (R1, R2) order.
      - Two files with the same read number (possible duplicates), or with
        missing/unexpected read numbers, are logged and not paired.
      - More than two files sharing a base ID are logged and not paired.
    Of the files that were not paired, those marked as interleaved are
    returned as interleaved and the rest as unpaired.

    Requires GZIP_DETECTION to contain every path (see input_handler()).

    Args:
        filepaths (list[str] | None): Paths to FASTQ files to inspect.
            None is treated as an empty list.

    Returns:
        tuple[list[tuple[str, str]], list[str], list[str]]:
            - Matched (R1, R2) file pairs.
            - Interleaved files (first two records share a base ID).
            - Remaining unpaired files.
        Files that cannot be read, or whose first header is empty or
        missing, are logged and left out of all three lists.
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
                    file = igzip.GzipFile(fileobj=raw)
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
        base_id_1, read_num_1 = read_info_from_header(header_1)
        header_2 = headers[1].strip().lstrip(b'@')
        base_id_2, read_num_2 = read_info_from_header(header_2)
    
        base_ids[filepath] = (base_id_1, read_num_1)
        interleaved_flag[filepath] = (
            base_id_1 is not None
            and base_id_1 == base_id_2
        )

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
    Extract the base read ID and read number from a FASTQ header.

    Handles two common conventions:
      - Modern Illumina (Casava 1.8+): "<id> 1:N:0:..." or "<id> 2:N:0:...".
        The base ID is everything before the first space; the read number
        is taken from the start of the second field.
      - Legacy Illumina, and MGI: "<id>/1" or "<id>/2". The base ID is the
        header without the trailing "/1" or "/2".

    If neither convention matches, the entire header is used as the base ID
    and the read number is set to None.

    Args:
        header (bytes): The FASTQ header line. Surrounding whitespace is
            stripped. A leading '@' is not removed, and is then part of the
            base ID.

    Returns:
        tuple[bytes, int | None]:
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
    Detect the Phred quality encoding offset (33 or 64) of a FASTQ file.

    If ``phred_offset`` is given, it is returned without reading the file.
    Otherwise up to ``reads_for_phred_offset`` records are sampled and the
    lowest and highest quality characters are recorded:
      - any character below ASCII 64 ('@') -> Phred+33 (Sanger/modern Illumina);
      - otherwise, if no character is above ASCII 104 ('h') -> Phred+64
        (older Illumina);
      - otherwise the encoding is ambiguous and an error is raised.

    Args:
        filepath (str): Path to the FASTQ file to inspect.
        reads_for_phred_offset (int): Maximum number of reads to sample.
        phred_offset (int | None): User-specified Phred offset, if known.

    Returns:
        int: 33 or 64, the detected or provided Phred offset.

    Raises:
        ValueError: If the FASTQ file cannot be read, or if the observed
            ASCII range fits neither encoding.
    """
    if phred_offset is not None:
        return phred_offset
    try:
        blob = next(lazy_fastq_blobs(filepath, reads_for_phred_offset), (b"",))[0]
    except (FileNotFoundError, IOError) as error:
        raise ValueError(f"Cannot read FASTQ file '{filepath}': {error}") from error
    qualities = b"".join(record[3] for record in parse_blob(blob))
    if not qualities:
        raise ValueError(f"No FASTQ records found in '{filepath}'; cannot detect Phred offset.")
    q_bytes = np.frombuffer(qualities, dtype=np.uint8)
    min_ascii = int(q_bytes.min())
    max_ascii = int(q_bytes.max())
    if min_ascii < 64:
        return 33
    if max_ascii <= 104:
        return 64
    min_char = chr(min_ascii) if 32 <= min_ascii <= 126 else '?'
    max_char = chr(max_ascii) if 32 <= max_ascii <= 126 else '?'
    raise ValueError(
        f"Ambiguous Phred encoding detected in {filepath} (ASCII range {min_ascii}-{max_ascii} ['{min_char}' - '{max_char}']). "
        f"Please specify the Phred offset (33/64) manually using the --phred-offset option."
    )

def validate_fastq(header, sequence, plus, quality, n_filter, min_length_input, max_length_input):
    """
    Check that a single FASTQ record is well-formed.

    A record passes if:
      - the header is longer than one character and starts with '@';
      - the plus line starts with '+';
      - the sequence and quality lines have the same length;
      - the sequence length lies within [min_length_input, max_length_input]
        (each bound is skipped if None);
      - the sequence contains only uppercase A, C, G and T, plus N unless
        ``n_filter`` is set (lowercase bases are rejected);
      - every quality character is printable ASCII (33-126).

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
        bool: True if the record passes all checks, False otherwise.
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
    Parse adapter sequences from a FASTA file.

    Each entry is a header line starting with ">" followed by one or more
    sequence lines. Sequence lines are uppercased and concatenated, so a
    sequence may be split over several lines. Trailing whitespace is
    stripped from every line.

    Note: parsing stops at the first empty (or whitespace-only) line,
    including one at the very top of the file, so the file should not
    contain blank lines.

    Args:
        fasta_file (str): Path to the FASTA file containing adapter sequences.

    Returns:
        list[tuple[str, str]]: A list of (name, sequence) tuples, where:
            - name (str): The header line without the leading ">" and
              surrounding whitespace.
            - sequence (str): The full, uppercased nucleotide sequence.

    Raises:
        ValueError: If an entry does not start with a ">" header line, or a
            header is not followed by any sequence lines.
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
    
def seq_to_array(sequence_list):
    """
    Convert a list of nucleotide sequences into a 2D numpy array of ASCII
    codes, and work out the chunk's padding layout for all downstream steps.

    If all reads have the same length, the sequences are joined directly.
    Otherwise each read is right-padded with null bytes (0x00) to the length
    of the longest read. The raw bytes are then reinterpreted as a
    (n_reads, max_len) matrix for vectorized processing.

    Args:
        sequence_list (list[bytes]): Sequence byte-strings (non-empty list).

    Returns:
        tuple: A tuple containing:
            - numpy.ndarray: Signed 8-bit integer array of shape
              (n_reads, max_len) of ASCII character codes.
            - bool: chunk_padding_bool, True if reads differ in length.
            - numpy.ndarray | int: row_tilde_count, (n_reads,) number of
              padding bytes per read, or 0 when not padded.
            - numpy.ndarray | None: padding_mask_bool, (n_reads, max_len)
              True at padding positions, or None when not padded.
            - numpy.ndarray: (n_reads,) real read lengths.
            - int: max_len, length of the longest read.
            - int: n_reads, number of reads.
    """
    sequence_arr_lengths = np.array([len(s) for s in sequence_list])
    n_reads = len(sequence_list)
    max_len = sequence_arr_lengths.max()
    min_len = sequence_arr_lengths.min()
    if max_len == min_len:
        joined = b''.join(sequence_list)
        array = np.frombuffer(joined, dtype=np.int8).reshape(n_reads, max_len)
        chunk_padding_bool = False
        row_tilde_count = 0
        padding_mask_bool = None
        return array, chunk_padding_bool, row_tilde_count, padding_mask_bool, sequence_arr_lengths, max_len, n_reads
    else:
        padded = [s.ljust(max_len, b'\x00') for s in sequence_list]
        joined = b''.join(padded)
        array = np.frombuffer(joined, dtype=np.int8).reshape(n_reads, max_len)
        padding_mask_bool = array == ord('\x00')
        row_tilde_count = np.sum(padding_mask_bool, axis=1)
        chunk_padding_bool = bool(np.any(row_tilde_count > 0))
        return array, chunk_padding_bool, row_tilde_count, padding_mask_bool, sequence_arr_lengths, max_len, n_reads

def qual_to_array(quality_list, phred_offset, chunk_padding_bool, padding_mask_bool, max_len, n_reads):
    """
    Convert a list of Phred quality strings into a 2D numeric numpy array,
    using the padding layout seq_to_array() determined for the same reads.

    Concatenates all quality strings (right-padded with null bytes when
    chunk_padding_bool is True), reinterprets the raw bytes as integers,
    subtracts the Phred offset (padding positions stay 0), and reshapes
    into a (n_reads, max_len) matrix.

    Args:
        quality_list (list[bytes]): Quality byte-strings, each the same length
            as its sequence (guaranteed by validate_fastq).
        phred_offset (int): The Phred encoding offset (33 or 64) to subtract.
        chunk_padding_bool (bool): From seq_to_array.
        padding_mask_bool (numpy.ndarray | None): From seq_to_array.
        max_len (int): From seq_to_array.
        n_reads (int): From seq_to_array.

    Returns:
        numpy.ndarray: Signed 8-bit integer array of shape
            (n_reads, max_len) with true quality scores.
    """
    if not chunk_padding_bool:
        joined = b''.join(quality_list)
        array = np.frombuffer(joined, dtype=np.int8).reshape(n_reads, max_len)
        return array - phred_offset
    padded = [q.ljust(max_len, b'\x00') for q in quality_list]
    joined = b''.join(padded)
    array = np.frombuffer(joined, dtype=np.int8).reshape(n_reads, max_len)
    return np.where(padding_mask_bool, array, array - phred_offset)

def header_mgi_to_illumina(mgi_header, barcode5, barcode7, instrument, run):
    """
    Convert an MGI/BGI FASTQ header into an Illumina (Casava 1.8+) header.

    MGI headers have the form "<flowcell>L<lane>C<column>R<row><tile>/<read>",
    where <row> is three digits. The Illumina fields missing from an MGI
    header (instrument, run and index sequences) are supplied by the caller.
    The result has the form:

        @<instrument>:<run>:<flowcell>:<lane>:<tile>:<column>:<row> <read>:N:0:<barcode7>+<barcode5>

    Headers that already look like Illumina headers are returned unchanged
    (with a leading '@'). Anything after "/<read>" in an MGI header is
    discarded.

    Args:
        mgi_header (bytes): MGI header, with or without leading '@'.
        barcode5 (str): The i5 index sequence or barcode string.
        barcode7 (str): The i7 index sequence or barcode string.
        instrument (str): The sequencing instrument identifier.
        run (str): The run identifier or run number.

    Returns:
        bytes: Illumina-style header (including leading '@').

    Raises:
        ValueError: If the header matches neither the MGI nor the Illumina format.
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
        barcode7.encode('ascii'),
        barcode5.encode('ascii')
    )
    return illumina_header

def build_pipeline(parameters, read_direction=None):
    """
    Build the list of enabled trimming and filtering modules.

    Each returned callable takes (sequence_arr, quality_arr,
    chunk_padding_bool, row_tilde_count, padding_mask_bool) for one chunk
    and returns a (left, right) pair of per-read boundary arrays, describing
    the [left, right) window of bases that module would keep. Callers merge
    all windows by taking the largest left and smallest right boundary per
    read, so the order of the modules does not affect the result.

    Modules included, when enabled in ``parameters``:
      - k-mer complexity filter ("kmer_filter_flag")
      - N end trimming ("n_trimming_flag")
      - homopolymer end trimming ("poly_filter_flag")
      - adapter trimming ("adapter_filter_flag")
      - quality-dependent end trimming ("endqual_filter_flag")
      - input average-quality filter ("minimum_average_qual_pre" > 0)
      - sliding-window quality trimming ("slider_filter_flag")
      - set-length end trimming ("cut_flag")
    Overlap trimming for paired reads is not part of this list; it is
    computed per pair and added separately in finish_reads().

    Args:
        parameters (dict): Run parameters (module flags and their settings).

    Returns:
        list[callable]: The enabled modules; empty if none are enabled.
    """
    pipeline = []
    #sequence_based
    if parameters.get("kmer_filter_flag"):
        pipeline.append(lambda seq, qual, chunk_padding_bool, row_tilde_count, padding_mask_bool: kmer_complexity_scan(seq, chunk_padding_bool, padding_mask_bool, kmer=parameters["kmer_size"], low_complex_cutoff=parameters["kmer_cutoff"], allow_n=parameters["allow_n_kmer"]))
    if parameters.get("n_trimming_flag"):
        pipeline.append(lambda seq, qual, chunk_padding_bool, row_tilde_count, padding_mask_bool: n_end_trimming(seq, padding_mask_bool, chunk_padding_bool, row_tilde_count))
    if parameters.get("poly_filter_flag"):
        pipeline.append(lambda seq, qual, chunk_padding_bool, row_tilde_count, padding_mask_bool: homopolymer_nucleotide_trimming(seq, padding_mask_bool, row_tilde_count, chunk_padding_bool, poly_length_both=parameters["poly_length_both"], poly_length_start=parameters["poly_length_start"], poly_length_end=parameters["poly_length_end"], poly_bases_both=parameters["poly_bases_both"], poly_bases_start=parameters["poly_bases_start"], poly_bases_end=parameters["poly_bases_end"]))
    if parameters.get("adapter_filter_flag"): 
        pipeline.append(lambda seq, qual, chunk_padding_bool, row_tilde_count, padding_mask_bool: adapter_trimming(seq, chunk_padding_bool, row_tilde_count, adapter_sequences=parameters["adapter_sequences"], mismatches=parameters["adapter_mismatch"], adapter_seed=parameters["adapter_seed"], read_direction=read_direction, adapter_group=None if parameters.get("adapter_fasta_excl") else parameters.get("adapter_group")))
        
    #quality_based
    if parameters.get("endqual_filter_flag"):
        pipeline.append(lambda seq, qual, chunk_padding_bool, row_tilde_count, padding_mask_bool: trim_ends_quality(qual, chunk_padding_bool, row_tilde_count, padding_mask_bool, min_quality_both=parameters["min_quality_both"], endqual_min_start=parameters["endqual_min_start"], endqual_min_end=parameters["endqual_min_end"]))
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
    Open a new output FASTQ file for writing.

    The output name is basename_file(filepath) followed by
    "_filtered.fastq", or "_filtered.fastq.gz" when gzip_output is True.
    If the resulting name contains "rejected", the "_filtered" part is
    dropped, so rejected-read files are named "<name>_rejected.fastq[.gz]".
    ``filepath`` does not need to be an existing file: input_handler() also
    passes plain name stems such as "<prefix>_R1_paired".

    The file is opened in exclusive binary write mode ("xb"). This function
    only chooses the name; the data itself is compressed by the workers.

    Args:
        filepath (str): Input path or name stem to derive the output name from.
        output_dir (str): Directory in which to create the output file.
        gzip_output (bool): If True, use a ".fastq.gz" extension.

    Returns:
        io.BufferedWriter: The opened binary file handle.

    Raises:
        FileExistsError: If the output file already exists.
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
    Whole-read average-quality filter (--min-average-qual-pre), in the
    pipeline's (left, right) boundary format.

    The mean quality is computed over each read's full, untrimmed length.
    Reads with a mean >= min_avg_qual keep their full window
    (left 0, right = read length). Reads below the threshold get right 0,
    i.e. a zero-length window, so they are discarded downstream.

    Args:
        quality_arr (numpy.ndarray): (n_reads, read_length) array of
            per-base quality scores.
        chunk_padding_bool (bool): Whether quality_arr contains padded
            (unequal-length) reads.
        row_tilde_count (numpy.ndarray | int): (n_reads,) count of padding
            bytes per read, used to recover each read's real length when
            chunk_padding_bool is True.
        min_avg_qual (float): Minimum average quality required to keep a read.

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: (left_cutoffs, right_cutoffs),
            each of shape (n_reads,). Left cutoffs are all 0 (int8); right
            cutoffs are int32.
    """
    n_reads, length = quality_arr.shape
    if not chunk_padding_bool:
        real_lengths = np.full(n_reads, length, dtype=np.int32)
    else:
        real_lengths = (length - row_tilde_count).astype(np.int32)
    avg_quals = average_quality_batch(quality_arr, lefts=0, rights=real_lengths)
    passed = avg_quals >= min_avg_qual
    right_cutoffs = np.where(passed, real_lengths, 0).astype(np.int32)
    return np.zeros(n_reads, dtype=np.int8), right_cutoffs

def trim_ends_quality(quality_arr, chunk_padding_bool, row_tilde_count, padding_mask_bool, min_quality_both, endqual_min_start, endqual_min_end):
    """
    Quality-dependent end trimming: per-read boundaries that trim each end
    up to the first base that meets a quality threshold.

    The left boundary is the first position (scanning from the 5' end) with
    quality >= endqual_min_start. The right boundary is one past the last
    position (scanning from the 3' end) with quality >= endqual_min_end.
    Only the ends are trimmed; low-quality bases inside the read are kept.
    Padding positions are ignored. Reads without any base meeting the
    threshold get left = row width and right = 0, i.e. an empty window.
    Vectorized across all reads in the chunk.

    Args:
        quality_arr (numpy.ndarray): (n_reads, read_length) array of
            per-base quality scores.
        chunk_padding_bool (bool): Whether quality_arr contains padded
            (unequal-length) reads.
        row_tilde_count (numpy.ndarray | int): (n_reads,) count of padding
            bytes per read, used when chunk_padding_bool is True.
        padding_mask_bool (numpy.ndarray | None): (n_reads, read_length)
            boolean mask marking right-padding positions, or None if the
            chunk isn't padded. Padded positions are excluded from the scan.
        min_quality_both (int | None): Threshold used for any end whose
            specific threshold is None. If this is also None, 0 is used
            (no trimming at that end).
        endqual_min_start (int | None): Minimum quality to keep bases from
            the start of the read, or None to use min_quality_both.
        endqual_min_end (int | None): Minimum quality to keep bases from
            the end of the read, or None to use min_quality_both.

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: (start_cutoffs, end_cutoffs),
            each of shape (n_reads,) and dtype int32, giving the left and
            right trim boundaries per read.
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

        real_lengths = length - row_tilde_count
        run_length = first_good_pos - row_tilde_count
        end_cutoffs = real_lengths - run_length

        zero_end_rows = first_good_pos == 0
        if zero_end_rows.any():
            end_good_pos = qual_mask[:, 0] | ~zero_end_rows
            end_cutoffs = np.where(end_good_pos, end_cutoffs, 0)

    return start_cutoffs.astype(np.int32), end_cutoffs.astype(np.int32)

def homopolymer_nucleotide_trimming(sequence_arr, padding_mask_bool, row_tilde_count, chunk_padding_bool, poly_length_both, poly_length_start, poly_length_end, poly_bases_both, poly_bases_start, poly_bases_end):
    """
    Homopolymer trimming: per-read boundaries that remove homopolymer runs
    from the start and/or end of each read.

    For each base to check at an end, the run of that base touching the end
    is measured; if it is at least the minimum length for that end, the whole
    run is trimmed. Each listed base is checked on its own (a mixed run such
    as "GGAGG" is not treated as one run), and the largest trim per end wins.

    Base selection: poly_bases_both applies to both ends, but
    poly_bases_start and poly_bases_end each replace it for their own end
    when they are not None. Length selection: poly_length_start and
    poly_length_end are used for their own end, and fall back to
    poly_length_both when 0.

    If no bases are given at all, or all three lengths are 0, no trimming is
    done. A read consisting entirely of a checked base ends up with an empty
    window.

    Args:
        sequence_arr (numpy.ndarray): (n_reads, read_length) array of
            per-base ASCII sequence codes.
        padding_mask_bool (numpy.ndarray | None): (n_reads, read_length)
            boolean mask marking right-padding positions, or None if the
            chunk isn't padded.
        row_tilde_count (numpy.ndarray | int): (n_reads,) count of padding
            bytes per read, used when chunk_padding_bool is True.
        chunk_padding_bool (bool): Whether sequence_arr contains padded
            (unequal-length) reads.
        poly_length_both (int): Minimum run length used for an end whose
            specific length is 0.
        poly_length_start (int): Minimum run length to trigger trimming at
            the start of the read (0 = use poly_length_both).
        poly_length_end (int): Minimum run length to trigger trimming at
            the end of the read (0 = use poly_length_both).
        poly_bases_both (str | None): Comma-separated bases to check at
            both ends.
        poly_bases_start (str | None): Comma-separated bases to check at
            the start; replaces poly_bases_both for the start if not None.
        poly_bases_end (str | None): Comma-separated bases to check at
            the end; replaces poly_bases_both for the end if not None.

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: (left_cutoffs, right_cutoffs),
            each of shape (n_reads,), giving the left and right trim
            boundaries per read.

    Raises:
        ValueError: If any comma-separated base entry is not a single character.
    """
    n_reads, length = sequence_arr.shape
    if not poly_bases_start and not poly_bases_end and not poly_bases_both:
        return np.zeros(n_reads, dtype=np.int8), np.full(n_reads, length, dtype=np.int32)
    if poly_length_both == poly_length_start == poly_length_end == 0:
        return np.zeros(n_reads, dtype=np.int8), np.full(n_reads, length, dtype=np.int32)

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

    right_cutoffs = np.full(n_reads, length, dtype=np.int32)
    left_cutoffs = np.zeros(n_reads, dtype=np.int8)

    for base in start_bases:
        base_code = ord(base)
        non_base_mask = sequence_arr != base_code
        first_non_pos = non_base_mask.argmax(axis=1)
        first_non_pos = np.where(non_base_mask.any(axis=1), first_non_pos, length)
        trim_amount = np.where(first_non_pos >= poly_length_start, first_non_pos, 0)
        left_cutoffs = np.maximum(left_cutoffs, trim_amount)

    if end_bases:
        rev_seq = np.ascontiguousarray(sequence_arr[:, ::-1])
        if not chunk_padding_bool:
            for base in end_bases:
                base_code = ord(base)
                non_base_mask = rev_seq != base_code
                first_non_pos = non_base_mask.argmax(axis=1)
                first_non_pos = np.where(non_base_mask.any(axis=1), first_non_pos, length)
                trim_amount = np.where(first_non_pos >= poly_length_end, first_non_pos, 0)
                base_right_cutoffs = length - trim_amount
                right_cutoffs = np.minimum(right_cutoffs, base_right_cutoffs)
        else:
            real_length = length - row_tilde_count
            for base in end_bases:
                non_base_mask = (rev_seq != ord(base)) & (rev_seq != ord('\x00'))
                first_non_pos = non_base_mask.argmax(axis=1)
                first_non_pos = np.where(non_base_mask.any(axis=1), first_non_pos, length)
                run_length = first_non_pos - row_tilde_count
                trim_amount = np.where(run_length >= poly_length_end, run_length, 0)
                base_right_cutoffs = real_length - trim_amount
                right_cutoffs = np.minimum(right_cutoffs, base_right_cutoffs)

    return left_cutoffs, right_cutoffs

def n_end_trimming(sequence_arr, padding_mask_bool, chunk_padding_bool, row_tilde_count):
    """
    N end trimming: remove runs of N bases (any length >= 1) from both ends
    of each read, leaving internal N's untouched.

    Implemented as homopolymer_nucleotide_trimming() with base "N" at both
    ends and a minimum run length of 1.

    Args:
        sequence_arr (numpy.ndarray): (n_reads, read_length) array of
            per-base ASCII sequence codes.
        padding_mask_bool (numpy.ndarray | None): (n_reads, read_length)
            boolean mask marking right-padding positions, or None if the
            chunk isn't padded.
        chunk_padding_bool (bool): Whether sequence_arr contains padded
            (unequal-length) reads.
        row_tilde_count (numpy.ndarray | int): (n_reads,) count of padding
            bytes per read, used when chunk_padding_bool is True.

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: (left_cutoffs, right_cutoffs),
            each of shape (n_reads,), giving the trim boundaries per read.
    """
    lefts, rights = homopolymer_nucleotide_trimming(sequence_arr, padding_mask_bool = padding_mask_bool, row_tilde_count = row_tilde_count, chunk_padding_bool = chunk_padding_bool, poly_length_both = 1, poly_length_start = 0, poly_length_end = 0, poly_bases_both = "N", poly_bases_start = None, poly_bases_end = None)
    return lefts, rights

def cut_set_ends(sequence_arr, chunk_padding_bool, row_tilde_count, cut_both, cut_start, cut_end):
    """
    Set-length end trimming: remove a fixed number of bases from each end of
    every read, regardless of sequence or quality (e.g. to hard-trim a known
    primer or a biased first/last few bases).

    `cut_start` and `cut_end` set the trim for their own end. A value of 0
    means "not given" and falls back to `cut_both`, so an explicit 0 cannot
    override a nonzero `cut_both`. If the two cuts meet or cross, the read's
    window becomes empty and the read is discarded downstream.

    For equal-length chunks, the boundaries are clamped to
    [0, read_length]. For padded chunks, the boundaries are computed from
    each read's real length and are not clamped; a negative right boundary
    has the same effect (the read is discarded).

    Args:
        sequence_arr (numpy.ndarray): (n_reads, read_length) array of
            per-base ASCII sequence codes, used only for its shape.
        chunk_padding_bool (bool): Whether sequence_arr contains padded
            (unequal-length) reads.
        row_tilde_count (numpy.ndarray | int): (n_reads,) count of padding
            bytes per read, used to recover each read's real length when
            chunk_padding_bool is True.
        cut_both (int): Number of bases to trim off both ends. Used for any
            end whose specific value is 0.
        cut_start (int): Number of bases to trim from the 5' end. Takes
            priority over cut_both if nonzero.
        cut_end (int): Number of bases to trim from the 3' end. Takes
            priority over cut_both if nonzero.

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: (left_cutoffs, right_cutoffs),
            each of shape (n_reads,), giving the left and right trim
            boundaries per read.
    """
    n_reads, length = sequence_arr.shape
    cut_start = cut_start if cut_start != 0 else cut_both
    cut_end = cut_end if cut_end != 0 else cut_both
    if not chunk_padding_bool:
        cut_end_pos = length - cut_end
        cut_end_pos = max(cut_end_pos, 0)
        cut_start_pos = min(cut_start, cut_end_pos)
        cut_start_pos = max(cut_start_pos, 0)
        return np.full(n_reads, cut_start_pos, dtype=np.int32), np.full(n_reads, cut_end_pos, dtype=np.int32)
    else:
        real_length = length - row_tilde_count
        end_positions = real_length - cut_end
        cut_start_pos = np.where(cut_start > end_positions, end_positions, cut_start)
        return np.full(n_reads, cut_start_pos, dtype=np.int32), end_positions

def sliding_window_quality(quality_arr, chunk_padding_bool, padding_mask_bool, row_tilde_count, slider_quality, slider_window, slider_step):
    """
    Sliding-window quality trimming: keep the best stretch of each read in
    which every window passes the quality threshold.

    Windows of `slider_window` bases start every `slider_step` positions
    (the last possible window is always included). A window fails if its
    mean quality is below `slider_quality`. Every base covered by at least
    one failing window is marked bad, and the longest run of good bases is
    kept, with ties broken by the higher mean quality. Because the kept
    stretch can lie in the middle of the read, both ends can be trimmed at
    once, and after a mid-read quality drop the longer side survives.

    Reads shorter than one window are not trimmed. In padded chunks,
    windows overlapping padding count as failing.

    Args:
        quality_arr (numpy.ndarray): (n_reads, read_length) array of
            per-base quality scores.
        chunk_padding_bool (bool): Whether quality_arr contains padded
            (unequal-length) reads.
        padding_mask_bool (numpy.ndarray | None): (n_reads, read_length)
            boolean mask marking right-padding positions, or None if the
            chunk isn't padded.
        row_tilde_count (numpy.ndarray | int): (n_reads,) count of padding
            bytes per read, used to recover each read's real length when
            chunk_padding_bool is True.
        slider_quality (int): Minimum mean quality for a window to pass.
        slider_window (int): Number of bases per sliding window.
        slider_step (int): Step size between successive window start positions.

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: (left_cutoffs, right_cutoffs),
            each of shape (n_reads,), giving the kept [left, right) region
            per read. Reads without failing windows keep their full length;
            reads in which every base is bad get an empty region (0, 0).
    """
    n_reads, length = quality_arr.shape

    if not chunk_padding_bool:
        if length < slider_window:
            return np.zeros(n_reads, dtype=np.int8), np.full(n_reads, length, dtype=np.int32)

        last_possible_start = length - slider_window
        window_starts = np.arange(0, last_possible_start + 1, slider_step)
        if window_starts[-1] != last_possible_start:
            window_starts = np.append(window_starts, last_possible_start)

        cumsum = np.cumsum(quality_arr, axis=1, dtype=np.int32)
        cumsum = np.concatenate([np.zeros((n_reads, 1), dtype=np.int32), cumsum], axis=1)
        window_sums = cumsum[:, window_starts + slider_window] - cumsum[:, window_starts]
        failed_mask = window_sums < (slider_quality * slider_window)

        bad_positions = np.zeros((n_reads, length), dtype=bool)

        if slider_step == 1:
            for offset in range(slider_window):
                bad_positions[:, offset:offset + len(window_starts)] |= failed_mask
        else:
            for j, start in enumerate(window_starts):
                bad_positions[:, start:start + slider_window] |= failed_mask[:, j:j + 1]

        real_lengths = np.full(n_reads, length, dtype=np.int32)
    else:
        real_lengths = length - row_tilde_count
        too_short = real_lengths < slider_window

        if length < slider_window:
            left_cutoffs = np.zeros(n_reads, dtype=np.int8)
            right_cutoffs = real_lengths.astype(np.int32)
            return left_cutoffs, right_cutoffs

        cumsum = np.cumsum(quality_arr, axis=1, dtype=np.int32)
        cumsum = np.concatenate([np.zeros((n_reads, 1), dtype=np.int32), cumsum], axis=1)
        pad_cumsum = np.cumsum(padding_mask_bool.astype(np.int32), axis=1)
        pad_cumsum = np.concatenate([np.zeros((n_reads, 1), dtype=np.int32), pad_cumsum], axis=1)

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
    left_cutoffs = np.zeros(n_reads, dtype=np.int32)
    right_cutoffs = np.zeros(n_reads, dtype=np.int32)
    right_cutoffs[no_bad] = real_lengths[no_bad].astype(np.int32)
    needs_stretch_search = ~no_bad & ~all_bad
    if not needs_stretch_search.any():
        if chunk_padding_bool:
            left_cutoffs[too_short] = 0
            right_cutoffs[too_short] = real_lengths[too_short].astype(np.int32)
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
        right_cutoffs[too_short] = real_lengths[too_short].astype(np.int32)

    return left_cutoffs, right_cutoffs

def adapter_trimming_overlap(seq_arr_1, seq_arr_2, chunk_padding_bool_1, row_tilde_count_1, chunk_padding_bool_2, row_tilde_count_2, min_overlap, max_mismatch, max_mismatch_frac):
    """
    Overlap-based adapter trimming for read pairs, independent of the
    adapter sequence.

    When the DNA insert is shorter than the reads, both mates read through
    into adapter. For a pair with insert size I, the first I bases of R1
    then equal the reverse complement of the first I bases of R2. For every
    candidate I from `min_overlap` upwards (only sizes that both mates cover
    and that at least one mate reads past), these two I-base regions are
    compared. Any position where either base is N counts as a mismatch,
    including N vs N. An insert size is accepted if it has at most
    `max_mismatch` mismatches and a mismatch fraction of at most
    `max_mismatch_frac` percent of the I overlap positions. If several
    insert sizes are accepted, the one with the lowest mismatch fraction
    wins (ties go to the larger insert). Both mates are then cut at I.

    Args:
        seq_arr_1 (numpy.ndarray): (n_reads, length_1) int8 ASCII array of R1.
        seq_arr_2 (numpy.ndarray): (n_reads, length_2) int8 ASCII array of R2,
            row-aligned with seq_arr_1 (row i of both is one pair).
        chunk_padding_bool_1 (bool): Whether seq_arr_1 contains padded
            (unequal-length) reads.
        row_tilde_count_1 (numpy.ndarray | int): (n_reads,) count of padding
            bytes per R1 read, used when chunk_padding_bool_1 is True.
        chunk_padding_bool_2 (bool): As chunk_padding_bool_1, for R2.
        row_tilde_count_2 (numpy.ndarray | int): As row_tilde_count_1, for R2.
        min_overlap (int): Shortest insert size tested.
        max_mismatch (int): Maximum mismatches allowed in the overlap
            (N counts as a mismatch).
        max_mismatch_frac (float): Maximum percentage (0-100) of overlap
            positions that may mismatch (N counts as a mismatch).

    Returns:
        tuple[tuple[numpy.ndarray, numpy.ndarray], tuple[numpy.ndarray, numpy.ndarray]]:
            ((left_cutoffs_1, right_cutoffs_1), (left_cutoffs_2, right_cutoffs_2)),
            each of shape (n_reads,). Left cutoffs are always 0 (3'-end trimming
            only, dtype int8); right cutoffs (int32) are the accepted insert
            size, or the read's real length if no overlap was accepted.

    Raises:
        ValueError: If seq_arr_1 and seq_arr_2 have different numbers of rows.
    """
    prefix_length = 2 * max_mismatch + 6
    max_mismatch_frac = max_mismatch_frac / 100

    n_base = np.int8(ord("N"))

    n_reads_1, length_1 = seq_arr_1.shape
    n_reads_2, length_2 = seq_arr_2.shape
    if n_reads_1 != n_reads_2:
        raise ValueError(f"R1 has {n_reads_1} reads, R2 has {n_reads_2}")

    rc_seq_arr_2 = np.frombuffer(seq_arr_2[:, ::-1].tobytes().translate(COMP_TABLE), dtype=np.int8).reshape(n_reads_2, length_2)    
    seq_arr_1 = np.where(seq_arr_1 == n_base, np.int8(-1), seq_arr_1)
    rc_seq_arr_2 = np.where(rc_seq_arr_2 == n_base, np.int8(-2), rc_seq_arr_2)
    
    best_insert = np.full(n_reads_1, -1, dtype=np.int32)
    max_insert_size = min(length_1, length_2)
    any_padding = chunk_padding_bool_1 or chunk_padding_bool_2

    real_lengths_1 = (length_1 - row_tilde_count_1 if chunk_padding_bool_1
                      else np.full(n_reads_1, length_1)).astype(np.int32)
    real_lengths_2 = (length_2 - row_tilde_count_2 if chunk_padding_bool_2
                      else np.full(n_reads_2, length_2)).astype(np.int32)

    if any_padding:
        shorter_mate = np.minimum(real_lengths_1, real_lengths_2)
        longer_mate = np.maximum(real_lengths_1, real_lengths_2)
        last_insert = int(np.minimum(shorter_mate, longer_mate - 1).max())
    else:
        last_insert = min(max_insert_size, max(length_1, length_2) - 1)

    for insert_size in range(min_overlap, last_insert + 1):
        offset = length_2 - insert_size
        prefix = min(prefix_length, insert_size)
        
        if any_padding:
            candidate = (insert_size <= shorter_mate) & (insert_size < longer_mate)
            candidate_rows = np.flatnonzero(candidate)
        else:
            candidate_rows = None

        if candidate_rows is None:
            prefix_arr_1 = seq_arr_1[:, :prefix]
            prefix_arr_2 = rc_seq_arr_2[:, offset:offset + prefix]
        else:
            prefix_arr_1 = seq_arr_1[candidate_rows, :prefix]
            prefix_arr_2 = rc_seq_arr_2[candidate_rows, offset:offset + prefix]
        prefix_mismatches = (prefix_arr_1 != prefix_arr_2).sum(axis=1)

        survive_bool = prefix_mismatches <= max_mismatch
        if not survive_bool.any():
            continue
        rows_survived = np.flatnonzero(survive_bool) if candidate_rows is None else candidate_rows[survive_bool]
        mismatches = prefix_mismatches[survive_bool]

        if insert_size > prefix:
            rest_arr_1 = seq_arr_1[rows_survived, prefix:insert_size]
            rest_arr_2 = rc_seq_arr_2[rows_survived, offset + prefix:]
            mismatches = mismatches + (rest_arr_1 != rest_arr_2).sum(axis=1)

        fraction_mismatch = mismatches / insert_size
        pass_bool = (mismatches <= max_mismatch) & (fraction_mismatch <= max_mismatch_frac)
        rows_passed = rows_survived[pass_bool]
        best_insert[rows_passed] = insert_size

    found = best_insert > 0
    right_cutoffs_1 = np.where(found, best_insert, real_lengths_1).astype(np.int32)
    right_cutoffs_2 = np.where(found, best_insert, real_lengths_2).astype(np.int32)
    return ((np.zeros(n_reads_1, dtype=np.int8), right_cutoffs_1), (np.zeros(n_reads_2, dtype=np.int8), right_cutoffs_2))

def select_adapters_for_read(adapter_sequences, adapter_group, read_direction):
    """
    For paired reads (read_direction "read_1"/"read_2"), drop the other mate's
    directional adapter (a DEFAULT_ADAPTERS entry named Read_1/Read_2) of the
    selected adapter groups. Sequences that another selected entry also needs
    are kept. With no direction or no group, the list is returned unchanged.
    """
    if read_direction not in ("read_1", "read_2") or not adapter_group:
        return adapter_sequences
    other_name = "Read_2" if read_direction == "read_1" else "Read_1"
    drop, keep = set(), set()
    for group, entries in DEFAULT_ADAPTERS:
        if group not in adapter_group:
            continue
        for name, seq in entries:
            data = seq.encode("utf-8")
            if name == other_name:
                drop.add(data)
            else:
                keep.add(data)
    return [a for a in adapter_sequences if a not in drop or a in keep]

def adapter_trimming(sequence_arr, chunk_padding_bool, row_tilde_count, adapter_sequences, mismatches, adapter_seed, read_direction=None, adapter_group=None):
    """
    Adapter trimming: cut each read at the earliest occurrence of any
    adapter sequence.

    For each read, the earliest position where any of the adapters occurs is
    found, and everything from that position onwards is removed (3'-end
    trimming). With mismatches == 0, exact substring search is used. With
    mismatches > 0, a vectorized Hamming-distance search allows up to that
    many substitutions (no insertions or deletions). The whole adapter must
    lie within the read, so partial adapters at the very 3' end are not
    detected (for paired reads, overlap trimming can catch those).

    Args:
        sequence_arr (numpy.ndarray): (n_reads, read_length) array of
            per-base ASCII sequence codes (int8).
        chunk_padding_bool (bool): Whether sequence_arr contains padded
            (unequal-length) reads.
        row_tilde_count (numpy.ndarray | int): (n_reads,) count of padding
            bytes per read, used to keep matches within each read's real
            length when chunk_padding_bool is True.
        adapter_sequences (list[bytes]): Adapter sequences to search for.
        mismatches (int): Number of allowed substitutions. 0 means exact
            matching.

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: (left_cutoffs, right_cutoffs),
            each of shape (n_reads,). Left cutoffs are always 0 (int8);
            right cutoffs (int32) are the adapter start position, or the
            read's real length if no adapter was found.
    """
    n_reads, length = sequence_arr.shape
    adapter_sequences = select_adapters_for_read(adapter_sequences, adapter_group, read_direction)
    if mismatches == 0:
        all_bytes = sequence_arr.tobytes()
        if not chunk_padding_bool:
            real_lengths = np.full(n_reads, length, dtype=np.int32)
        else:
            real_lengths = length - row_tilde_count
        right_cutoffs = np.array(real_lengths, dtype=np.int32)
        for adapter_bytes in adapter_sequences:
            adapter_len = len(adapter_bytes)
            seed_len = min(adapter_seed, adapter_len)
            seed = adapter_bytes[:seed_len]
            start = 0
            while True:
                pos = all_bytes.find(seed, start)
                if pos == -1:
                    break
                i, pos_in_row = divmod(pos, length)
                avail = real_lengths[i] - pos_in_row
                if avail >= adapter_len:
                    if all_bytes.startswith(adapter_bytes, pos) and pos_in_row < right_cutoffs[i]:
                        right_cutoffs[i] = pos_in_row
                    start = (i + 1) * length
                    continue
                if avail == seed_len:
                    if pos_in_row < right_cutoffs[i]:
                        right_cutoffs[i] = pos_in_row
                    start = (i + 1) * length
                    continue
                if avail > seed_len:
                    if all_bytes[pos:pos + avail] == adapter_bytes[:avail] and pos_in_row < right_cutoffs[i]:
                        right_cutoffs[i] = pos_in_row
                    start = (i + 1) * length
                    continue
                start = pos + 1
    else:
        if not chunk_padding_bool:
            right_cutoffs = np.full(n_reads, length, dtype=np.int32)
            for adapter_bytes in adapter_sequences:
                adapter_arr = np.frombuffer(adapter_bytes, dtype=np.int8)
                adapter_len = len(adapter_bytes)
                if adapter_len > length:
                   continue
                n_windows = length - adapter_len + 1
                mismatch_matrix = np.zeros((n_reads, n_windows), dtype=np.uint8)
                for j in range(adapter_len):
                    mismatch_matrix += sequence_arr[:, j:j + n_windows] != adapter_arr[j]
                valid_mask = mismatch_matrix <= mismatches
                has_match = valid_mask.any(axis=1)

                if has_match.any():
                    first_match_col = valid_mask.argmax(axis=1)
                    right_cutoffs[has_match] = np.minimum(
                        right_cutoffs[has_match],
                        first_match_col[has_match].astype(np.int32)
                    )
        else:
            real_lengths = length - row_tilde_count
            right_cutoffs = np.full(n_reads, real_lengths, dtype=np.int32)
            for adapter_bytes in adapter_sequences:
                adapter_arr = np.frombuffer(adapter_bytes, dtype=np.int8)
                adapter_len = len(adapter_bytes)
                if adapter_len > length:
                   continue
                n_windows = length - adapter_len + 1
                mismatch_matrix = np.zeros((n_reads, n_windows), dtype=np.uint8)
                for j in range(adapter_len):
                    mismatch_matrix += sequence_arr[:, j:j + n_windows] != adapter_arr[j]
                window_ends = np.arange(adapter_len, length + 1)
                within_bounds = window_ends <= real_lengths[:, None]
                valid_mask = (mismatch_matrix <= mismatches) & within_bounds
                has_match = valid_mask.any(axis=1)

                if has_match.any():
                    first_match_col = valid_mask.argmax(axis=1)
                    right_cutoffs[has_match] = np.minimum(
                        right_cutoffs[has_match],
                        first_match_col[has_match].astype(np.int32)
                    )
    return np.zeros(n_reads, dtype=np.int8), right_cutoffs
    
def average_quality_batch(quality_arr, lefts, rights):
    """
    Compute the mean quality score within a per-read [left, right) window,
    vectorized across all reads at once.

    Args:
        quality_arr (numpy.ndarray): (n_reads, read_length) array of
            per-base quality scores.
        lefts (numpy.ndarray | int): Per-read left boundary (inclusive),
            shape (n_reads,), or a single value for all reads.
        rights (numpy.ndarray | int): Per-read right boundary (exclusive),
            shape (n_reads,), or a single value for all reads.

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
    Low-complexity filter: flag reads with too few distinct k-mers for removal.

    For each k-mer length, the number of distinct k-mers in a read is divided
    by the maximum possible number, min(number of k-mer windows,
    alphabet_size ** k). A read fails if this ratio is below
    `low_complex_cutoff` percent for any of the k-mer lengths. Failing reads
    get an empty window; passing reads are not trimmed.

    Bases are encoded with 3 bits each, which limits k to 21 (63 bits). The
    alphabet size is 5 (A, C, G, T, N) when `allow_n` is True, and 4
    otherwise. In padded chunks, k-mer windows overlapping padding are
    excluded from both the distinct count and the maximum.

    Args:
        sequence_arr (numpy.ndarray): A 2D array of ASCII sequence codes
            of shape (n_reads, length).
        chunk_padding_bool (bool): Whether sequence_arr contains padded
            (unequal-length) reads.
        padding_mask_bool (numpy.ndarray | None): (n_reads, length) boolean
            mask marking right-padding positions, or None if the chunk
            isn't padded.
        kmer (int | str | iterable): The k-mer length(s) to evaluate. A
            string may contain comma-separated values, e.g. "4,6".
        low_complex_cutoff (float): Minimum percentage of distinct k-mers,
            relative to the maximum possible, required to pass.
        allow_n (bool): If True, N is counted as a fifth letter of the
            alphabet.

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: (left_cutoffs, right_cutoffs),
            each of shape (n_reads,). Left cutoffs are all 0 (int8); right
            cutoffs (int32) are the chunk's row width for passing reads and
            0 for low-complexity reads.

    Raises:
        ValueError: If any k-mer length is greater than the chunk's row
            width, or greater than 21.
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
        
        n_bits = bits_per_base * k
        kmer_dtype = np.int16 if n_bits <= 15 else (np.int32 if n_bits <= 31 else np.int64)
        kmer_ints = np.zeros((n_reads, max_kmers), dtype=kmer_dtype)
        
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
            unique_counts = unique_counts - any_pad_in_row.astype(np.int32)

            denom = np.minimum(valid_kmer_counts, alphabet_size**k)

            ratio = np.where(denom > 0, unique_counts / denom, 0.0)
            global_passed &= (ratio >= (low_complex_cutoff / 100))

    second_array = np.where(global_passed, length, 0).astype(np.int32)
    return np.zeros(n_reads, dtype=np.int8), second_array

##### Unpaired reads workflow functions #####
def process_unpaired_chunk(chunk, phred_offset, minimum_average_qual_post, gzip_output, gzip_level, min_length_output, max_length_output, min_length_output_perc, max_length_output_perc, write_rejected, parameters):
    """
    Validate, trim, filter and format one chunk of unpaired FASTQ reads.

    Steps:
      1. Each record is checked with validate_fastq(); invalid records are
         rejected.
      2. If MGI conversion is enabled, headers are converted to Illumina
         format and plus lines are reset to "+".
      3. All enabled modules from build_pipeline() are run. Their windows
         are merged by taking the largest left and smallest right boundary
         per read (the most stringent result per end).
      4. A read is kept if its merged window is non-empty, lies within the
         output length limits, and (if set) has a mean quality of at least
         minimum_average_qual_post.
      5. Kept reads are cut to their window and, if requested, their
         quality strings are re-encoded to parameters["phred_out"].

    Rejected records are written untrimmed, with their original quality
    encoding.

    Args:
        chunk (Iterable[tuple[bytes, bytes, bytes, bytes]]): (header,
            sequence, plus, quality) records, as produced by parse_blob().
        phred_offset (int): Phred encoding offset (33 or 64).
        minimum_average_qual_post (float): Minimum acceptable mean quality
            after trimming. 0 disables this filter.
        gzip_output (bool): If True, compresses output records with isal gzip.
        gzip_level (int): isal gzip compression level (0-3).
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
        write_rejected (bool): If True, also collect rejected records
            (invalid on input or filtered out after trimming).
        parameters (dict): Run parameters, used to build the trimming
            pipeline and for N-filtering, input length limits,
            MGI-to-Illumina header conversion and Phred re-encoding.

    Returns:
        tuple[bytes, int, int, bytes | list]: A tuple containing:
            - Kept records as FASTQ text (gzip-compressed if gzip_output).
            - Count of kept reads.
            - Count of rejected reads (invalid input records plus records
              filtered out after trimming).
            - Rejected records as FASTQ text (gzip-compressed if
              gzip_output) if write_rejected is True, otherwise an empty list.

    Raises:
        RuntimeError: If gzip compression fails.
    """
    valid_headers = []
    valid_sequences = []
    valid_pluses = []
    valid_qualities = []
    rejected_reads = []
    rejected = 0
    for header, sequence, plus, quality in chunk:
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
        rejected_data = b"".join(rejected_reads) if write_rejected else []
        if gzip_output and rejected_data:
            rejected_data = igzip.compress(rejected_data, compresslevel=gzip_level)
        return b"", 0, rejected, rejected_data
    if parameters["mgi_convert_flag"]:
        valid_pluses = [b"+"] * len(valid_headers)
        valid_headers = [header_mgi_to_illumina(header, parameters["mgi_bc5"], parameters["mgi_bc7"], parameters["mgi_instrument"], parameters["mgi_run"]) for header in valid_headers]
    sequence_arr, chunk_padding_bool, row_tilde_count, padding_mask_bool, raw_lengths, max_len, n_reads = seq_to_array(sequence_list = valid_sequences)
    quality_arr = qual_to_array(quality_list = valid_qualities, phred_offset = phred_offset, chunk_padding_bool = chunk_padding_bool, padding_mask_bool = padding_mask_bool, max_len = max_len, n_reads = n_reads)
    n_reads, length = sequence_arr.shape
    left_list = [np.zeros(n_reads, dtype=np.int8)]
    right_list = [np.full(n_reads, length, dtype=np.int32) - row_tilde_count]
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
            results = igzip.compress(results, compresslevel=gzip_level)
            if write_rejected:
                rejected_reads = igzip.compress(rejected_reads, compresslevel=gzip_level)
        except Exception as e:
            raise RuntimeError(f"Compression failed inside worker: {str(e)}") from None
    return results, len_results, rejected, rejected_reads

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
        if GZIP_DETECTION[filepath]:
            logger.info("%s: gzip format detected.", os.path.basename(filepath))
            if ESTIMATED_ZIP_RATIO.get(filepath) is not None:
                logger.info("%s: estimated gzip compression ratio: %s.", os.path.basename(filepath), ESTIMATED_ZIP_RATIO.get(filepath))
        else:
            logger.info("%s: text format detected.", os.path.basename(filepath))
        logger.info("%s: file size: %s bytes.", os.path.basename(filepath), os.path.getsize(filepath))
        phred_offset = detect_phred_offset(
            filepath=filepath,
            reads_for_phred_offset=parameters["reads_for_phred_offset"],
            phred_offset=parameters["phred_offset"]
        )
        logger.info("%s: Phred offset of %s detected.", os.path.basename(filepath), phred_offset)
        if ESTIMATED_READ_COUNTS.get(filepath) is not None:
            logger.info("%s: estimated bytes per read: %s.", os.path.basename(filepath), ESTIMATED_BYTE_PER_READ[filepath])
            logger.info("%s: (estimated) read count: %s.", os.path.basename(filepath), ESTIMATED_READ_COUNTS[filepath])
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
                    "%s: finished processing %d reads in %.2fs (%.0f pairs/sec).",
                    os.path.basename(filepath), total_reads, elapsed, total_reads / elapsed
                )
            else:
                logger.info("%s: finished processing in %.2fs.", os.path.basename(filepath), elapsed)
        else:
            raise TypeError(f"file {filepath} returned filetype {filetype}, which is not recognized.")

def process_unpaired_task_flat(task, parameters):
    """
    Worker wrapper for an unpaired task: parses the task's blob
    (parse_blob()), runs process_unpaired_chunk() on it and tags the result
    with the task type and filepath.

    Args:
        task (dict): An "unpaired" task from generate_unpaired_tasks().
        parameters (dict): Run parameters.

    Returns:
        tuple: (type, filepath, chunk_results, kept, rejected, rejected_reads),
            where chunk_results and rejected_reads are the kept and rejected
            record data, and kept and rejected are counts (see
            process_unpaired_chunk()).
    """
    chunk_results, kept, rejected, rejected_reads = process_unpaired_chunk(
        chunk=parse_blob(task["chunk"]),
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
    Worker wrapper for a paired task (two-file or interleaved): parses the
    task's blob(s) (parse_blob()), runs process_paired_chunk() on the R1/R2
    records and tags the result with the task type and filepaths. For an
    interleaved task, the records of its single blob are split alternately
    into R1 (records 1, 3, 5, ...) and R2 (records 2, 4, 6, ...).

    Args:
        task (dict): A "paired" task from generate_paired_tasks() or
            generate_unpaired_tasks(filetype="interleaved").
        parameters (dict): Run parameters.

    Returns:
        tuple: (type, file1, file2, paired_out_1, paired_out_2,
            R1_singles_out, R2_singles_out, num_paired, num_R1_singles,
            num_R2_singles, rejected_1, rejected_2, rejected_R1, rejected_R2).
            rejected_1/rejected_2 are rejected-read counts and
            rejected_R1/rejected_R2 the rejected record data, kept separate
            per mate since each mate is trimmed and filtered independently
            before reconciliation (see process_paired_chunk()).
    """
    file1 = task["file1"]
    file2 = task["file2"]
    if task.get("interleaved"):
        records = list(parse_blob(task["chunk1"]))
        chunk1, chunk2 = records[0::2], records[1::2]
    else:
        chunk1, chunk2 = parse_blob(task["chunk1"]), parse_blob(task["chunk2"])
    paired_out_1, paired_out_2, R1_singles_out, R2_singles_out, num_paired, num_R1_singles, num_R2_singles, rejected_1, rejected_2, rejected_R1, rejected_R2 = process_paired_chunk(
        chunks = (chunk1, chunk2),
        phred_offset_1=task["phred_offset_1"],
        phred_offset_2=task["phred_offset_2"],
        gzip_output = task["gzip_output"],
        gzip_level = task["gzip_level"],
        discard_singles=task["discard_singles"],
        parameters=parameters
    )
    return task["type"], file1, file2, paired_out_1, paired_out_2, R1_singles_out, R2_singles_out, num_paired, num_R1_singles, num_R2_singles, rejected_1, rejected_2, rejected_R1, rejected_R2

def generate_paired_tasks(files, chunk_size, parameters):
    """
    Lazily split paired-end file pairs into worker tasks.

    For each (R1, R2) pair, logs basic file information, detects the Phred
    offset of each file once, and then yields one task per chunk of
    ``chunk_size`` read pairs, as one raw FASTQ blob per file
    (lazy_fastq_blobs()), so that a single global pool can process chunks
    from many files. Mates are matched by position here; within a
    chunk they are reconciled by base read ID in process_paired_chunk().
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
            if GZIP_DETECTION[f]:
                logger.info("%s: gzip format detected.", os.path.basename(f))
                if ESTIMATED_ZIP_RATIO.get(f) is not None:
                    logger.info("%s: estimated gzip compression ratio: %s.", os.path.basename(f), ESTIMATED_ZIP_RATIO.get(f))
            else:
                logger.info("%s: text format detected.", os.path.basename(f))
            logger.info("%s: file size: %s bytes.", os.path.basename(f), os.path.getsize(f))
        
        phred_offset_1 = detect_phred_offset(filepath=file1, reads_for_phred_offset=parameters["reads_for_phred_offset"], phred_offset=parameters["phred_offset"])
        phred_offset_2 = detect_phred_offset(filepath=file2, reads_for_phred_offset=parameters["reads_for_phred_offset"], phred_offset=parameters["phred_offset"])
        logger.info("%s: Phred offset of %s detected.", os.path.basename(file1), phred_offset_1)
        logger.info("%s: Phred offset of %s detected.", os.path.basename(file2), phred_offset_2)
        logger.info("%s: (estimated) read count: %s.", os.path.basename(file1), ESTIMATED_READ_COUNTS.get(file1, "unknown"))
        logger.info("%s: (estimated) read count: %s.", os.path.basename(file2), ESTIMATED_READ_COUNTS.get(file2, "unknown"))
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

def prepare_reads(records, phred_offset, write_rejected, parameters):
    """
    Validate one mate's batch of FASTQ records and build the arrays every
    later step needs. Runs once per mate per chunk.

    Invalid records (see validate_fastq()) are rejected. If MGI conversion
    is enabled, the valid records' headers are converted to Illumina format
    and their plus lines reset to "+".

    Args:
        records (Iterable[tuple[bytes, bytes, bytes, bytes]]): (header,
            sequence, plus, quality) records, as produced by parse_blob().
        phred_offset (int): Phred encoding offset (33 or 64).
        write_rejected (bool): If True, also collect invalid records.
        parameters (dict): Run parameters (N-filtering, input length
            limits, MGI-to-Illumina header conversion).

    Returns:
        tuple[dict | None, int, list[bytes]]:
            - The batch, or None if no record passed validation. The batch
              is a dict with:
                - "valid_reads" (numpy.ndarray): index in ``records`` of the
                  record in each array row.
                - "headers", "sequences", "pluses", "qualities" (list[bytes]):
                  the valid records' fields, one entry per array row.
                - "sequence_arr", "quality_arr", "chunk_padding_bool",
                  "row_tilde_count", "padding_mask_bool", "raw_lengths":
                  outputs of seq_to_array() / qual_to_array().
            - Count of records rejected during validation.
            - Rejected records, populated only if write_rejected is True.
    """
    valid_reads = []
    valid_headers = []
    valid_sequences = []
    valid_pluses = []
    valid_qualities = []
    rejected_reads = []
    rejected = 0
    for k, (header, sequence, plus, quality) in enumerate(records):
        if validate_fastq(header, sequence, plus, quality, n_filter = parameters["n_filter"], min_length_input = parameters["min_length_input"], max_length_input = parameters["max_length_input"]):
            valid_reads.append(k)
            valid_headers.append(header)
            valid_sequences.append(sequence)
            valid_pluses.append(plus)
            valid_qualities.append(quality)
        else:
            rejected += 1
            if write_rejected:
                rejected_reads.append(b"\n".join((header, sequence, plus, quality)) + b"\n")
    if not valid_headers:
        return None, rejected, rejected_reads
    if parameters["mgi_convert_flag"]:
        valid_pluses = [b"+"] * len(valid_headers)
        valid_headers = [header_mgi_to_illumina(header, parameters["mgi_bc5"], parameters["mgi_bc7"], parameters["mgi_instrument"], parameters["mgi_run"]) for header in valid_headers]
    sequence_arr, chunk_padding_bool, row_tilde_count, padding_mask_bool, raw_lengths, max_len, n_reads = seq_to_array(sequence_list = valid_sequences)
    quality_arr = qual_to_array(quality_list = valid_qualities, phred_offset = phred_offset, chunk_padding_bool = chunk_padding_bool, padding_mask_bool = padding_mask_bool, max_len = max_len, n_reads = n_reads)
    batch = {
        "valid_reads": np.array(valid_reads),
        "headers": valid_headers,
        "sequences": valid_sequences,
        "pluses": valid_pluses,
        "qualities": valid_qualities,
        "sequence_arr": sequence_arr,
        "quality_arr": quality_arr,
        "chunk_padding_bool": chunk_padding_bool,
        "row_tilde_count": row_tilde_count,
        "padding_mask_bool": padding_mask_bool,
        "raw_lengths": raw_lengths,
    }
    return batch, rejected, rejected_reads


def paired_overlap_cutoffs(batch_1, batch_2, parameters):
    """
    Run adapter_trimming_overlap() on the reads of a chunk that are valid in
    both mates, and map the results back onto each mate's own batch rows.

    Args:
        batch_1 (dict): Output of prepare_reads() for R1 of the chunk.
        batch_2 (dict): Output of prepare_reads() for R2 of the same chunk
            (record k of one is the mate of record k of the other). Neither
            batch may be None.
        parameters (dict): Run parameters; uses "overlap_min_length",
            "overlap_mismatches" and "overlap_portion_mismatch".

    Returns:
        tuple: ((left_1, right_1), (left_2, right_2)), one entry per batch
            row. Rows without a valid mate get left 0 and right int32 max,
            so the other modules decide their trimming alone.
    """
    no_cut = np.iinfo(np.int32).max
    n_1, n_2 = len(batch_1["headers"]), len(batch_2["headers"])
    left_1 = np.zeros(n_1, dtype=np.int32)
    right_1 = np.full(n_1, no_cut, dtype=np.int32)
    left_2 = np.zeros(n_2, dtype=np.int32)
    right_2 = np.full(n_2, no_cut, dtype=np.int32)
    _, rows_1, rows_2 = np.intersect1d(batch_1["valid_reads"], batch_2["valid_reads"], assume_unique=True, return_indices=True)
    if rows_1.size == 0:
        return (left_1, right_1), (left_2, right_2)
    seq_arr_1 = batch_1["sequence_arr"][rows_1]
    seq_arr_2 = batch_2["sequence_arr"][rows_2]
    pad_1, pad_2 = batch_1["chunk_padding_bool"], batch_2["chunk_padding_bool"]
    tilde_1 = batch_1["row_tilde_count"][rows_1] if pad_1 else 0
    tilde_2 = batch_2["row_tilde_count"][rows_2] if pad_2 else 0
    (l1, r1), (l2, r2) = adapter_trimming_overlap(seq_arr_1, seq_arr_2, pad_1, tilde_1, pad_2, tilde_2, min_overlap = parameters["overlap_min_length"], max_mismatch = parameters["overlap_mismatches"], max_mismatch_frac = parameters["overlap_portion_mismatch"])
    left_1[rows_1], right_1[rows_1] = l1, r1
    left_2[rows_2], right_2[rows_2] = l2, r2
    return (left_1, right_1), (left_2, right_2)    

def finish_reads(batch, phred_offset, minimum_average_qual_post, min_length_output, max_length_output, min_length_output_perc, max_length_output_perc, write_rejected, overlap_cutoffs, parameters, read_direction=None):
    """
    Trim and length/quality-filter one prepared batch (one mate of a paired
    chunk) and format the survivors, keyed by their base (mate-independent)
    read ID.

    Runs all enabled modules from build_pipeline(), plus the overlap
    cutoffs if given, merges their windows (largest left, smallest right
    boundary), and applies the same keep rules as process_unpaired_chunk():
    non-empty window, within the output length limits, and (if set) mean
    quality of at least minimum_average_qual_post. Kept reads are cut to
    their window and optionally Phred re-encoded; rejected records are
    collected untrimmed.

    Args:
        batch (dict | None): Output of prepare_reads (None = nothing valid).
        phred_offset (int): Phred encoding offset (33 or 64).
        minimum_average_qual_post (float): Minimum acceptable mean quality
            after trimming. 0 disables this filter.
        min_length_output, max_length_output (int | None): Absolute length
            limits after trimming.
        min_length_output_perc, max_length_output_perc (int | None): Length
            limits as a percentage of input length; ignored if the absolute
            limit is given.
        write_rejected (bool): If True, also collect filtered-out records.
        overlap_cutoffs (tuple | None): (left, right) per batch row from
            paired_overlap_cutoffs, or None.
        parameters (dict): Run parameters (pipeline, Phred re-encoding).

    Returns:
        tuple[dict[bytes, bytes], int, list[bytes]]:
            - Survivors as {base read ID: formatted FASTQ record}, in batch
              order. If two survivors share a base ID, only the last is kept.
            - Count of records filtered out here.
            - Rejected records (populated only if write_rejected is True).
    """
    rejected_reads = []
    rejected = 0
    if batch is None:
        return {}, rejected, rejected_reads
    valid_headers, valid_sequences = batch["headers"], batch["sequences"]
    valid_pluses, valid_qualities = batch["pluses"], batch["qualities"]
    sequence_arr, quality_arr = batch["sequence_arr"], batch["quality_arr"]
    chunk_padding_bool, row_tilde_count = batch["chunk_padding_bool"], batch["row_tilde_count"]
    padding_mask_bool, raw_lengths = batch["padding_mask_bool"], batch["raw_lengths"]
    n_reads, length = sequence_arr.shape
    left_list = [np.zeros(n_reads, dtype=np.int8)]
    right_list = [np.full(n_reads, length, dtype=np.int32) - row_tilde_count]
    for step in build_pipeline(parameters, read_direction):
        left, right = step(sequence_arr, quality_arr, chunk_padding_bool, row_tilde_count, padding_mask_bool)
        left_list.append(left)
        right_list.append(right)
    if overlap_cutoffs is not None:
        overlap_left, overlap_right = overlap_cutoffs
        left_list.append(overlap_left)
        right_list.append(overlap_right)
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
        
def process_paired_chunk(chunks, phred_offset_1, phred_offset_2, gzip_output, gzip_level, discard_singles, parameters):
    """
    Trim and filter one chunk of paired reads, then sort the survivors into
    intact pairs and orphaned singletons.

    Both mates are validated (prepare_reads()), optionally overlap-trimmed
    as a pair (paired_overlap_cutoffs(), when "overlap_filter_flag" is set),
    and then trimmed and filtered independently (finish_reads()). Survivors
    are matched by base read ID: records whose mate also survived are output
    as pairs (in R1 order), the others as singletons.

    If parameters["interleaved_out"] is True (always the case with
    --stdout), pairs are written as alternating R1/R2 records into the first
    paired output, and both mates' rejected records are merged into the
    first rejected output; the second outputs are then empty.

    Args:
        chunks (tuple[Iterable, Iterable]): (chunk1, chunk2): R1 and R2
            records as (header, sequence, plus, quality) tuples (see
            parse_blob()), covering the same reads in the same order.
        phred_offset_1 (int): Phred encoding offset (33 or 64) for R1.
        phred_offset_2 (int): Phred encoding offset (33 or 64) for R2.
        gzip_output (bool): If True, compresses output records with isal gzip.
        gzip_level (int): isal gzip compression level (0-3).
        discard_singles (bool): If True, singleton records are dropped
            (their counts are still reported).
        parameters (dict): Run parameters, including "interleaved_out",
            "write_rejected", "overlap_filter_flag" and all module settings.

    Returns:
        tuple: An 11-element tuple containing:
            - paired_out_1 (bytes): R1 records of surviving pairs (all pair
              records, alternating R1/R2, when interleaving).
            - paired_out_2 (bytes): R2 records of surviving pairs (empty when
              interleaving).
            - singles_out_1 (bytes): Surviving R1 records whose mate did not
              survive (empty if discard_singles).
            - singles_out_2 (bytes): Surviving R2 records whose mate did not
              survive (empty if discard_singles).
            - num_paired (int): Count of surviving read pairs.
            - num_R1_singles (int): Count of surviving R1 singletons.
            - num_R2_singles (int): Count of surviving R2 singletons.
            - rejected_1 (int): Count of rejected R1 reads (invalid plus
              filtered out).
            - rejected_2 (int): Count of rejected R2 reads.
            - rejected_R1 (bytes): Rejected R1 records, untrimmed (both
              mates' rejected records when interleaving); empty unless
              "write_rejected" is set.
            - rejected_R2 (bytes): Rejected R2 records (empty when
              interleaving).
        Non-empty byte outputs are gzip-compressed when gzip_output is True.

    Raises:
        RuntimeError: If gzip compression fails.
    """
    chunk1, chunk2 = chunks
    batch_1, invalid_1, rejected_R1 = prepare_reads(chunk1, phred_offset_1, write_rejected = parameters["write_rejected"], parameters = parameters)
    batch_2, invalid_2, rejected_R2 = prepare_reads(chunk2, phred_offset_2, write_rejected = parameters["write_rejected"], parameters = parameters)
    overlap_1 = overlap_2 = None
    if parameters["overlap_filter_flag"]:
        overlap_1, overlap_2 = paired_overlap_cutoffs(batch_1, batch_2, parameters)
    survivors_1, filtered_1, filtered_R1 = finish_reads(batch_1, phred_offset_1, minimum_average_qual_post = parameters["minimum_average_qual_post"], min_length_output = parameters["min_length_output"], max_length_output = parameters["max_length_output"], min_length_output_perc = parameters["min_length_output_perc"], max_length_output_perc = parameters["max_length_output_perc"], write_rejected = parameters["write_rejected"], overlap_cutoffs = overlap_1, read_direction="read_1", parameters = parameters)
    survivors_2, filtered_2, filtered_R2 = finish_reads(batch_2, phred_offset_2, minimum_average_qual_post = parameters["minimum_average_qual_post"], min_length_output = parameters["min_length_output"], max_length_output = parameters["max_length_output"], min_length_output_perc = parameters["min_length_output_perc"], max_length_output_perc = parameters["max_length_output_perc"], write_rejected = parameters["write_rejected"], overlap_cutoffs = overlap_2, read_direction="read_2", parameters = parameters)
    rejected_1 = invalid_1 + filtered_1
    rejected_2 = invalid_2 + filtered_2
    rejected_R1 = rejected_R1 + filtered_R1
    rejected_R2 = rejected_R2 + filtered_R2

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
    
    if discard_singles:
        singles_out_1 = b""
        singles_out_2 = b""
        
    interleave = parameters["interleaved_out"]   # see point 4
    if interleave:
        paired_out_1 = b"".join(r1 + r2 for r1, r2 in zip(paired_out_1, paired_out_2))
        paired_out_2 = b""
        rejected_R1 = b"".join(rejected_R1) + b"".join(rejected_R2)
        rejected_R2 = b""
    else:
        paired_out_1 = b"".join(paired_out_1)
        paired_out_2 = b"".join(paired_out_2)
        rejected_R1 = b"".join(rejected_R1)
        rejected_R2 = b"".join(rejected_R2)
    singles_out_1 = b"".join(singles_out_1)
    singles_out_2 = b"".join(singles_out_2)

    if gzip_output:
        try:
            if paired_out_1:
                paired_out_1 = igzip.compress(paired_out_1, compresslevel=gzip_level)
            if paired_out_2:
                paired_out_2 = igzip.compress(paired_out_2, compresslevel=gzip_level)
            if rejected_R1:
                rejected_R1 = igzip.compress(rejected_R1, compresslevel=gzip_level)
            if rejected_R2:
                rejected_R2 = igzip.compress(rejected_R2, compresslevel=gzip_level)
            if singles_out_1:
                singles_out_1 = igzip.compress(singles_out_1, compresslevel=gzip_level)
            if singles_out_2:
                singles_out_2 = igzip.compress(singles_out_2, compresslevel=gzip_level)
        except Exception as e:
            raise RuntimeError(f"Compression failed inside worker: {str(e)}") from None
    return paired_out_1, paired_out_2, singles_out_1, singles_out_2, num_paired, num_R1_singles, num_R2_singles, rejected_1, rejected_2, rejected_R1, rejected_R2

##### Input handler functions #####
def worker_initilizer(parameters):
    """
    Pool initializer: runs once in each worker process at startup.

    Makes the worker ignore SIGINT, so Ctrl+C is handled by the main
    process only, and stores `parameters` in the WORKER_PARAMETERS global,
    so it is pickled and sent once per worker instead of with every task.

    Args:
        parameters (dict): Run parameters.
    """
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    global WORKER_PARAMETERS
    WORKER_PARAMETERS = parameters

def unified_worker(task):
    """
    Pool entry point: route a task to the paired or unpaired handler, using
    the run parameters stored by worker_initilizer().

    Args:
        task (dict): A task from generate_unpaired_tasks() or
            generate_paired_tasks(). Its "type" field ("paired" or
            "unpaired") selects the handler.

    Returns:
        tuple: The result of process_paired_task_flat() or
            process_unpaired_task_flat().

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

def input_handler(unspecified_files, unpaired_files, paired_files, interleaved_files, output_dir, threads, chunk_size, show_progress, stdout, interleaved_out, write_rejected, discard_singles, parameters):
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

    Returns:
        dict: Summary statistics per input. Unpaired files are keyed by
            their filepath, with "kept" and "rejected" counts. Paired and
            interleaved inputs are keyed by their name prefix, with
            "kept_pairs", "kept_R1_singletons", "kept_R2_singletons",
            "rejected_R1" and "rejected_R2" counts, since each mate is
            trimmed and filtered independently before reconciliation.
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
    if show_progress:    
        for file in unpaired:
            ESTIMATED_READ_COUNTS[file] = count_reads_estimated(file)
        for file in interleaved:
            ESTIMATED_READ_COUNTS[file] = count_reads_estimated(file)
        for f1, f2 in paired:
            ESTIMATED_READ_COUNTS[f1] = count_reads_estimated(f1)
            ESTIMATED_READ_COUNTS[f2] = count_reads_estimated(f2)
        total_reads = sum(ESTIMATED_READ_COUNTS.values())
        tracker = ProgressTracker(total_reads)
    else:
        class _NullTracker:
            def update(self, n): pass
            def close(self, interrupted=False): pass
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
    
    def open_pair_outputs(prefix):
        """Open every output file one paired/interleaved input needs (see pair_output_keys)."""
        for key in pair_output_keys(prefix).values():
            file_writing_handles[key] = open_fastq_writer(key, output_dir, gzip_output=parameters["gzip_output"])

    for file in interleaved:
        prefix = unique_prefix(basename_file(file))
        pair_keys[(file, file)] = prefix
        file_stats[prefix] = {"kept_pairs": 0, "kept_R1_singletons": 0, "kept_R2_singletons": 0, "rejected_R1": 0, "rejected_R2": 0}
        open_pair_outputs(prefix)
    for file1, file2 in paired:
        prefix = unique_prefix(common_name_parts([os.path.basename(file1), os.path.basename(file2)]))
        pair_keys[(file1, file2)] = prefix
        file_stats[prefix] = {"kept_pairs": 0, "kept_R1_singletons": 0, "kept_R2_singletons": 0, "rejected_R1": 0, "rejected_R2": 0}
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

    backpressure = threading.Semaphore(threads * 5)
    def bounded_chunk_stream():
        """
        Yield tasks from chunk_stream, blocking while threads * 5 tasks are
        still unprocessed (the semaphore is released for each result).
        """
        for item in chunk_stream:
            backpressure.acquire()
            yield item
            
    finished = False
    sys.setswitchinterval(0.0005)
    try:
        with mp.Pool(threads, initializer=worker_initilizer, initargs=(parameters,)) as pool:
            submit = pool.imap if parameters["ordered_output"] else pool.imap_unordered
            for result in submit(unified_worker, bounded_chunk_stream(), chunksize=1):
                backpressure.release()
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
                    tracker.update(num_paired * 2 + num_R1_singles + num_R2_singles + rejected_1 + rejected_2)
            finished = True
    finally:
        tracker.close(interrupted = not finished) 
        ACTIVE_PROGRESS_TRACKER = None
        for handle in file_writing_handles.values():
            if not handle.closed:
                if parameters["gzip_output"] and handle.tell() == 0:
                    handle.write(igzip.compress(b""))
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
    Print the help for --full-auto/-GO (called for --help combined with
    --full-auto): which options are still respected, and the options whose
    value full-auto mode changes (FULL_AUTO_OVERRIDES), grouped by the
    parser's argument groups. All other options use their regular defaults.

    Args:
        parser (argparse.ArgumentParser): The fully built Readzor parser.
    """
    print(
        "\n"
        "When --full-auto/-GO is specified, only input parameters --input-files/-i, --input-paired/-ip, --input-interleaved/-ii, and --input-unpaired/-iu are respected."
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
    Parse Readzor's command-line arguments into a resolved parameters dict.

    Input can come from any combination of --input-files (pairing is
    auto-detected), --input-paired, --input-unpaired and
    --input-interleaved. If none of these is given, FASTQ files
    (.fastq/.fq, optionally .gz/.gzip) in the current working directory are
    used with --full-auto; otherwise data piped into stdin is spooled to a
    temporary file (resolve_stdin_input()).

    With --full-auto, every option except the input options is reset to its
    default, apart from the values listed in FULL_AUTO_OVERRIDES.

    After parsing, some settings are resolved or adjusted:
      - --gzip turns --stdout off.
      - --stdout turns off verbose output, the progress bar and
        --write-rejected, and turns on --discard-singles and
        --interleaved-out.
      - --n-filter turns off N end trimming (there would be nothing left
        to trim).
      - With adapter trimming on, the adapter list is taken from
        --adapter-fasta-excl alone, or else from the selected
        --adapter-group(s) plus --adapter-fasta-add, and de-duplicated.
      - The thread count is resolved (worker_determination())

    Returns:
        dict: Run parameters, keyed one-for-one with the resolved CLI
            options (inputs, module flags and settings, output, threading
            and run-mode settings). See the "--- Store parameters ---" block
            in this function for the authoritative list of keys;
            "adapter_sequences" (list[bytes]) is added after that block.

    Side Effects:
        Sets up console logging. When run without any arguments, prints the
        help to stderr and exits with status 1. With --help (optionally
        combined with --full-auto), --version or --list-adapters, prints
        the corresponding output and exits. Calls parser.error() (which
        exits with status 2) if no input can be found. May create a
        temporary file holding stdin data.

    Raises:
        ValueError: If an adapter FASTA file is malformed, or
            --adapter-fasta-excl contains no sequences.
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
        help="Stream resulting FASTQ reads to stdout. Forces --interleaved-out and --discard-singles for paired and interleaved files. Overrides --verbose and --progress. Overridden to 'off' by --gzip. Note: all files will be streamed on end, without any seperators."
    )
    output_group.add_argument(
        "--interleaved-out", action="store_true", default=False,
        help="Interleave surviving paired FASTQ reads of paired and interleaved input files. Single reads are written to seperate files."
    )    
    output_group.add_argument(
        "--write-rejected", action="store_true", default=False,
        help="Write rejected reads to file. Either one (for unpaired and when --interleaved-out is set), or two (for forward and reverse reads) are produced. Overridden to 'off' when --stdout is set. Default: off."
    )  
    output_group.add_argument(
        "--discard-singles", action="store_true", default=False,
        help="Discard single leftover reads. For paired and interleaved reads, single surviving reads will be discarded instead of written to a seperate file. No effect on unpaired read filtering. Default: off."
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
                                                 "Trim reads for adapter sequences. Standard sequences included are TruSeq3, Nextera, and Illumina adapters. Adapter trimming is performed independent of quality.")
    adapter_trimming.add_argument(
        "--adapter-filter-flag", "-af", action="store_true", default = False,
        help='[FLAG] Turn on adapter trimming module. Default: off.'
    )
    adapter_trimming.add_argument(
        "--adapter-seed", "-as", type = int, default = 8, metavar="", choices=range(1,1000),
        help="Minimal 5'-end match length. Shorter seeds can result in more partial hits found at end of reads. Default: 8."
    )
    adapter_trimming.add_argument(
        "--adapter-group", "-ag", nargs = "+", default = ["TruSeq"], choices = [name for name, _ in DEFAULT_ADAPTERS], metavar="",
        help='Specify the group(s) of adapters to be used. Ignored if --adapter-fasta-excl is set. Choices: ' + ", ".join(["all"] + [name for name, _ in DEFAULT_ADAPTERS]) + '. Default: TruSeq.'
        )
    adapter_trimming.add_argument(
        "--adapter-mismatch", "-am", type = int, default = 0, metavar="",
        help="Number of mismatches allowed in adapter finding. Default: 0."
    )
    adapter_trimming.add_argument(
        "--adapter-fasta-add", "-ad", type = str, default = None, metavar="",
        help="Fasta file with adapter sequences to trim for, in addition to predefined sequences."
    )
    adapter_trimming.add_argument(
        "--adapter-fasta-excl", "-ax", type = str, default = None, metavar="",
        help="Fasta file with adapter sequences to trim for, excluding predefined and additional sequences specified."
    )
    

    overlap_trimming = parser.add_argument_group("Overlap trimming",
                                                 "Perform overlap analysis to find adapter sequences. Works independently of adapter sequence. Only available for paired and interleaved reads. Can be combined with --adapter-filter-flag.")
    overlap_trimming.add_argument(
        "--overlap-filter-flag", "-of", action="store_true", default = False,
        help="[FLAG] Turn on overlap trimming module. Default: off."
    )
    overlap_trimming.add_argument(
        "--overlap-mismatches", "-om", type = int, default = 5, metavar="",
        help='Maximum mismatches allowed in the overlapping region. Default: 5.'
    )
    overlap_trimming.add_argument(
        "--overlap-min-length", "-ol", type = int, default = 30, metavar="",
        help="Minimum overlap length between the paired reads. Default: 30."
    )
    overlap_trimming.add_argument(
        "--overlap-portion-mismatch", "-op", type = int, default = 10, metavar="",
        help="Maximum percentage of mismatched bases allowed in the overlapping region. Default: 5"
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
                                                  "Convert read header from MGI (BGI) format to Illumina format. Original header will be discarded.")
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
        "--chunk-size", type=int, default = 1500, metavar="",
        help="Number of reads per chunk sent to each worker. Changing can alter processing speed. Default: 1500. "
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
        if args.full_auto:
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
        elif not sys.stdin.isatty():
            args.input_files = [resolve_stdin_input(parser)]
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
    parameters["adapter_seed"] = args.adapter_seed
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
    parameters["adapter_group"] = (
        [name for name, _ in DEFAULT_ADAPTERS] if "all" in args.adapter_group
        else list(dict.fromkeys(args.adapter_group))
    )
    parameters["discard_singles"] = args.discard_singles
    parameters["overlap_filter_flag"] = args.overlap_filter_flag
    parameters["overlap_mismatches"] = args.overlap_mismatches
    parameters["overlap_min_length"] = args.overlap_min_length
    parameters["overlap_portion_mismatch"] = args.overlap_portion_mismatch
    
    if parameters["gzip_output"]:
        parameters["stdout"] = False
        
    if parameters["stdout"]:
        parameters["verbose"] = False
        parameters["show_progress"] = False
        parameters["write_rejected"] = False
        parameters["discard_singles"] = True
        parameters["interleaved_out"] = True
    
    if parameters["n_filter"]:
        parameters["n_trimming_flag"] = False

    if parameters.get("adapter_filter_flag"):
        if parameters.get("adapter_fasta_excl"):
            raw_adapters = load_adapters_from_fasta(parameters["adapter_fasta_excl"])
            if raw_adapters == []:
                raise ValueError(f"No sequences in file '{parameters['adapter_fasta_excl']}' detected.")
        else:
            selected_groups = set(parameters["adapter_group"])
            flat_default_adapters = [entry for name, entries in DEFAULT_ADAPTERS if name in selected_groups for entry in entries]
            if parameters.get("adapter_fasta_add"):
                raw_adapters = flat_default_adapters + load_adapters_from_fasta(parameters["adapter_fasta_add"])
            else:
                raw_adapters = flat_default_adapters
        parameters["adapter_sequences"] = list({seq.encode('utf-8') for _, seq in raw_adapters})
    else:
        parameters["adapter_sequences"] = []

    parameters["threads"] = worker_determination(parameters["threads"])
    
    return parameters

##### Wrap up functions #####
def write_summary_and_statistics(summary_results, parameters, output_dir):
    """
    Write per-input read counts to results_summary.txt in output_dir.

    The tab-separated file contains a "[Paired Reads]" table (one row per
    paired or interleaved input: name prefix, kept pairs, kept R1
    singletons, kept R2 singletons, rejected R1, rejected R2) and an
    "[Unpaired Reads]" table (filename, kept, rejected). Each table is only
    written if there are inputs of that kind. Run parameters are not
    written here; they are recorded in the log file (see log_parameters()).

    Args:
        summary_results (dict): Output of input_handler(): input name ->
            counts. Entries contain either ``kept`` and ``rejected``, or
            ``kept_pairs``, ``kept_R1_singletons``, ``kept_R2_singletons``,
            ``rejected_R1`` and ``rejected_R2`` for paired entries.
        parameters (dict): Run parameters (currently unused).
        output_dir (str): Directory in which to create the file.

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
        
def print_final_message(stdout):
    """
    Print Readzor's completion message (a citation request and a randomly
    chosen sign-off) to stdout, and write the same text to the log file
    only. Called at the end of a successful run.

    Args:
        stdout (bool): If True (reads are being streamed to stdout),
            nothing is printed or logged, to keep the FASTQ stream clean.

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
    Readzor entry point.

    Selects the multiprocessing start method (fork if available, otherwise
    forkserver, otherwise spawn) and parses the command line. With
    --testrun, hands over to test_run(). Otherwise creates the timestamped
    results folder, starts file logging, logs all parameters, processes
    all inputs (input_handler()), writes results_summary.txt and prints the
    completion message. On Ctrl+C, exits with status 130. Temporary stdin
    files are always removed.
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
        summary_results = input_handler(unspecified_files = parameters["unspecified_files"], unpaired_files = parameters["unpaired_files"], paired_files = parameters["paired_files"], interleaved_files = parameters["interleaved_files"], output_dir = created_output_dir, threads = parameters["threads"], chunk_size = parameters["chunk_size"], show_progress = parameters["show_progress"], stdout = parameters["stdout"], interleaved_out = parameters["interleaved_out"], discard_singles = parameters["discard_singles"], write_rejected = parameters["write_rejected"], parameters = parameters)
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
    Run the full pipeline on a small sample of the input to validate the
    settings, without keeping any output. Called when --testrun is set.

    Works like a normal run, except that the results folder (including the
    log file and summary) is created inside a temporary directory that is
    deleted afterwards, console output is always verbose, and only the
    first ``threads`` chunks of each input are processed (one chunk per
    worker). Temporary stdin files are always removed.

    Args:
        parameters (dict): Run parameters.
    """
    try:
        with tempfile.TemporaryDirectory(prefix="readzor_testrun_") as tmp_dir:
            created_output_dir = create_folder_structure(tmp_dir)
            setup_logging(output_dir = created_output_dir, verbose = True, parameters = parameters)
            log_parameters(parameters)
            summary_results = input_handler(unspecified_files = parameters["unspecified_files"], unpaired_files = parameters["unpaired_files"], paired_files = parameters["paired_files"], interleaved_files = parameters["interleaved_files"], output_dir = created_output_dir, threads = parameters["threads"], chunk_size = parameters["chunk_size"], show_progress = parameters["show_progress"], stdout = parameters["stdout"], interleaved_out = parameters["interleaved_out"],discard_singles = parameters["discard_singles"], write_rejected = parameters["write_rejected"], parameters = parameters)
            write_summary_and_statistics(summary_results, parameters, output_dir = created_output_dir)
            logger.info("All files deleted.")
            logger.info("Testrun completed!")
    finally:
        cleanup_stdin_temp_files()
    
if __name__ == "__main__":
    main()