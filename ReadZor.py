#Ideas:
    # 5* in worker determination still neded on slurm now that stream is continuous?
    # check barcodes
    # logging
    # if --gzip-output; compress and write chunks on workers. Later; concatenate chunks. n chunks = n workers, maybe 3x for paired files.
    
#!/usr/bin/env python3

##### Import packages #####
import argparse
from collections import defaultdict
import datetime
import gzip
import io
import itertools
from multiprocessing import Pool
import os
import random
import re
import shlex
import shutil
import sys
import time

import numpy as np

##### Definition of constant values #####
VERSION = "0.0.3"
STRICT_NUCLEOTIDE_REGEX = re.compile(r'^[ATCG]+$')
LENIENT_NUCLEOTIDE_REGEX = re.compile(r'^[ATCGN]+$')
PHRED_REGEX = re.compile(r'^[!-~]+$')
DEFAULT_ADAPTERS = [
    ("TruSeq3", "AGATCGGAAGAGC"), #12x in human genome
    ("TruSeq3_R1", "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"), 
    ("TruSeq3_R2", "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"),
    ("Nextera", "CTGTCTCTTATACACATCT"),
    ("Nextera_transposas_R1", "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"),
    ("Nextera_transposas_R2", "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"),
    ("Nextera_PCR_i7", "GTCTCGTGGGCTCGG"),
    ("Nextera_PCR_i5", "TCGTCGGCAGCGTC"),
    ("Illumina_RNA","ACTGTCTCTTATACACATCT"),
    ("TruSeq_small_RNA","TGGAATTCTCGGGTGCCAAGG")
]
FULL_AUTO_PRESERVED_DESTS = {"input_files", "input_paired", "input_unpaired", "full_auto"}
FULL_AUTO_OVERRIDES = {
    "endqual_filter_flag": True,
    "adapter_trim_flag": True,
    "nucl_filter": True,
    "progress": True,
    "min_length": 100,
}

##### Progress tracker #####
def estimate_bytes_per_read(filepath):
    """
    Estimates the average on-disk (uncompressed) bytes consumed by one FASTQ
    record, by sampling the first record's header, sequence, plus-line, and
    quality lengths and adding back the newline stripped by `lazy_fastq`.

    Args:
        filepath (str): Path to the FASTQ file (.fastq, .fq, or gzip-compressed).

    Returns:
        int: Estimated number of bytes one record occupies on disk,
            including line-ending characters.

    Raises:
        ValueError: If the file contains no readable FASTQ records.
    """
    for header, sequence, plus, quality in lazy_fastq(filepath):
        header_bytes = len(header.encode('utf-8')) + 1
        seq_bytes = len(sequence.encode('utf-8')) + 1 
        plus_bytes = len(plus.encode('utf-8')) + 1
        qual_bytes = len(quality.encode('utf-8')) + 1
        return header_bytes + seq_bytes + plus_bytes + qual_bytes
    raise ValueError(f"No FASTQ records found in '{filepath}'; cannot estimate bytes per read.")

def estimate_gzip_ratio(filepath, sample_bytes=4 * 1024 * 1024):
    """
    Estimate a gzip file's compression ratio (uncompressed / compressed) by
    decompressing a leading sample, rather than the whole file.
 
    Returns None if the sample is too small to give a stable ratio
    (e.g. the whole file is tiny) -- caller should fall back to a fixed
    default ratio in that case.
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
 
def count_reads_estimated(filepath, default_gzip_ratio=4):
    """
    Estimate the number of reads in a FASTQ file from file size and a
    per-read byte estimate, without a full parse.
 
    For gzip files, estimates the uncompressed size via a sampled
    compression ratio (falls back to `default_gzip_ratio` if the file is
    too small to sample reliably).
    """
    file_size = os.path.getsize(filepath)
    bytes_per_read = estimate_bytes_per_read(filepath)
 
    if is_gz_file(filepath):
        ratio = estimate_gzip_ratio(filepath)
        if ratio is None:
            ratio = default_gzip_ratio
        estimated_uncompressed_size = file_size * ratio
    else:
        estimated_uncompressed_size = file_size
 
    return max(1, round(estimated_uncompressed_size / bytes_per_read))

class ProgressTracker:
    """
    Tracks completed reads against a known total and renders a single-line
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
        """Formats seconds into human-readable duration strings."""
        if seconds_val == float('inf') or seconds_val < 0 or seconds_val is None:
            return "?"
        
        total_seconds = int(round(seconds_val))
        hours, remainder = divmod(total_seconds, 3600)
        minutes, seconds = divmod(remainder, 60)

        # Concise output for live ETA (e.g., '1h 05m 12s', '05m 12s', '12s')
        if concise:
            if hours > 0:
                return f"{hours}h {minutes:02d}m {seconds:02d}s"
            if minutes > 0:
                return f"{minutes}m {seconds:02d}s"
            return f"{seconds}s"

        # Verbose output for completion message
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
        if self.done >= self.total:
            return
        if now - self._last_render >= self.min_interval:
            self._render()
            self._last_render = now

    def _render(self):
        term_width = shutil.get_terminal_size(fallback=(80, 24)).columns

        frac = min(self.done / self.total, 1.0)
        elapsed = time.time() - self._start_time
        rate = self.done / elapsed if elapsed > 0 else 0
        
        # Calculate ETA
        eta = (self.total - self.done) / rate if rate > 0 else float('inf')
        eta_str = self._format_duration(eta, concise=True)

        # Construct status text with formatted ETA
        stats = (
            f" {frac*100:5.1f}% "
            f"({self.done:,}/{self.total:,} reads). "
            f"Rate: {rate:,.0f} reads/s. Estimated time remaining: {eta_str}"
        )

        max_bar_len = term_width - len(stats) - 3

        if max_bar_len >= 5:
            effective_bar_width = min(self.bar_width, max_bar_len)
            filled = int(effective_bar_width * frac)
            bar = "#" * filled + "-" * (effective_bar_width - filled)
            line = f"\r[{bar}]{stats}"
        else:
            line = f"\r{stats.strip()}"

        # Clamp and pad line to prevent line-wrapping
        line = line[: term_width - 1].ljust(term_width - 1)
        sys.stderr.write(line)
        sys.stderr.flush()

    def close(self):
        term_width = shutil.get_terminal_size(fallback=(80, 24)).columns
        elapsed = time.time() - self._start_time
        rate = self.done / elapsed if elapsed > 0 else 0
        time_str = self._format_duration(elapsed, concise=False)

        stats = (
            f" 100% ({self.done:,} reads analyzed). "
            f"Average rate: {rate:,.0f} reads/s. Total time: {time_str}."
        )

        max_bar_len = term_width - len(stats) - 3
        if max_bar_len >= 5:
            effective_bar_width = min(self.bar_width, max_bar_len)
            bar = "#" * effective_bar_width
            line = f"\r[{bar}]{stats}"
        else:
            line = f"\r{stats.strip()}"

        line = line[: term_width - 1].ljust(term_width - 1)
        sys.stderr.write(line + "\n")
        sys.stderr.flush()
        
class _NullTracker:
    """No-op progress tracker used as a drop-in stand-in for a real progress
    bar when progress display is disabled. Mirrors the minimal
    interface (`update`, `close`) so calling code doesn't need conditional
    branches to check whether tracking is active.
    """
    def update(self, n): pass
    def close(self): pass

##### Helper functions #####
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
    
    Generates a new subfolder named "ReadZor_results_<YYYY-MM-DD_HH-MM-SS>"
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
    created_output_dir = os.path.join(output_dir, "ReadZor_results_" + datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S"))
    os.makedirs(created_output_dir, exist_ok=True)
    return created_output_dir

def worker_determination(threads):
    """
    Determine the number of worker processes for parallel task execution.
    
    The worker count is calculated based on the execution environment (Slurm vs.
    local) and whether the user provided an explicit thread request:
    
    - **On a Slurm-managed cluster** (SLURM_CPUS_PER_TASK detected):
        - If `threads` is an int: returns `threads * 10`.
        - Otherwise: returns `SLURM_CPUS_PER_TASK * 10`.
    - **On a local/non-Slurm system**:
        - If `threads` is an int: returns `threads` as-is.
        - Otherwise: returns `(available CPUs - 1)`, minimum of 1.
    
    All return values are clamped to a minimum of 1 worker.

    Args:
        threads (int | None): User-requested number of base threads/CPUs.
            Ignored if not a valid int (e.g., None, float, bool). Defaults to None.

    Returns:
        int: The number of worker processes to spawn, always >= 1.
    
    Notes:
        - On Slurm systems, the 10x multiplier allows for efficient handling of
          I/O-bound operations alongside CPU-bound work.
        - Boolean values are explicitly excluded from int detection to prevent
          True/False being treated as 1/0.
        - If `os.cpu_count()` returns None, defaults to 1 CPU.
    """
    slurm_cpus = os.environ.get('SLURM_CPUS_PER_TASK')
    threads_is_int = isinstance(threads, int) and not isinstance(threads, bool)
    if slurm_cpus is not None:
        if threads_is_int:
            return max(1, threads * 10)
        else:
            try:
                return max(1, int(slurm_cpus) * 10)
            except ValueError:
                available_cpu = os.cpu_count() or 1
                return max(1, available_cpu - 1)
    else:
        if threads_is_int:
            return max(1, threads)
        else:
            available_cpu = os.cpu_count() or 1
            return max(1, available_cpu - 1)

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

def query_read_length(filepath):
    """
    Retrieve the sequence length of the first read in a FASTQ file.
    
    Lazily reads the file and returns the length of the first read's sequence,
    assuming all reads in the file have uniform length. Stops reading after
    the first record, making this efficient even for large files.

    Args:
        filepath (str): Path to the FASTQ file to inspect.

    Returns:
        int | None: The sequence length of the first read, or None if the file
        is empty or contains no valid records.
    """
    for _, sequence, _, _ in lazy_fastq(filepath):
        return len(sequence)
    return None

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
    
    Detects gzip compression via magic bytes (independent of file extension),
    uses a 10 MB read buffer for efficient I/O, and yields each 4-line FASTQ
    record as (header, sequence, plus, quality) tuples with line endings stripped.

    Args:
        filepath (str): Path to the FASTQ file (.fastq, .fq, or gzip-compressed).

    Yields:
        tuple[str, str, str, str]: (header, sequence, plus, quality) for each read.
    """
    if is_gz_file(filepath):
        raw = open(filepath, 'rb', buffering=10 * 1024 * 1024)
        try:
            fastq_file = io.TextIOWrapper(gzip.GzipFile(fileobj=raw), encoding='utf-8')
        except Exception:
            raw.close()
            raise
    else:
        fastq_file = open(filepath, 'r', buffering=10 * 1024 * 1024, encoding='utf-8')
    with fastq_file:
        lines = iter(fastq_file)
        for header in lines:
            header = header.strip()
            sequence = next(lines).strip()
            plus = next(lines).strip()
            quality = next(lines).strip()
            yield header, sequence, plus, quality

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
        tuple[list[tuple[str, str]], list[str]]:
            - A list of (file_1, file_2) tuples representing matched pairs (ordered R1, R2).
            - A list of filepaths that could not be matched with a valid pair.
    """
    base_ids = {}
    for filepath in filepaths:
        try:
            if is_gz_file(filepath):
                raw = open(filepath, 'rb', buffering=1024)
                try:
                    file = gzip.GzipFile(fileobj=raw)
                    with file:
                        first_line = file.readline()
                finally:
                    raw.close()
            else:
                with open(filepath, 'rb', buffering=1024) as file:
                    first_line = file.readline()
        except (IOError, OSError):
            continue
        header = first_line.decode('utf-8', errors='replace').strip().lstrip('@')
        if not header:
            continue
        base_id, read_num = read_info_from_header(header)
        base_ids[filepath] = (base_id, read_num)
    pairs_dict = defaultdict(list)
    for filepath, (base_id, read_num) in base_ids.items():
        pairs_dict[base_id].append((filepath, read_num))
    pairs = []
    used = set()
    for base_id, file_list in pairs_dict.items():
        if len(file_list) == 2:
            (file_1, read_1), (file_2, read_2) = file_list
            if read_1 == 1 and read_2 == 2:
                pairs.append((file_1, file_2))
            elif read_1 == 2 and read_2 == 1:
                pairs.append((file_2, file_1))
            used.add(file_1)
            used.add(file_2)
        elif len(file_list) > 2:
            print(
                    f"Warning: Found {len(file_list)} files matching base ID '{base_id}' "
                    f"(expected max 2 for paired-end data). Skipping automatic pairing for these files.",
                    file=sys.stderr
                )
    unpaired = [f for f in base_ids if f not in used]
    return pairs, unpaired

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
    if ' ' in header:
        parts = header.split(' ', 1)
        base_id = parts[0]
        m = re.match(r'([12]):', parts[1])
        if m:
            read_num = int(m.group(1))
    else:
        m = re.search(r'/([12])$', header)
        if m:
            read_num = int(m.group(1))
            base_id = re.sub(r'/[12]$', '', header)
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
        for _, _, _, quality in reader:
            if count >= reads_for_phred_offset:
                break
            q_bytes = np.frombuffer(quality.encode('utf-8'), dtype=np.uint8)
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

def validate_fastq(header, sequence, plus, quality, min_raw_read_length, nucleotide_regex, read_length):
    """
    Validates that a single FASTQ record is well-formed.

    Checks that the header starts with "@" and has content, the plus-line 
    starts with "+" (and optionally matches the header description), the 
    sequence and quality strings are of equal length, meet a minimum length 
    requirement, and match the expected nucleotide and Phred character sets. 
    Also enforces an exact expected read length if specified.

    Args:
        header (str): The FASTQ header line (including leading "@").
        sequence (str): The nucleotide sequence line.
        plus (str): The separator line (starting with "+").
        quality (str): The Phred quality string.
        min_raw_read_length (int): Minimum acceptable sequence length.
        nucleotide_regex (re.Pattern): Compiled regex for valid nucleotides.
        read_length (int | None): Exact expected read length, if enforced.

    Returns:
        bool: True if the record is valid, False otherwise.

    Raises:
        ValueError: If a read length mismatch occurs against the expected length,
            or if the quality string length does not match the sequence length.
    """
    if header[:1] != "@":
        return False
    if len(header) <= 1 or not header[1:].strip():
        return False
    if plus[:1] != '+':
        return False
    if len(plus) > 1 and plus[1:] != header[1:]:
        return False
    seq_len = len(sequence)
    if read_length is not None and seq_len != read_length:
        raise ValueError(f"FASTQ file contain uneven read lengths. Found {seq_len}, expected {read_length}.")
    if seq_len < min_raw_read_length or seq_len != len(quality):
        return False
    if not nucleotide_regex.match(sequence):
        return False
    if not PHRED_REGEX.match(quality):
        return False
    if not len(quality) == seq_len:
        raise ValueError(f"FASTQ file malformed. Found quality line length of {len(quality)}, expected {read_length}.")
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
                output_sequence.append(sequence)
                sequence = fastafile.readline().rstrip()
            if len(output_sequence) == 0:
                raise ValueError(f"No sequence detected for header {header}")
            name = header.lstrip(">").strip()
            joined_sequence = "".join(output_sequence)
            result.append((name, joined_sequence))
            header = sequence
        return result

def chunk_size_setter(chunk_size):
    """
    Resolves the chunk size to use for parallel processing.

    If a chunk size is explicitly given, it's returned unchanged.
    Otherwise, detects whether Slurm tools are available on the system
    (via `sinfo`/`sbatch` on PATH) and picks a larger default chunk size
    for Slurm-managed (typically higher-resource) systems, or a smaller
    default otherwise.

    Args:
        chunk_size (int | None): User-specified chunk size. Defaults to
            None (auto-detect).

    Returns:
        int: The resolved chunk size — 20000 on Slurm systems, 1000
            otherwise, unless overridden.
    """
    if chunk_size is not None:
        return chunk_size
    if shutil.which('sinfo') is not None or shutil.which('sbatch') is not None:
        chunk_size = 20000
    else:
        chunk_size = 1000
    return chunk_size

def qual_to_bin(quality_list, phred_offset):
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
    joined = ''.join(quality_list).encode('utf-8')
    array = np.frombuffer(joined, dtype=np.uint8).astype(np.int8).reshape(len(quality_list), len(quality_list[0])) - phred_offset
    return array

def seq_to_bin(sequence_list):
    """
    Converts a list of nucleotide sequence strings into a 2D numpy array.

    Concatenates all sequences and reinterprets the raw bytes as an
    array of ASCII codes, reshaping into a (n_reads, read_length) matrix 
    for vectorized downstream processing.

    Args:
        sequence_list (list[str]): Sequence strings, all of equal length.

    Returns:
        numpy.ndarray: Signed 8-bit integer array of shape
            (len(sequence_list), read_length) of ASCII character codes.
    """
    joined = ''.join(sequence_list).encode('utf-8')
    array = np.frombuffer(joined, dtype = np.uint8).astype(np.int8).reshape(len(sequence_list), len(sequence_list[0]))
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
    mgi_header = mgi_header.lstrip("@").strip()
    strings = re.search(r"^(\w+)L(\d+)C(\d+)R(\d{3})(\d+)\/([12])", mgi_header)
    if strings is None:
        raise ValueError(f"Header does not match expected MGI format: {mgi_header!r}")
    illumina_header = "@%s:%s:%s:%d:%d:%d:%d %d:N:0:%s+%s"%(
    instrument, run, strings.group(1), int(strings.group(2)), int(strings.group(5)), int(strings.group(3)), int(strings.group(4)), int(strings.group(6)), barcode5, barcode7)
    return illumina_header
    
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
    out_filepath = os.path.join(output_dir, basename_for_write + extension)
    return open(out_filepath, "xb")

##### Processing reads functions #####
def trim_ends_quality(quality_arr, min_quality_both, endqual_min_start, endqual_min_end, endqual_filter_flag):
    """
    Determines per-read trim boundaries based on quality thresholds at each end.

    For each read, finds the position meeting `endqual_min_start` at the start,
    and the position meeting `endqual_min_end` at the end. Vectorized
    across all reads in a chunk. Reads with no position meeting the threshold
    are assigned zero-length cutoffs. If `endqual_filter_flag` is False,
    returns full-length arrays without trimming.

    Args:
        quality_arr (numpy.ndarray): (n_reads, read_length) array of per-base quality scores.
        min_quality_both (int | None): Minimum quality fallback to keep bases from both ends.
        endqual_min_start (int | None): Minimum quality to keep bases from the start of the read.
        endqual_min_end (int | None): Minimum quality to keep bases from the end of the read.
        endqual_filter_flag (bool): Flag to enable or disable end quality trimming.

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: A tuple of `(start_cutoffs, end_cutoffs)` 
            arrays of shape `(n_reads,)` and dtype `int32`, giving the left and right trim 
            boundaries per read.
    """
    n_reads, length = quality_arr.shape
    if not endqual_filter_flag:
        return np.zeros(n_reads, dtype = np.int32), np.full(n_reads, length, dtype = np.int32)
    if endqual_min_start is None:
        endqual_min_start = min_quality_both if min_quality_both is not None else 0
    if endqual_min_end is None:
        endqual_min_end = min_quality_both if min_quality_both is not None else 0
    qual_mask = quality_arr >= endqual_min_start
    start_good_pos = qual_mask.any(axis=1)
    start_cutoffs = qual_mask.argmax(axis=1)
    start_cutoffs = np.where(start_good_pos, start_cutoffs, 0)
    if endqual_min_end == endqual_min_start:
        end_cutoffs = length - qual_mask[:, ::-1].argmax(axis=1)
        end_cutoffs = np.where(start_good_pos, end_cutoffs, 0)
        return start_cutoffs.astype(np.int32), end_cutoffs.astype(np.int32)
    else:
        qual_mask = quality_arr >= endqual_min_end
        any_good_right = qual_mask.any(axis=1)
        end_cutoffs = length - qual_mask[:, ::-1].argmax(axis=1)
        end_cutoffs = np.where(any_good_right, end_cutoffs, 0)
        return start_cutoffs.astype(np.int32), end_cutoffs.astype(np.int32)
        
def homopolymer_nucleotide_trimming(sequence_arr, poly_length_both, poly_length_start, poly_length_end, poly_bases_both, poly_bases_start, poly_bases_end, poly_filter_flag):
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
        poly_filter_flag (bool): Flag to enable or disable homopolymer trimming.

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: (left_cutoffs, right_cutoffs),
            each of shape (n_reads,) and dtype int32, giving the left and right 
            trim boundaries per read.
    """
    n_reads, length = sequence_arr.shape
    if not poly_filter_flag:
        return np.zeros(n_reads, dtype=np.int32), np.full(n_reads, length, dtype=np.int32)
    if not poly_bases_start and not poly_bases_end and not poly_bases_both:
        return np.zeros(n_reads, dtype=np.int32), np.full(n_reads, length, dtype=np.int32)
    if poly_length_both == poly_length_start == poly_length_end == 0:
        return np.zeros(n_reads, dtype=np.int32), np.full(n_reads, length, dtype=np.int32)
    
    start_bases = []
    end_bases = []
    if poly_bases_both is not None:
        bases = [b.strip() for b in poly_bases_both.split(",") if b.strip()]
        start_bases = bases
        end_bases = bases
    if poly_bases_start is not None:
        start_bases = [b.strip() for b in poly_bases_start.split(",") if b.strip()]
    if poly_bases_end is not None:
        end_bases = [b.strip() for b in poly_bases_end.split(",") if b.strip()]
    
    poly_length_start = poly_length_start if poly_length_start != 0 else poly_length_both
    poly_length_end = poly_length_end if poly_length_end != 0 else poly_length_both
    
    right_cutoffs = np.full(n_reads, length, dtype=np.int32)
    left_cutoffs = np.zeros(n_reads, dtype=np.int32)
    
    for base in start_bases:
            base_code = ord(base)
            non_base_mask = sequence_arr != base_code
            padded_mask = np.column_stack([non_base_mask, np.ones(n_reads, dtype=bool)])
            first_non_pos = padded_mask.argmax(axis=1)
            trim_amount = np.where(first_non_pos >= poly_length_start, first_non_pos, 0)
            left_cutoffs = np.maximum(left_cutoffs, trim_amount)
    
    for base in end_bases:
            base_code = ord(base)
            rev_seq = sequence_arr[:, ::-1]
            non_base_mask = rev_seq != base_code
            padded_mask = np.column_stack([non_base_mask, np.ones(n_reads, dtype=bool)])
            first_non_pos = padded_mask.argmax(axis=1)
            trim_amount = np.where(first_non_pos >= poly_length_end, first_non_pos, 0)
            base_right_cutoffs = length - trim_amount
            right_cutoffs = np.minimum(right_cutoffs, base_right_cutoffs)
    
    return left_cutoffs, right_cutoffs

def n_end_trimming(sequence_arr, n_trimming_flag):
    """
    Removes leading and trailing N-bases from each read.

    This function detects and strips runs of N's from both the ends
    of each read (any length >= 1), leaving internal N's untouched.
    Internally reuses `homopolymer_nucleotide_trimming` under the hood
    configured for N bases.

    Args:
        sequence_arr (numpy.ndarray): (n_reads, read_length) array of
            per-base ASCII sequence codes.
        n_trimming_flag (bool): Whether N-end trimming is enabled. If False, returns
            boundaries that trim nothing (left=0, right=read_length) for
            every read.

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: (left_cutoffs, right_cutoffs),
            each of shape (n_reads,) and dtype int32, giving the trim
            boundaries per read.
    """
    n_reads, sequence_length = sequence_arr.shape
    if not n_trimming_flag: 
        return np.zeros(n_reads, dtype=np.int32), np.full(n_reads, sequence_length, dtype=np.int32)
    lefts, rights = homopolymer_nucleotide_trimming(sequence_arr, poly_length_both = 1, poly_length_start = 0, poly_length_end = 0, poly_bases_both = "N", poly_bases_start = None, poly_bases_end = None, poly_filter_flag = True)
    return lefts, rights

def cut_set_ends(sequence_arr, cut_both, cut_start, cut_end, cut_flag):
    """
    Produces fixed, user-specified trim boundaries applied uniformly to
    every read (e.g. hard-trimming a known number of adapter/primer bases
    off each end, regardless of quality).

    `cut_both` sets a symmetric default trim for both ends, but is overridden
    on either side individually if `cut_start` and/or `cut_end` are also
    given — so a user can specify a general trim amount while still
    customizing one end specifically. If `cut_flag` is False, returns
    boundaries that trim nothing (left=0, right=read_length) for every read.

    Args:
        sequence_arr (numpy.ndarray): (n_reads, read_length) array of
            per-base ASCII sequence codes.
        cut_both (int | None): Number of bases to trim off both ends.
            Used as a fallback for any side not given explicitly via
            cut_start/cut_end.
        cut_start (int): Number of bases to trim from the 5'
            end. Takes priority over cut_both if given.
        cut_end (int): Number of bases to trim from the 3'
            end. Takes priority over cut_both if given.
        cut_flag (bool): Flag to enable or disable fixed end cutting.

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: (left_cutoffs, right_cutoffs),
            each of shape (n_reads,) and dtype int32, identical across all
            reads, giving the left and right trim boundaries per read.
    """

    n_reads, length = sequence_arr.shape
    if not cut_flag:
        return np.zeros(n_reads, dtype = np.int32), np.full(n_reads, length, dtype = np.int32)
    cut_start = cut_start if cut_start != 0 else cut_both
    cut_end = cut_end if cut_end != 0 else cut_both
    cut_end = length - cut_end
    if cut_start > cut_end:
        cut_start = cut_end
    return np.full(n_reads, cut_start, dtype=np.int32), np.full(n_reads, cut_end, dtype=np.int32)

def sliding_window_quality(quality_arr, slider_quality, slider_window, slider_step, slider_filter_flag):
    """
    Determines per-read trim boundaries using a sliding-window quality scan,
    finding the longest (and highest-average-quality, as tiebreaker) stretch
    of the read where every window of `slider_window` bases has a mean
    quality above `slider_quality`.

    Rather than trimming only from the ends inward, it
    identifies the best surviving internal stretch of acceptable quality
    and reports its boundaries. If `slider_filter_flag` is False,
    returns boundaries that trim nothing (left=0, right=read_length) for
    every read.

    Args:
        quality_arr (numpy.ndarray): (n_reads, read_length) array of
            per-base quality scores.
        slider_quality (int): Minimum acceptable mean quality within a window.
        slider_window (int): Number of bases per sliding window.
        slider_step (int): Step size between successive window start positions.
        slider_filter_flag (bool): Flag to enable or disable sliding window quality filtering.

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: (left_cutoffs, right_cutoffs),
            each of shape (n_reads,) and dtype int32, giving the best surviving 
            [left, right) region per read. Reads with no failing windows keep
            their full length; reads that fail everywhere get a zero-length region.
    """
    n_reads, length = quality_arr.shape
    if not slider_filter_flag:
        return np.zeros(n_reads, dtype = np.int32), np.full(n_reads, length, dtype = np.int32)
    if length < slider_window:
        return np.zeros(n_reads, dtype = np.int32), np.full(n_reads, length, dtype = np.int32)

    cumsum = np.cumsum(quality_arr, axis=1, dtype=np.int32)
    cumsum = np.concatenate([np.zeros((n_reads, 1), dtype=np.int32), cumsum], axis=1)

    window_starts = np.arange(0, length - slider_window + 1, slider_step)
    window_sums = cumsum[:, window_starts + slider_window] - cumsum[:, window_starts]
    window_means = window_sums / slider_window
    failed_mask = window_means <= slider_quality

    bad_positions = np.zeros((n_reads, length), dtype=bool)
    for j, start in enumerate(window_starts):
        rows_failed = failed_mask[:, j]
        bad_positions[rows_failed, start:start + slider_window] = True

    good_positions = ~bad_positions
    no_bad = ~bad_positions.any(axis=1)
    all_bad = bad_positions.all(axis=1)

    left_cutoffs = np.zeros(n_reads, dtype=np.int32)
    right_cutoffs = np.zeros(n_reads, dtype=np.int32)
    right_cutoffs[no_bad] = length

    needs_stretch_search = ~no_bad & ~all_bad
    if not needs_stretch_search.any():
        return left_cutoffs, right_cutoffs

    sub = good_positions[needs_stretch_search]
    m = sub.shape[0]
    padded = np.zeros((m, length + 2), dtype=bool)
    padded[:, 1:-1] = sub
    diffs = np.diff(padded.astype(np.int8), axis=1)

    row_idx, start_cols = np.where(diffs == 1)
    _, end_cols = np.where(diffs == -1)

    run_lengths = end_cols - start_cols
    row_cumsum = cumsum[needs_stretch_search]
    run_sums = row_cumsum[row_idx, end_cols] - row_cumsum[row_idx, start_cols]
    run_means = run_sums / run_lengths

    order = np.lexsort((-run_means, -run_lengths, row_idx))
    sorted_row_idx = row_idx[order]
    group_first_idx = np.concatenate([[0], np.nonzero(np.diff(sorted_row_idx))[0] + 1])

    best_order_positions = order[group_first_idx]
    best_rows = row_idx[best_order_positions]
    best_starts = start_cols[best_order_positions]
    best_ends = end_cols[best_order_positions]

    needs_row_indices = np.where(needs_stretch_search)[0]
    left_cutoffs[needs_row_indices[best_rows]] = best_starts
    right_cutoffs[needs_row_indices[best_rows]] = best_ends

    return left_cutoffs, right_cutoffs

def adapter_trimming(sequence_arr, trim_sequences, adapter_trim_flag):
    """
    Determines per-read trim boundaries to remove specific adapter sequences.
    Detects adapter read-through from the end of each read.

    For each read, finds the earliest position where any provided adapter
    sequence appears as an exact substring, and trims the read at that position.
    If `adapter_trim_flag` is False, returns boundaries that trim nothing 
    (left=0, right=read_length) for every read.

    Args:
        sequence_arr (numpy.ndarray): (n_reads, read_length) array of
            per-base ASCII sequence codes.
        trim_sequences (list[tuple[str, str]]): List of (name, sequence) tuples 
            representing the adapter sequences to search for.
        adapter_trim_flag (bool): Flag to enable or disable adapter trimming.

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: (left_cutoffs, right_cutoffs),
            each of shape (n_reads,) and dtype int32, giving the left and right 
            trim boundaries per read. Left cutoffs are always 0 (3'-end trimming only).
    """
    n_reads, length = sequence_arr.shape
    if not adapter_trim_flag:
        return np.zeros(n_reads, dtype = np.int32), np.full(n_reads, length, dtype = np.int32)
    right_cutoffs = np.full(n_reads, length, dtype=np.int32)
    
    
    adapter_byte_strs = list({seq.encode('utf-8') for _, seq in trim_sequences})
    sequence_arr = sequence_arr.astype(np.uint8)

    for i in range(n_reads):
        row_bytes = sequence_arr[i].tobytes()
        best = length
        for adapter_bytes in adapter_byte_strs:
            pos = row_bytes.find(adapter_bytes)
            if pos != -1 and pos < best:
                best = pos
        right_cutoffs[i] = best

    return np.zeros(n_reads, dtype=np.int32), right_cutoffs

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
    n_reads = quality_arr.shape[0]
    cumsum = np.concatenate([np.zeros((n_reads, 1), dtype = np.int32), np.cumsum(quality_arr, axis = 1, dtype=np.int32)], axis = 1)
    row_indices = np.arange(n_reads)
    sums = cumsum[row_indices, rights] - cumsum[row_indices, lefts]
    counts = rights - lefts
    avg_qualities = np.divide(sums, counts, out = np.zeros_like(sums, dtype = np.float64), where=counts > 0)
    return avg_qualities

def kmer_complexity_scan(sequence_arr, kmer_scan, kmer, low_complex_cutoff, allow_n):
    """
    Counts k-mers per read and flags reads falling below the specified complexity
    cutoff for removal.

    Args:
        sequence_arr (numpy.ndarray): A 2D array of ASCII sequence codes
            of shape (n_reads, length).
        kmer_scan (bool): If False, bypasses the k-mer calculation entirely 
            and returns default placeholder arrays.
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
    if not kmer_scan:
        return np.zeros(n_reads, dtype=np.int32), np.full(n_reads, length, dtype=np.int32)

    if isinstance(kmer, str):
        kmer_list = [int(k.strip()) for k in kmer.split(',')]
    elif hasattr(kmer, '__iter__'):
        kmer_list = [int(k) for k in kmer]
    else:
        kmer_list = [int(kmer)]

    global_passed = np.ones(n_reads, dtype=bool)
    for k in kmer_list:
        if k > length:
            raise ValueError(f"k-mer length {k} is greater than sequence length {length}")
        mapping = np.zeros(256, dtype=np.int8)
        if allow_n:
            mapping[ord('A')] = 0
            mapping[ord('C')] = 1
            mapping[ord('G')] = 2
            mapping[ord('T')] = 3
            mapping[ord('N')] = 4
            bits_per_base = 3
            max_val = 1 << (3 * k)
        else:
            mapping[ord('A')] = 0
            mapping[ord('C')] = 1
            mapping[ord('G')] = 2
            mapping[ord('T')] = 3
            bits_per_base = 2
            max_val = 1 << (2 * k)
        int_matrix = mapping[sequence_arr]
        max_kmers = length - k + 1
        kmer_ints = np.zeros((n_reads, max_kmers), dtype=np.int64)
        for i in range(k):
            kmer_ints = (kmer_ints << bits_per_base) | int_matrix[:, i:i+max_kmers]

        if max_val <= 4096:
            per_read_counts = np.zeros((n_reads, max_val), dtype=np.int64)
            for i in range(n_reads):
                per_read_counts[i] = np.bincount(kmer_ints[i], minlength=max_val)
                unique_count = np.count_nonzero(per_read_counts[i])
                if (unique_count / max_kmers) < (low_complex_cutoff/100):
                    global_passed[i] = False
        else:
            for i in range(n_reads):
                unique_count = np.unique(kmer_ints[i]).size
                if (unique_count / max_kmers) < (low_complex_cutoff/100):
                    global_passed[i] = False

    second_array = np.where(global_passed, length, 0).astype(np.int32)
    return np.zeros(n_reads, dtype=np.int32), second_array

##### Unpaired reads workflow functions #####
def process_unpaired_chunk(chunk, phred_offset, minimum_length, maximum_length, minimum_average_qual, read_length, gzip_output, gzip_level, parameters):
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
        chunk (list[tuple]): A list of (header, sequence, plus, quality) record tuples.
        phred_offset (int): Phred encoding offset (33 or 64).
        minimum_length (int): Minimum acceptable read length after trimming.
        maximum_length (int): Maximum acceptable read length after trimming.
        minimum_average_qual (float): Minimum acceptable mean quality after trimming.
        read_length (int): Expected original read length (for validation).
        gzip_output (bool): If True, compresses output records with gzip.
        gzip_level (int): Gzip compression level (1-9).
        parameters (dict): Dictionary of configuration parameters.

    Returns:
        tuple[bytes | str, int, int]: A tuple containing:
            - Formatted and optionally gzipped FASTQ records for surviving reads.
            - Count of kept reads.
            - Count of rejected reads.
    """
    valid_headers = []
    valid_sequences = []
    valid_pluses = []
    valid_qualities = []
    for r in chunk:
        if validate_fastq(*r, min_raw_read_length = parameters["min_raw_read_length"], nucleotide_regex = parameters["nucleotide_regex"], read_length = read_length):
            valid_headers.append(r[0])
            valid_sequences.append(r[1])
            valid_pluses.append(r[2])
            valid_qualities.append(r[3])
    rejected = len(chunk) - len(valid_headers)
    if not valid_headers:
        return [], 0, rejected
    if parameters["mgi_convert_flag"]:
        valid_pluses = [f"{plus}_OriginalHeader:{header}" for plus, header in zip(valid_pluses, valid_headers)]
        valid_headers = [header_mgi_to_illumina(header, parameters["mgi_bc5"], parameters["mgi_bc7"], parameters["mgi_instrument"], parameters["mgi_run"]) for header in valid_headers]
    quality_arr = qual_to_bin(quality_list = valid_qualities, phred_offset = phred_offset)
    sequence_arr = seq_to_bin(sequence_list = valid_sequences)
    left_list = []
    right_list = []
    for left, right in [
        trim_ends_quality(quality_arr, min_quality_both = parameters["min_quality_both"], endqual_min_start = parameters["endqual_min_start"], endqual_min_end = parameters["endqual_min_end"], endqual_filter_flag = parameters["endqual_filter_flag"]),
        sliding_window_quality(quality_arr, slider_quality = parameters["slider_quality"], slider_window = parameters["slider_window"], slider_step = parameters["slider_step"], slider_filter_flag = parameters["slider_filter_flag"]),
        homopolymer_nucleotide_trimming(sequence_arr, poly_length_both = parameters["poly_length_both"], poly_length_start = parameters["poly_length_start"], poly_length_end = parameters["poly_length_end"], poly_bases_both = parameters["poly_bases_both"], poly_bases_start = parameters["poly_bases_start"], poly_bases_end = parameters["poly_bases_end"], poly_filter_flag = parameters["poly_filter_flag"]),
        adapter_trimming(sequence_arr, parameters["adapter_sequences"], adapter_trim_flag=parameters["adapter_trim_flag"]),
        cut_set_ends(sequence_arr, cut_both = parameters["cut_both"], cut_start = parameters["cut_start"], cut_end = parameters["cut_end"], cut_flag = parameters["cut_flag"]),
        n_end_trimming(sequence_arr, n_trimming_flag = parameters["n_trimming_flag"]),
        kmer_complexity_scan(sequence_arr, kmer_scan = parameters["kmer_filter_flag"], kmer = parameters["kmer_size"], low_complex_cutoff = parameters["kmer_cutoff"], allow_n = parameters["allow_n_kmer"])
        ]:
        left_list.append(left)
        right_list.append(right)
    lefts = np.maximum.reduce(left_list)
    rights = np.minimum.reduce(right_list)
    if maximum_length == minimum_length == sequence_arr.shape[1]:
        length_mask = np.ones(len(lefts), dtype=bool)
    else:
        lengths_out = rights - lefts
        length_mask = (lengths_out <= maximum_length) & (lengths_out >= minimum_length)
    if minimum_average_qual > 0:
        average_quals = average_quality_batch(quality_arr, lefts, rights)
        qual_mask = average_quals >= minimum_average_qual
    else:
        qual_mask = np.ones(len(lefts), dtype=bool)
    keep_mask = length_mask & qual_mask
    results = []
    for header, sequence, plus_line, quality, left, right, keep in zip(
            valid_headers, valid_sequences, valid_pluses, valid_qualities, lefts, rights, keep_mask):
        if not keep:
            rejected += 1
            continue
        results.append(f"{header}\n{sequence[left:right]}\n{plus_line}\n{quality[left:right]}\n")
    len_results = len(results)
    results = "".join(results)
    if gzip_output:
        results = gzip.compress(results.encode('utf-8'), compresslevel=gzip_level)
    else:
        results = results.encode('utf-8')
    return results, len_results, rejected

def generate_unpaired_tasks(filepaths, chunk_size, parameters):
    """
    Lazily yields individual chunks alongside their filepaths and precomputed parameters,
    allowing a single global pool to process chunks from multiple files concurrently.

    Args:
        filepaths (list[str]): List of paths to unpaired FASTQ files.
        chunk_size (int): Number of reads per chunk.
        parameters (dict): Dictionary of configuration parameters.

    Yields:
        dict: A task dictionary containing the task type, filepath, chunk data,
            and precomputed metadata.
    """
    for filepath in filepaths:
        phred_offset = detect_phred_offset(
            filepath=filepath,
            reads_for_phred_offset=parameters["reads_for_phred_offset"],
            phred_offset=parameters["phred_offset"]
        )
        read_length = query_read_length(filepath)
        minimum_len = parameters["minimum_length"]
        maximum_len = parameters["maximum_length"]
        if minimum_len is None:
            minimum_len = int(read_length // 3)
        if maximum_len is None:
            maximum_len = read_length
        reads_iter = lazy_fastq(filepath)
        while True:
            chunk = list(itertools.islice(reads_iter, chunk_size))
            if not chunk:
                break
            yield {
                "type": "unpaired",
                "filepath": filepath,
                "chunk": chunk,
                "phred_offset": phred_offset,
                "read_length": read_length,
                "minimum_length": minimum_len,
                "maximum_length": maximum_len,
                "parameters": parameters
            }

def process_unpaired_task_flat(task, parameters):
    """
    Worker wrapper for flat task queue execution on unpaired files.

    Args:
        task (dict): A task dictionary containing chunk data, filepaths, and metadata.
        parameters (dict): Dictionary of configuration parameters.

    Returns:
        tuple: A tuple containing (type, filepath, chunk_results, kept, rejected).
    """
    chunk_results, kept, rejected = process_unpaired_chunk(
        chunk=task["chunk"],
        phred_offset=task["phred_offset"],
        minimum_length=task["minimum_length"],
        maximum_length=task["maximum_length"],
        minimum_average_qual=parameters["minimum_average_qual"],
        read_length=task["read_length"],
        gzip_output = parameters["gzip_output"], 
        gzip_level = parameters["gzip_level"], 
        parameters=parameters
    )
    return task["type"], task["filepath"], chunk_results, kept, rejected

##### Paired reads workflow funtions #####
def process_paired_task_flat(task, parameters):
    """
    Worker wrapper for flat task queue execution on paired-end files.

    Args:
        task (dict): A task dictionary containing chunks for R1 and R2, filepaths, and metadata.
        parameters (dict): Dictionary of configuration parameters.

    Returns:
        tuple: A tuple containing (type, file1, file2, paired_out_1, paired_out_2, 
            singles_out, num_paired, num_singles, rejected).
    """
    file1 = task["file1"]
    file2 = task["file2"]
    paired_out_1, paired_out_2, singles_out, num_paired, num_singles, rejected = process_paired_chunk(
        chunks = (task["chunk1"],task["chunk2"] ),
        phred_offset=task["phred_offset"],
        minimum_length=task["minimum_length"],
        maximum_length=task["maximum_length"],
        read_length=task["read_length"],
        gzip_output = task["gzip_output"], 
        gzip_level = task["gzip_level"], 
        parameters=parameters
    )
    return task["type"], file1, file2, paired_out_1, paired_out_2, singles_out, num_paired, num_singles, rejected

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
        phred_offset_1 = detect_phred_offset(filepath = file1, reads_for_phred_offset = parameters["reads_for_phred_offset"], phred_offset = parameters["phred_offset"])
        read_length_1 = query_read_length(file1)
        phred_offset_2 = detect_phred_offset(filepath = file2, reads_for_phred_offset = parameters["reads_for_phred_offset"], phred_offset = parameters["phred_offset"])
        read_length_2 = query_read_length(file2)
        if read_length_1 != read_length_2 or phred_offset_1 != phred_offset_2:
            print(
                f"Paired files must match in read length, read count, and Phred offset. "
                f"offsets ({phred_offset_1}, {phred_offset_2})."
                f"read lengths ({read_length_1}, {read_length_2})."
            )
            continue
        minimum_length = parameters["minimum_length"]
        maximum_length = parameters["maximum_length"]
        if minimum_length is None:
            minimum_length = read_length_1 // 3
        if maximum_length is None:
            maximum_length = read_length_1
        reads_iter_1 = lazy_fastq(file1)
        reads_iter_2 = lazy_fastq(file2)
        while True:
            chunk1 = list(itertools.islice(reads_iter_1, chunk_size))
            chunk2 = list(itertools.islice(reads_iter_2, chunk_size))
            if not chunk1 and not chunk2:
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
                "phred_offset": phred_offset_1,
                "read_length": read_length_1,
                "minimum_length": minimum_length,
                "maximum_length": maximum_length,
                "gzip_output": parameters["gzip_output"], 
                "gzip_level": parameters["gzip_level"], 
                "parameters": parameters
            }

def trim_reads(records, phred_offset, minimum_length, maximum_length, read_length, minimum_average_qual, parameters):
    """
    Validates, quality-trims, and length/quality-filters a batch of FASTQ
    reads, keyed by their base (mate-independent) read ID.

    Same trimming logic as `process_unpaired_chunk`, but returns a dict
    keyed by base read ID rather than a flat list of formatted strings —
    this allows the paired workflow to later match up surviving R1/R2 mates
    by ID.

    Args:
        records (list[tuple]): (header, sequence, plus, quality) tuples.
        phred_offset (int): Phred encoding offset (33 or 64).
        minimum_length (int): Minimum acceptable read length after trimming.
        maximum_length (int): Maximum acceptable read length after trimming.
        read_length (int): Expected original read length (for validation).
        minimum_average_qual (float): Minimum acceptable mean quality after trimming.
        parameters (dict): Dictionary of configuration parameters.

    Returns:
        tuple[dict[str, str], int]: A tuple containing:
            - A dictionary mapping each surviving read's base ID to its formatted FASTQ record string.
            - Total count of rejected reads.
    """
    valid_headers = []
    valid_sequences = []
    valid_pluses = []
    valid_qualities = []
    rejected = 0
    for r in records:
        if validate_fastq(*r,min_raw_read_length = parameters["min_raw_read_length"],nucleotide_regex = parameters["nucleotide_regex"], read_length = read_length):
            valid_headers.append(r[0])
            valid_sequences.append(r[1])
            valid_pluses.append(r[2])
            valid_qualities.append(r[3])
        else:
            rejected += 1
    if not valid_headers:
        return {}, len(records) // 4
    if parameters["mgi_convert_flag"]:
        valid_pluses = [f"{plus}_OriginalHeader:{header}" for plus, header in zip(valid_pluses, valid_headers)]
        valid_headers = [header_mgi_to_illumina(header, parameters["mgi_bc5"], parameters["mgi_bc7"], parameters["mgi_instrument"], parameters["mgi_run"]) for header in valid_headers]
    quality_arr = qual_to_bin(quality_list = valid_qualities, phred_offset = phred_offset)
    sequence_arr = seq_to_bin(sequence_list = valid_sequences)
    left_list, right_list = [], []
    for left, right in [
        trim_ends_quality(quality_arr, min_quality_both = parameters["min_quality_both"], endqual_min_start = parameters["endqual_min_start"], endqual_min_end = parameters["endqual_min_end"], endqual_filter_flag = parameters["endqual_filter_flag"]),
        sliding_window_quality(quality_arr, slider_quality = parameters["slider_quality"], slider_window = parameters["slider_window"], slider_step = parameters["slider_step"], slider_filter_flag = parameters["slider_filter_flag"]),
        homopolymer_nucleotide_trimming(sequence_arr, poly_length_both = parameters["poly_length_both"], poly_length_start = parameters["poly_length_start"], poly_length_end = parameters["poly_length_end"], poly_bases_both = parameters["poly_bases_both"], poly_bases_start = parameters["poly_bases_start"], poly_bases_end = parameters["poly_bases_end"], poly_filter_flag = parameters["poly_filter_flag"]),
        adapter_trimming(sequence_arr, parameters["adapter_sequences"], adapter_trim_flag=parameters["adapter_trim_flag"]),
        cut_set_ends(sequence_arr, cut_both = parameters["cut_both"], cut_start = parameters["cut_start"], cut_end = parameters["cut_end"], cut_flag = parameters["cut_flag"]),
        n_end_trimming(sequence_arr, n_trimming_flag = parameters["n_trimming_flag"]),
        kmer_complexity_scan(sequence_arr, kmer_scan = parameters["kmer_filter_flag"], kmer = parameters["kmer_size"], low_complex_cutoff = parameters["kmer_cutoff"], allow_n = parameters["allow_n_kmer"])
    ]:
        left_list.append(left)
        right_list.append(right)
    lefts = np.maximum.reduce(left_list)
    rights = np.minimum.reduce(right_list)
    if maximum_length == minimum_length == sequence_arr.shape[1]:
        length_mask = np.ones(len(lefts), dtype=bool)
    else:
        lengths_out = rights - lefts
        length_mask = (lengths_out <= maximum_length) & (lengths_out >= minimum_length)
    if minimum_average_qual > 0:
        avg_quals = average_quality_batch(quality_arr, lefts, rights)
        qual_mask = avg_quals >= minimum_average_qual
    else:
        qual_mask = np.ones(len(lefts), dtype=bool)
    
    keep_mask = length_mask & qual_mask
    survivors = {}
    for i, keep in enumerate(keep_mask):
        if not keep:
            rejected += 1
            continue
        left, right = int(lefts[i]), int(rights[i])
        seq_out = valid_sequences[i][left:right]
        qual_out = valid_qualities[i][left:right]
        base_id, _ = read_info_from_header(valid_headers[i])
        survivors[base_id] = f"{valid_headers[i]}\n{seq_out}\n{valid_pluses[i]}\n{qual_out}\n"

    return survivors, rejected

def process_paired_chunk(chunks, phred_offset, read_length, minimum_length, maximum_length, gzip_output, gzip_level, parameters):
    """
    Trims and filters one paired chunk of R1/R2 reads, then reconciles the
    two mates by base read ID to determine which reads survive as intact
    pairs versus as orphaned singletons.

    Args:
        chunks (tuple[list, list]): (chunk1, chunk2) — record lists for R1
            and R2 respectively, covering the same reads in the same order.
        phred_offset (int): Phred encoding offset (33 or 64).
        read_length (int): Expected original read length (for validation).
        minimum_length (int): Minimum acceptable read length after trimming.
        maximum_length (int): Maximum acceptable read length after trimming.
        gzip_output (bool): If True, compresses output records with gzip.
        gzip_level (int): Gzip compression level (1-9).
        parameters (dict): Dictionary of configuration parameters.

    Returns:
        tuple[bytes, bytes, bytes, int, int, int]: A tuple containing:
            - Surviving R1 records whose R2 mate also survived (optionally gzipped).
            - Surviving R2 records whose R1 mate also survived (optionally gzipped).
            - Surviving records whose mate did not survive, treated as unpaired singletons (optionally gzipped).
            - Count of surviving read pairs.
            - Count of surviving singleton reads.
            - Total count of rejected reads.
    """
    chunk1, chunk2 = chunks
    survivors_1, rejected_1 = trim_reads(chunk1, phred_offset, minimum_length, maximum_length, minimum_average_qual = parameters["minimum_average_qual"], read_length = read_length, parameters = parameters)
    survivors_2, rejected_2 = trim_reads(chunk2, phred_offset, minimum_length, maximum_length, minimum_average_qual = parameters["minimum_average_qual"], read_length = read_length, parameters = parameters)
    paired_out_1 = []
    paired_out_2 = []
    paired_ids = set()
    for bid, record in survivors_1.items():
        if bid in survivors_2:
            paired_out_1.append(record)
            paired_out_2.append(survivors_2[bid])
            paired_ids.add(bid)
    singles_out = []
    singles_out.extend(v for bid, v in survivors_1.items() if bid not in paired_ids)
    singles_out.extend(v for bid, v in survivors_2.items() if bid not in paired_ids)
    num_singles = len(singles_out)
    num_paired = len(paired_out_1)
    paired_out_1 = "".join(paired_out_1)
    paired_out_2 = "".join(paired_out_2)
    singles_out = "".join(singles_out)

    if gzip_output:
        paired_out_1 = gzip.compress(paired_out_1.encode('utf-8'), compresslevel=gzip_level)
        paired_out_2 = gzip.compress(paired_out_2.encode('utf-8'), compresslevel=gzip_level)
        singles_out = gzip.compress(singles_out.encode('utf-8'), compresslevel=gzip_level)
    else:
        paired_out_1 = paired_out_1.encode('utf-8')
        paired_out_2 = paired_out_2.encode('utf-8')
        singles_out = singles_out.encode('utf-8')
    return paired_out_1, paired_out_2, singles_out, num_paired, num_singles, rejected_1 + rejected_2

##### Input handler functions #####
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
    parameters = task.get("parameters")

    if task_type == "paired":
        return process_paired_task_flat(task, parameters=parameters)
    elif task_type == "unpaired":
        return process_unpaired_task_flat(task, parameters=parameters)
    else:
        raise ValueError(f"Unknown or missing task type: {task_type}")

def input_handler(unspecified_files, unpaired_files, paired_files, output_dir, threads, chunk_size, parameters):
    """
    Top-level orchestrator that separates input files into paired and
    unpaired groups, sets up file writers, and runs the multiprocessing pool
    workflow across all inputs.

    Args:
        unspecified_files (list[str]): Files with unknown pairing status to auto-detect.
        unpaired_files (list[str]): Explicitly provided unpaired input FASTQ files.
        paired_files (list[tuple[str, str]]): Explicitly provided paired FASTQ file pairs.
        output_dir (str): Directory to write output files to.
        threads (int | None): Number of worker threads/processes for the pool.
        chunk_size (int): Number of reads per processing chunk.
        parameters (dict): Dictionary of configuration parameters.

    Returns:
        dict: A dictionary containing statistics and counts for kept and rejected reads per file.
    """
    auto_paired, auto_unpaired = find_paired_files(unspecified_files)
    unpaired = auto_unpaired + (unpaired_files or [])
    paired = auto_paired + (paired_files or [])
    
    total_reads = 0
    for file in unpaired:
        total_reads += count_reads_estimated(file)
    for f1, f2 in paired:
        total_reads += count_reads_estimated(f1) + count_reads_estimated(f2)
        
    class _NullTracker:
        def update(self, n): pass
        def close(self): pass
    tracker = ProgressTracker(total_reads) if parameters["show_progress"] else _NullTracker()
    
    file_stats = {}
    file_writing_handles = {}
    for file in unpaired:
        file_stats[file] = {"kept": 0, "rejected": 0}
        file_writing_handles[file] = open_fastq_writer(file, output_dir, gzip_output=parameters["gzip_output"])
    for pair in paired:
        file1, file2 = pair
        common_prefix = common_name_parts([os.path.basename(file1), os.path.basename(file2)])
        file_stats[common_prefix] = {"kept_pairs": 0, "kept_singletons": 0, "rejected": 0}
        for suffix in ["_R1_paired", "_R2_paired", "_unpaired"]:
                    key = f"{common_prefix}{suffix}"
                    # Assign the returned handle directly to your dictionary
                    file_writing_handles[key] = open_fastq_writer(
                        key, 
                        output_dir, 
                        gzip_output=parameters["gzip_output"]
                    )
    def unified_chunk_streamer():
        yield from generate_unpaired_tasks(filepaths=unpaired, chunk_size=chunk_size, parameters=parameters)
        yield from generate_paired_tasks(files=paired, chunk_size=chunk_size, parameters=parameters)
    try:
        with Pool(threads) as pool:
            for result in pool.imap_unordered(unified_worker, unified_chunk_streamer(), chunksize=1):
                if result[0] == "unpaired":  # Unpaired result: (type, filepath, chunk_results, kept, rejected)
                    _, filepath, chunk_results, kept, rejected = result
                    if chunk_results:
                        file_writing_handles[filepath].write(chunk_results)
                    file_stats[filepath]["kept"] += kept
                    file_stats[filepath]["rejected"] += rejected
                    tracker.update(kept + rejected)   

                elif result[0] == "paired":  # Paired result: (type, file1, file2, p1, p2, singles, rejected)
                    _, file1, file2, paired_out_1, paired_out_2, singles_out, num_paired, num_singles, rejected = result
                    common_prefix = common_name_parts([os.path.basename(file1), os.path.basename(file2)])
                    writes = [
                        (f"{common_prefix}_R1_paired", paired_out_1),
                        (f"{common_prefix}_R2_paired", paired_out_2),
                        (f"{common_prefix}_unpaired", singles_out)
                    ]
                    for handle_key, records in writes:
                        if records:
                            file_writing_handles[handle_key].write(records)                     
                    file_stats[common_prefix]["kept_pairs"] += num_paired
                    file_stats[common_prefix]["kept_singletons"] += num_singles
                    file_stats[common_prefix]["rejected"] += rejected
                    tracker.update(num_paired * 2 + num_singles + rejected)   
    finally:
        tracker.close()
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
        # Match: description line, blank line, then an indented option line
        # (option lines start with 2+ spaces followed by "-")
        help_text = re.sub(r"\n\n(?=  -)", "\n", help_text)
        return help_text
    
def print_full_auto_help(parser):
    """
    Prints the effective settings --full-auto/-GO applies, grouped by the
    parser's argument groups, using each argument's registered default
    (or its full-auto override, where one exists).
    """
    print(
        "\n"
        "When --full-auto/-GO is specified, only input parameters --input-files/-i, --input-paired/-ip, and --input-unpaired/-iu are respected."
        "\nIf none of these are provided, ReadZor will auto-detect FASTQ files in the current working directory."
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
    Parses ReadZor's command-line arguments into a resolved parameters dict.
    
    Supports three main input modes: fully automatic operation (--full-auto,
    which autodetects files and ignores all other options), a flat list of
    FASTQ files (--files), or an explicit R1/R2 pair (--paired). The latter
    two are mutually exclusive.
    
    All trimming, quality, and runtime parameters default to None if not
    specified, signaling to downstream code that an adaptive/platform-
    dependent default should be resolved later (e.g. based on read length
    or Slurm detection) rather than being hardcoded here.
    
    Returns:
        dict: A parameters dictionary with the following keys:
            - full_auto (bool): Whether full-auto mode was requested.
            - paired (bool): Whether --paired was used (vs --files).
            - files (list[str] or None): Resolved input file(s), or None
              if full-auto (to be autodetected downstream).
            - minimum_length, maximum_length (int or None): Read length bounds.
            - min_quality_both (int or None): Phred score for end trimming
              (applied to both ends unless overridden).
            - endqual_min_start, endqual_min_end (int or None): Phred score
              thresholds for trimming the 5' and 3' ends specifically;
              override min_quality_trim for their respective end if given.
            - minimum_average_qual (int or None): Minimum average read quality.
            - cut_start, cut_end (int or None): Fixed number of bases to trim
              from the start and end of the read specifically.
            - cut_both (int or None): Fixed number of bases to trim from both
              ends; mutually exclusive with cut_start/cut_end.
            - slider_window (int): Sliding window size for quality trimming.
              Defaults to 5.
            - slider_step (int): Step size between successive sliding windows.
              Defaults to 1.
            - poly_bases (str): Base to check for a trailing homopolymer run
              (e.g. poly-G trimming). Defaults to "G".
            - poly_length_start (int): Minimum homopolymer run length of poly_bases
              required to trigger trimming. Defaults to 10.
            - min_raw_read_length (int): Minimum length a raw (untrimmed) read
              must have to be considered valid. Defaults to 0.
            - reads_for_phred_offset (int): Number of reads to sample when
              auto-detecting the Phred quality encoding offset. Defaults to 1000.
            - gzip_output (bool): Whether to gzip-compress the output file(s).
            - gzip_level (int): Gzip compression level (1-9). Defaults to 6.
            - threads (int or None): Requested CPU count.
            - chunk_size (int or None): Requested chunk size.
            - output_dir (str): Resolved output directory (defaults to the
              current working directory if not specified).
 
    Side Effects:
        If --version is passed, prints the version string and exits the
        program immediately (via `exit()`). If neither --files nor --paired
        nor --full-auto is given, calls `parser.error(...)`, which prints a
        usage message to stderr and exits with a non-zero status.
    """
    parser = argparse.ArgumentParser(
        prog="readzor",
        description="ReadZor: a modular FASTQ quality trimming pipeline.\n\n"
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
        help="[FLAG] Show ReadZor version and exit."
    )
    general_group.add_argument(
        "--full-auto", "-GO", action = "store_true", default = False,
        help="[FLAG] Run ReadZor in fully automatic mode. Combine with --help/-h for more information on fully automatic mode."
    )
    general_group.add_argument(
    "--progress", action="store_true", default=False,
    help="[FLAG] Show a live progress bar and estimated time remaining during processing, based on estimated read counts. Default: off."
    )
 
    input_group = parser.add_argument_group(
        "Input options",
        "Specify input FASTQ files using any combination of --input-files, --input-paired, and --input-unpaired. "
        "Lists with any combination of regular (fastq/fq), and gzipped (fastq.gz/fq.gz) files accepted."
    )
    input_group.add_argument(
        "--input-files", "-i", nargs='+', default = None, metavar = "",
        help="FASTQ files of unspecified pairing. Paired and unpaired files will be auto-detected."
    )
    input_group.add_argument(
        "--input-paired", "-ip", nargs='+', default=None, metavar = "",
        help="Paired-end FASTQ files, given as one or more R1/R2 pairs, e.g. --input-paired sample1_R1 sample1_R2 sample2_R1 sample2_R2"
    )
    input_group.add_argument(
        "--input-unpaired", "-iu", nargs='+', default = None, metavar = "",
        help="Unpaired FASTQ files."
    )
 
    output_group = parser.add_argument_group("Output options")
    output_group.add_argument(
        "--output", "-o", type=str, default = None, metavar = "",
        help="Path to directory in which the timestamped results folder will be created. Default: current working directory."
    )
    output_group.add_argument(
        "--gzip", action="store_true", default = False,
        help="[FLAG] Compress filtered FASTQ files in gzip format. Default: off."
    )
    output_group.add_argument(
        "--gzip-level", type=int, default = 4, metavar = "", choices=range(1, 10),
        help="Set gzip compression level. Higher compression decreases processing speed. Possible levels: 1-9. Default: 4."
    )
 
    general_quality_group = parser.add_argument_group("General output filter options")
    general_quality_group.add_argument(
        "--min-average-qual", type=int, default = 0, metavar = "", choices=range(0, 127),
        help="Minimum average quality of output read. Default: 0."
    )
    general_quality_group.add_argument(
        "--min-length", type=int, default = 0, metavar = "", 
        help="Minimum length of output read. Default: 0."
    )
    general_quality_group.add_argument(
        "--max-length", type = int, default = None, metavar = "", 
        help="Maximum length of output read. Default: off."
    )
    general_quality_group.add_argument(
        "--nucl-filter", action="store_true", default=False,
        help="[FLAG] Reject raw reads containing N bases anywhere in read. Default: off."
    )
 
    trim_ends_group = parser.add_argument_group("Set-length end trimming",
                                                "Trim a set number of bases of the ends of each read, independent of sequence or quality."
                                                )
    trim_ends_group.add_argument(
        "--cut-flag", "-clf",action="store_true", default = False,
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
        "--endqual-filter-flag", "-eff", action="store_true", default = False,
        help="[FLAG] Turn on quality-dependent end trimming. Default: off."
    )
    quality_ends_group.add_argument(
        "--endqual-min-start", "-ems", type = int, default = 25, metavar="", choices=range(0, 127),
        help="Specific phred score threshold for the start of the read. Default: 25."
    )
    quality_ends_group.add_argument(
        "--endqual-min-end", "-eme", type=int, default = 25, metavar="",choices=range(0, 127),
        help="Specific phred score threshold for the end of the read.  Default: 25."
    )
    quality_ends_group.add_argument(
        "--endqual-min-both", "-emb", type=int, default = 25, metavar="", choices=range(0, 127),
        help="Phred score threshold for the quality trimming of read ends. Overwritten by --endqual-min-start and --endqual-min-end. Default: 25."
    )
    
    n_ends_group = parser.add_argument_group("N nucleotide end-trimming",
                                             "Trim the ends of each read for N bases. Redundant when --nucl-filter is set.")
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
        "--poly-length-start", "-pls", type = int, default = 0, metavar = "",
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
                                                 "Trim reads for Illumina adapter sequences. Standard sequences included are TruSeq3 universal and index adapters, and Nextera adapters. Only perfectly matching sequences are trimmed. Indepedent of quality.")
    adapter_trimming.add_argument(
        "--adapter-trim-flag", "-af", action="store_true", default = False,
        help='[FLAG] Turn on adapter trimming module. Default: off.'
    )
    adapter_trimming.add_argument(
        "--adapter-fasta", "-ad", type = str, default = None, metavar="",
        help="Fasta file with additional adapter sequences to trim for."
    )
 
    low_complexity_group = parser.add_argument_group("Low complexity filtering",
                                                     "Detect complexity of reads using kmer-based nucleotide frequencies. Low complex reads discarded entirely.")
    low_complexity_group.add_argument(
        "--kmer-filter-flag", action="store_true", default=False,
        help="[FLAG] Turn on the kmer-based complexity filtering module. Default: off."
    )
    low_complexity_group.add_argument(
        "--kmer-size", type = int, default = 4, metavar="",
        help="Kmer length for kmer-based complexity filtering. Comma-separated values are checked independently. Default: 4."
    )
    low_complexity_group.add_argument(
        "--kmer-cutoff", type = int, default = 50, metavar="",
        help="Minimum percentage of unique k-mers (relative to the maximum possible for the read) required to pass the complexity filter. Higher values are stricter. Default: 50."
    )
    mgi_convert_group = parser.add_argument_group("MGI header conversion",
                                                  "Convert read header from MGI (BGI) format to Illumina format. Original header will be stored in the placeholder line. Conversion is necessary for downstream analysis with tools such as samtools")
    mgi_convert_group.add_argument(
        "--mgi-convert-flag", action="store_true", default=False,
        help="[FLAG] Turn on the MGI-to-Illumina header conversion module. Default: off."
    )
    mgi_convert_group.add_argument(
        "--mgi-bc5", type = str, default = "PLACEHOLDERi5", metavar="",
        help="Input a i5 barcode for Illumina header conversion. Default: 'PLACEHOLDERi5'."
    )
    mgi_convert_group.add_argument(
        "--mgi-bc7", type = str, default = "PLACEHOLDERi7", metavar="",
        help="Input a i7 barcode for Illumina header conversion. Default: 'PLACEHOLDERi7'."
    )
    mgi_convert_group.add_argument(
        "--mgi-instrument", type = str, default = "PLACEHOLDERinstrument", metavar="",
        help="Instrument name for Illumina header conversion. Default: 'PLACEHOLDERinstrument'."
    )
    mgi_convert_group.add_argument(
        "--mgi-run", type = str, default = "PLACEHOLDERrun", metavar="",
        help="Run ID for Illumina header conversion. Default: 'PLACEHOLDERrun'."
    )    
 
    advanced_group = parser.add_argument_group("Advanced options",
                                               "Further options that can be specified to alter the behaviour of ReadZor.")
    advanced_group.add_argument(
        "--threads", "-t", type = int, default = None, metavar="",
        help="Number of threads to use (defaults: platform-dependent through auto-detection: detection of assigned CPUs on Slurm-managed systems, all-1 otherwise. Fallback: 1."
    )
    advanced_group.add_argument(
        "--min-raw-read-length", type = int, default = 0, metavar="",
        help="Minimum length a raw (untrimmed) read must have to be considered valid. Default: 0."
    )
    advanced_group.add_argument(
        "--reads-for-phred-offset", type = int, default = 500, metavar="",
        help="Number of reads to sample per file for detection of Phred quality encoding offset. Default: 500."
    )
    advanced_group.add_argument(
        "--chunk-size", type=int, default = 1000, metavar="",
        help="Number of reads per chunk sent to each worker (default: platform-dependent: 20.000 for Slurm-managed systems, 1000 otherwise)."
    )
    advanced_group.add_argument(
        "--phred-offset", type = int, default = None, metavar="",
        help="Define phred offset for all FASTQ files. Default: off (auto-detection per file)."
    )
    
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
 
    args = parser.parse_args()
 
    if args.help and args.full_auto:
        print_full_auto_help(parser)
        sys.exit()
    if args.help:
        parser.print_help()
        sys.exit()
 
    if args.version:
        print(f"ReadZor version: {VERSION}")
        sys.exit()
 
    if args.full_auto:
        if not (args.input_files or args.input_paired or args.input_unpaired):
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
 
        # --- Full-auto overrides ---
        print("[WARNING] --full-auto/-GO specified; ignoring all other input parameters (except input file parameters).")
 
        for action in parser._actions:
            dest = action.dest
            if dest == "help" or dest in FULL_AUTO_PRESERVED_DESTS:
                continue
            reset_value = FULL_AUTO_OVERRIDES.get(dest, action.default)
            setattr(args, dest, reset_value)
    elif not (args.input_files or args.input_paired or args.input_unpaired):
        parser.error(
            "You must specify any combination of --input-files, --input-paired, and/or --input-unpaired (unless using --full-auto)."
        )
 
 
    # --- Store parameters ---
    parameters = {}
    parameters["full_auto"] = args.full_auto
    parameters["unspecified_files"] = args.input_files
    parameters["unpaired_files"] = args.input_unpaired
    parameters["paired_files"] = group_paired_input_into_pairs(files = args.input_paired, parser = parser)
    parameters["minimum_length"] = args.min_length
    parameters["maximum_length"] = args.max_length
    parameters["min_quality_both"] = args.endqual_min_both
    parameters["endqual_min_start"] = args.endqual_min_start
    parameters["endqual_min_end"] = args.endqual_min_end
    parameters["minimum_average_qual"] = args.min_average_qual
    parameters["cut_start"] = args.cut_start
    parameters["cut_end"] = args.cut_end
    parameters["cut_both"] = args.cut_both
    parameters["n_trimming_flag"] = args.n_trimming_flag
    parameters["slider_window"] = args.slider_window
    parameters["slider_step"] = args.slider_step
    parameters["slider_quality"] = args.slider_quality
    parameters["gzip_output"] = args.gzip
    parameters["gzip_level"] = args.gzip_level
    parameters["min_raw_read_length"] = args.min_raw_read_length
    parameters["reads_for_phred_offset"] = args.reads_for_phred_offset
    parameters["adapter_trim_flag"] = args.adapter_trim_flag
    parameters["adapter_fasta"] = args.adapter_fasta
    parameters["nucl_filter"] = args.nucl_filter
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
    parameters["allow_n_kmer"] = args.nucl_filter
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
 
    if parameters["nucl_filter"]:
        parameters["n_trimming_flag"] = False
 
    if parameters["adapter_trim_flag"] and parameters.get("adapter_fasta"):
        parameters["adapter_sequences"] = DEFAULT_ADAPTERS + load_adapters_from_fasta(parameters["adapter_fasta"])
    else:
        parameters["adapter_sequences"] = DEFAULT_ADAPTERS
 
    parameters["threads"] = worker_determination(parameters["threads"])    
    parameters["chunk_size"] = chunk_size_setter(parameters["chunk_size"])
 
    if parameters["nucl_filter"]:
        parameters["nucleotide_regex"] = STRICT_NUCLEOTIDE_REGEX
    else:
        parameters["nucleotide_regex"] = LENIENT_NUCLEOTIDE_REGEX
        
    return parameters
    
##### Wrap up functions #####
def write_summary_and_statistics(summary_results, parameters, used_command, output_dir):
    """
    Write summary statistics and parameters to output text files.

    The function creates two files in the specified output directory:
    ``results_summary.txt`` containing per-file summary counts and
    ``parameters.txt`` containing the parameter key-value pairs used
    for the analysis.
    
    Args:
        summary_results (dict): Dictionary mapping file names to dictionaries
            of summary statistics. Entries may contain either ``kept`` and
            ``rejected`` counts or ``kept_pairs``, ``kept_singletons``, and
            ``rejected`` counts.
        used_command (str): The exact shell-quoted command line used to invoke
            this run.
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
        rejected = counts["rejected"]
        if "kept" in counts:
            unpaired_data.append((filename, counts["kept"], rejected))
        else:
            paired_data.append((filename, counts["kept_pairs"], counts["kept_singletons"], rejected))
    summary_path = os.path.join(output_dir, "results_summary.txt")
    with open(summary_path, "w", encoding="utf-8") as f:
        if paired_data:
            f.write("[Paired Reads]\n")
            f.write("Pair with common prefix\tKept_Pairs\tKept_Singletons\tRejected\n")
            for item in paired_data:
                f.write(f"{item[0]}\t{item[1]}\t{item[2]}\t{item[3]}\n")
            f.write("\n")
        if unpaired_data:
            f.write("[Unpaired Reads]\n")
            f.write("Filename\tKept\tRejected\n")
            for item in unpaired_data:
                f.write(f"{item[0]}\t{item[1]}\t{item[2]}\n")
            
    params_path = os.path.join(output_dir, "parameters.txt")
    with open(params_path, "w", encoding="utf-8") as f:
        for key, value in sorted(parameters.items()):
            f.write(f"{key}\t{value}\n")
        f.write("\n")
        f.write(f"Entered command: {used_command}")

def print_final_message():
    """
    Prints ReadZor's completion message, including a citation request and
    a randomly selected sign-off phrase. Called at the end of a successful run.

    Returns:
        None
    """
    sign_off_messages = [
    "Please come again!",
    "Thanks for trimming with ReadZor!",
    "May your reads be long and your adapters be gone!",
    "Until next time, happy analyzing!",
    "See you soon!",
    "Thanks for using ReadZor!",
    "Happy analyzing!",
    "Good luck with your data!",
    "Base-ically, we're done here. See you!",
    "ReadZor, signing off!"
    ]
    print("Analysis successfully completed!")
    print("\nIf you find ReadZor useful, please consider citing:")
    print("\nAxel B. Janssen \n2026 \nReadZor: A modular and user-friendly Swiss-army knife approach to short-read sequencing preprocessing.")
    print(f"\n{random.choice(sign_off_messages)}\n")

##### Main #####
if __name__ == "__main__":
    """
    Main execution block for ReadZor.

    Initializes command-line argument parsing and parameter configuration, creates the output
    directory structure, handles input stream processing and multithreaded trimming
    via the top-level orchestrator, and exports summary metrics and final parameters.
    Finishes with a sign-off message.
    """
    parameters = parse_args()
    created_output_dir = create_folder_structure(parameters["output_dir"])
    used_command = " ".join(map(shlex.quote, [sys.executable] + sys.argv))
    summary_results = input_handler(unspecified_files = parameters["unspecified_files"], unpaired_files = parameters["unpaired_files"], paired_files = parameters["paired_files"], output_dir = created_output_dir, threads = parameters["threads"], chunk_size = parameters["chunk_size"], parameters = parameters)
    write_summary_and_statistics(summary_results, parameters, used_command, output_dir = created_output_dir)
    print_final_message()
