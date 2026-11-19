"""FASTQ input/output: reading, framing, pairing, validation and file naming."""

from collections import defaultdict
import logging
import os
import re
import shutil
import sys
import tempfile

from isal import igzip, igzip_threaded
import numpy as np

from readzor import constants
from readzor.constants import NUCL_ATCG, NUCL_ATCGN, PHRED_ALLOWED

logger = logging.getLogger(__name__)


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
    constants.STDIN_TEMP_FILES.append(tmp_path)

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
    for tmp_path in constants.STDIN_TEMP_FILES:
        try:
            os.remove(tmp_path)
        except OSError:
            pass
    constants.STDIN_TEMP_FILES.clear()


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
    if constants.GZIP_DETECTION[filepath]:
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

    Carriage returns are removed (CRLF input); other whitespace is kept.
    The blob always ends on a record boundary; record contents are checked
    later by validate_fastq().

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
            if constants.GZIP_DETECTION[filepath]:
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
        if len(headers) > 1:
            header_2 = headers[1].strip().lstrip(b'@')
            base_id_2, read_num_2 = read_info_from_header(header_2)
        else:
            base_id_2, _ = None, None
    
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
                logger.info("Paired '%s' (R1) with '%s' (R2).", os.path.basename(file_2), os.path.basename(file_1))
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
        Used whenever the header contains a space. The base ID is everything
        before the first space; the read number is taken from the start of
        the second field (None if that field does not start with 1 or 2).
      - Legacy Illumina, and MGI: "<id>/1" or "<id>/2", for headers without
        a space. The base ID is the header without the trailing "/1" or "/2".

    If a header without a space does not end in "/1" or "/2", the entire
    header is used as the base ID and the read number is set to None.

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
        m = re.match(rb'([12]):.:', parts[1])
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

    When sampling, requires GZIP_DETECTION[filepath] to be set (see
    input_handler()).

    Args:
        filepath (str): Path to the FASTQ file to inspect.
        reads_for_phred_offset (int): Maximum number of reads to sample.
        phred_offset (int | None): User-specified Phred offset, if known.

    Returns:
        int: 33 or 64, the detected or provided Phred offset.

    Raises:
        ValueError: If the FASTQ file cannot be read, contains no complete
            FASTQ records, or if the observed ASCII range fits neither
            encoding.
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
