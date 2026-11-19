"""Progress bar and read-count estimation."""

import os
import random
import shutil
import sys
import time

from isal import igzip

from readzor import constants
from readzor.fastq_io import lazy_fastq_blobs, parse_blob


def count_reads_estimated(filepath, sample_size=50, default_gzip_ratio=4, gzip_sample_bytes=50 * 1024 * 1024):
    """
    Estimate the number of reads in a FASTQ file.

    The average number of (uncompressed) bytes per read is estimated from
    the first ``sample_size`` records (or all records, if the file has
    fewer), read as one raw blob with lazy_fastq_blobs(): the size of the
    blob, newline characters included, divided by the number of records in
    it. The estimated uncompressed file size is then divided by this value.

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
            bytes-per-read estimate. Defaults to 50.
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
    constants.ESTIMATED_BYTE_PER_READ[filepath] = bytes_per_read

    file_size = os.path.getsize(filepath)
    if constants.GZIP_DETECTION[filepath]:
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
        constants.ESTIMATED_ZIP_RATIO[filepath] = ratio
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
            Defaults to 0.1.
    """
    def __init__(self, total_reads, bar_width=50, min_interval=0.1):
        self.total = max(total_reads, 1)
        self.done = 0
        self.bar_width = bar_width
        self.min_interval = min_interval
        self.last_render = time.time()
        self.start_time = time.time()

    def format_duration(self, seconds_val, concise=False):
        """
        Format time duration for display.

        Args:
            seconds_val (float | None): Duration in seconds. None, infinite
                or negative values are shown as "Calculating...".
            concise (bool): If True, return a compact form ("1h 02m 03s.");
                otherwise a spelled-out form ("1 hour, 2 minutes, 3 seconds.").

        Returns:
            str: The formatted duration, ending with a period.
        """
        if seconds_val is None or seconds_val == float('inf') or seconds_val < 0:
            return "Calculating..."
        total_seconds = int(round(seconds_val))
        hours, remainder = divmod(total_seconds, 3600)
        minutes, seconds = divmod(remainder, 60)
        def plural(value, unit):
            if value != 1:
                return f"{value} {unit}{'s'}"
            else:
                return f"{value} {unit}"
        
        if concise:
            if hours > 0:
                return f"{hours}h {minutes:02d}m {seconds:02d}s."
            if minutes > 0:
                return f"{minutes}m {seconds:02d}s."
            return f"{seconds}s."
        
        else:
            parts = []
            if hours > 0:
                parts.append(plural(hours, 'hour'))
            if minutes > 0:
                parts.append(plural(minutes, 'minute'))
            if seconds > 0 or total_seconds < 60:
                parts.append(plural(seconds, 'second'))
            return ", ".join(parts) + "."

    def update(self, n):
        """
        Add ``n`` processed reads and redraw the bar if at least
        ``min_interval`` seconds have passed since the last redraw.

        Args:
            n (int): Number of reads processed since the previous call.
        """
        self.done += n
        now = time.time()
        if now - self.last_render >= self.min_interval:
            self.render()
            self.last_render = now

    def render(self):
        """
        Draw the current progress line on stderr, overwriting the previous
        one. Once the processed count reaches the estimated total, the time
        remaining is replaced by a "finishing up" message. Does nothing if
        stderr is not a terminal.
        """
        if not sys.stderr.isatty():
            return
        term_width = shutil.get_terminal_size(fallback=(80, 24)).columns
        frac = min(self.done / self.total, 0.999)
        elapsed = time.time() - self.start_time
        rate = self.done / elapsed
        
        if self.done >= self.total:
            eta_messages = [
            "Finishing up...",
            "Almost done...",
            "Nearly there...",
            "Wrapping up..."
            ]
            eta_str = random.choice(eta_messages)
        else:
            if rate > 0:
                eta = (self.total - self.done) / rate
                if eta < 1.0:
                    eta = 1.0
            else:
                eta = float('inf')
            eta_str = self.format_duration(eta, concise=False)
            
        stats = (
            f" {frac*100:4.1f}% "
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
        the bar. The bar is redrawn by the next `update()` call that falls
        after ``min_interval``.
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
                an error). The line then says where it stopped ("Interrupted
                at 6.1%") and the bar is only filled up to that point, instead
                of claiming 100%. As the total is an estimate, the
                percentage is capped at 99.9%. Defaults to False.
        """
        elapsed = time.time() - self.start_time
        rate = self.done / elapsed
        time_str = self.format_duration(elapsed, concise=False)
        if interrupted:
            frac = min(self.done / self.total, 0.999)
            progress = f"Interrupted at {frac*100:.1f}%"
        else:
            frac = 1.0
            progress = "100%"
        stats = (
            f" {progress} ({self.done:,} reads analyzed). "
            f"Average rate: {rate:,.0f} reads/s. Total time: {time_str}"
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
