"""Logging setup: console and log-file handlers, aware of the progress bar."""

import logging
import os
import shlex
import sys

from readzor import constants
from readzor.constants import VERSION

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
        if constants.ACTIVE_PROGRESS_TRACKER is not None:
            constants.ACTIVE_PROGRESS_TRACKER.clear_line()
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
