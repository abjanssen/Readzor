"""Command-line interface and entry point."""

import argparse
import datetime
import logging
import multiprocessing as mp
import os
import random
import re
import sys
import tempfile

from readzor.constants import DEFAULT_ADAPTERS, FULL_AUTO_OVERRIDES, FULL_AUTO_PRESERVED_DESTS, VERSION
from readzor.adapter_setup import load_adapters_from_arg
from readzor.fastq_io import cleanup_stdin_temp_files, resolve_stdin_input
from readzor.log import file_only_logger, setup_logging
from readzor.runner import input_handler, log_parameters, write_summary_and_statistics

logger = logging.getLogger(__name__)

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
        parser.error(f"argument --input-paired: expected an even number of files (R1 R2 pairs), got {len(files)}.")
    pairs = []
    for i in range(0, len(files), 2):
        pairs.append(tuple(files[i:i+2]))
    return pairs

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

class CleanHelpFormatter(argparse.HelpFormatter):
    """
    Custom argparse HelpFormatter that cleans up comma spacing, starts the
    help explanations at column 50 at most (instead of argparse's 24),
    removes the blank line between a group's description and its first
    option, suppresses empty metavar placeholder artifacts (like '[ ...]'),
    and keeps manual paragraph breaks ('\\n\\n') in descriptions.
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
        """
        Return an empty string instead of the argument placeholder (e.g.
        '[ ...]') when the action's metavar is empty.
        """
        if action.metavar == "" or action.metavar == ("",):
            return ""
        return super()._format_args(action, default_metavar)

    def _format_action_invocation(self, action):
        """
        Remove the space before commas in option lists, e.g.
        '--cut-both , -cb' -> '--cut-both, -cb'.
        """
        invocation = super()._format_action_invocation(action)
        return invocation.replace(" ,", ",")

    def _fill_text(self, text, width, indent):
        """
        Wrap each '\\n\\n'-separated paragraph of a description separately,
        so manual paragraph breaks are preserved instead of the whole text
        being collapsed into one wrapped block.
        """
        paragraphs = text.split("\n\n")
        return "\n\n".join(
            super(CleanHelpFormatter, self)._fill_text(p, width, indent)
            if p.strip()
            else ""
            for p in paragraphs
        )

    def format_help(self):
        """
        Remove the blank line argparse inserts between a group's description
        and its first option, while keeping the blank lines between groups
        and manual '\\n\\n' breaks within descriptions.
        """
        help_text = super().format_help()
        help_text = re.sub(r"\n\n(?=  -)", "\n", help_text)
        return help_text

def print_adapters():
    """
    Prints all built-in adapter groups and their sequences to stdout,
    formatted and aligned by group. Called for --list-adapters.
    """
    print("\nBuilt-in adapter groups and sequences.")
    print("For matched paired reads, forward reads will only be trimmed for 'Read 1', and reverse reads for 'Read 2' adapters (when applicable).\n")
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
      - --input-paired files are grouped into (R1, R2) tuples
        (group_paired_input_into_pairs()).
      - --output defaults to the current working directory.
      - --gzip turns --stdout off.
      - --stdout turns off verbose output, the progress bar and
        --write-rejected, and turns on --discard-singles and
        --interleaved-out.
      - The start/end settings of set-length, quality-dependent and
        homopolymer end trimming take their specific value if given, and
        otherwise the matching "both" value (--cut-both,
        --endqual-min-both, --poly-length-both, --poly-bases-both).
      - --n-filter means N is no longer counted as a separate letter in
        k-mer complexity filtering ("allow_n_kmer").
      - --adapter-group "all" selects every built-in group.
      - With adapter trimming on: if --adapter-read-1 and/or
        --adapter-read-2 is given, those sequences are used per mate
        ("adapter_per_read") and all other adapter options are ignored.
        Otherwise the adapter list is taken from --adapter-excl alone, or
        else from the selected --adapter-group(s) plus --adapter-add, and
        de-duplicated.
      - The thread count is resolved (worker_determination()).

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
        the corresponding output and exits. With --full-auto, prints a
        warning to stdout that other options are ignored (--verbose is
        reset by full-auto too, so this is always printed). Calls
        parser.error() (which exits with status 2) if no input can be
        found, if --input-paired has an odd number of files, or if piped
        stdin is empty. May create a temporary file holding stdin data.

    Raises:
        ValueError: If an adapter FASTA file is malformed or empty, or an
            adapter option holds neither an existing file nor valid
            sequences (see load_adapters_from_arg()).
    """
    setup_logging()

    parser = argparse.ArgumentParser(
        prog="readzor",
        description="Readzor: a modular FASTQ quality trimming pipeline.\n\n"
                     "All modules are off by default. To use a module, specify a module flag. "
                     "Further specifications with module settings possible.",
        formatter_class=CleanHelpFormatter,
        add_help=False
    )

    general_group = parser.add_argument_group(
        "General settings",
        "General options for help, version information, automatic mode, and console output.")
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
        "Specify input FASTQ files using any combination of --input-files, --input-paired, --input-unpaired, and --input-interleaved. "
        "Any combination of regular (fastq/fq), and gzipped (fastq.gz/fq.gz) files is accepted. "
        "If no input flags are provided, data will be read directly from standard input (stdin). "
        "If --full-auto/-GO is specified, suitable files will be detected in the current working directory."
    )
    input_group.add_argument(
        "--input-files", "-i", nargs='+', default = None, metavar = "",
        help="[STR] FASTQ files of unspecified pairing (or data from stdin). Paired and unpaired files will be auto-detected."
    )
    input_group.add_argument(
        "--input-paired", "-ip", nargs='+', default=None, metavar = "",
        help="[STR] Paired-end FASTQ files, given as one or more R1/R2 pairs, e.g. --input-paired sample1_R1 sample1_R2 sample2_R1 sample2_R2"
    )
    input_group.add_argument(
        "--input-interleaved", "-ii", nargs='+', default = None, metavar = "",
        help="[STR] Interleaved FASTQ files. Note: these files will be split into forward and reverse reads."
    )
    input_group.add_argument(
        "--input-unpaired", "-iu", nargs='+', default = None, metavar = "",
        help="[STR] Unpaired FASTQ files."
    )

    output_group = parser.add_argument_group(
        "Output options",
        "Filtered sequences can be either written to file, or output to stdout. Setting --stdout overrides other settings."
        )
    output_group.add_argument(
        "--gzip", action="store_true", default = False,
        help="[FLAG] Compress filtered FASTQ files in gzip format. Default: off."
    )
    output_group.add_argument(
        "--gzip-level", type=int, default = 1, metavar = "", choices=range(0, 4),
        help="[INT] Set gzip compression level. Higher compression decreases processing speed. Possible values: 0-3. Default: 1."
    )
    output_group.add_argument(
        "--output", "-o", type=str, default = None, metavar = "",
        help="[STR] Path to directory in which the timestamped results folder will be created. Default: current working directory."
    )
    output_group.add_argument(
        "--stdout", action="store_true", default=False,
        help="[FLAG] Stream resulting FASTQ reads to stdout. Forces --interleaved-out and --discard-singles for paired and interleaved files. Overrides --verbose and --progress. Overridden to 'off' by --gzip. Note: all files will be streamed end to end, without any separators. Default: off."
    )
    output_group.add_argument(
        "--interleaved-out", action="store_true", default=False,
        help="[FLAG] Interleave surviving paired FASTQ reads of paired and interleaved input files. Single reads are written to separate files. Default: off."
    )    
    output_group.add_argument(
        "--write-rejected", action="store_true", default=False,
        help="[FLAG] Write rejected reads to file. Either one (for unpaired and when --interleaved-out is set), or two (for forward and reverse reads) are produced. Overridden to 'off' when --stdout is set. Default: off."
    )  
    output_group.add_argument(
        "--discard-singles", action="store_true", default=False,
        help="[FLAG] For paired and interleaved reads, single surviving reads will be discarded instead of written to a separate file. No effect on unpaired read filtering. Default: off."
    )  
    output_group.add_argument(
        "--stats", action="store_true", default=False,
        help="[FLAG] Write detailed summary read trimming statistics. Default: off."
    )

    general_quality_group = parser.add_argument_group(
        "General output filter options",
        "Filters raw input or trimmed output sequences for overall parameters. When a read does not pass these filters, it is discarded entirely."
        )
    general_quality_group.add_argument(
        "--min-length-input", type=int, default = None, metavar = "",
        help="[INT] Minimum length for input read. Default: off."
    )
    general_quality_group.add_argument(
        "--max-length-input", type=int, default = None, metavar = "",
        help="[INT] Maximum length for input read. Default: off."
    )
    general_quality_group.add_argument(
        "--min-length-output", type=int, default = None, metavar = "", 
        help="[INT] Minimum length of output read. Default: off."
    )
    general_quality_group.add_argument(
        "--max-length-output", type = int, default = None, metavar = "", 
        help="[INT] Maximum length of output read. Default: off."
    )
    general_quality_group.add_argument(
        "--min-length-output-perc", type=int, default = None, metavar = "", 
        help="[INT] Minimum length of output read as percentage of its input read. Overridden by --min-length-output. Default: off."
    )
    general_quality_group.add_argument(
        "--max-length-output-perc", type = int, default = None, metavar = "", 
        help="[INT] Maximum length of output read as percentage of its input read. Overridden by --max-length-output. Default: off."
    )
    general_quality_group.add_argument(
        "--min-average-qual-pre", type=int, default = 0, metavar = "", choices=range(0, 128),
        help="[INT] Minimum average quality of input read. Default: off."
    )
    general_quality_group.add_argument(
        "--min-average-qual-post", type=int, default = 0, metavar = "", choices=range(0, 128),
        help="[INT] Minimum average quality of output read. Default: off."
    )
    general_quality_group.add_argument(
        "--n-filter", action="store_true", default=False,
        help="[FLAG] Reject input reads containing N bases. Default: off."
    )

    trim_ends_group = parser.add_argument_group(
        "Set-length end trimming",
        "Trim a defined number of bases from the ends of each read, independent of sequence or quality."
        )
    trim_ends_group.add_argument(
        "--cut-flag", "-cf",action="store_true", default = False,
        help="[FLAG] Turn on set-length end trimming module. Default: off."
    )
    trim_ends_group.add_argument(
        "--cut-start", "-cs", type=int, default = None, metavar="", 
        help="[INT] Number of bases to trim from the start of the read. Default: 0."
    )
    trim_ends_group.add_argument(
        "--cut-end", "-ce", type=int, default = None, metavar="", 
        help="[INT] Number of bases to trim from the end of the read. Default: 0."
    )
    trim_ends_group.add_argument(
        "--cut-both", "-cb", type=int, default = 0, metavar="", 
        help="[INT] Number of bases to trim from both ends of the read. Overridden by --cut-start and --cut-end. Default: 0."
    )

    quality_ends_group = parser.add_argument_group(
        "Quality-dependent end trimming",
        "Trim the ends of each read for quality, up to the first position that fulfills quality requirements."
        )
    quality_ends_group.add_argument(
        "--endqual-filter-flag", "-ef", action="store_true", default = False,
        help="[FLAG] Turn on quality-dependent end trimming. Default: off."
    )
    quality_ends_group.add_argument(
        "--endqual-min-start", "-ems", type = int, default = None, metavar="", choices=range(0, 94),
        help="[INT] Phred score threshold for the start of the read. Default: 20."
    )
    quality_ends_group.add_argument(
        "--endqual-min-end", "-eme", type=int, default = None, metavar="",choices=range(0, 94),
        help="[INT] Phred score threshold for the end of the read. Default: 20."
    )
    quality_ends_group.add_argument(
        "--endqual-min-both", "-emb", type=int, default = 20, metavar="", choices=range(0, 94),
        help="[INT] Phred score threshold for the quality trimming of read ends. Overridden by --endqual-min-start and --endqual-min-end. Default: 20."
    )

    sliding_window_group = parser.add_argument_group(
        "Sliding window quality trimming",
        "Trim the reads for quality based on a sliding window. "
        "In case of mid-read quality dropoff, longest portion survives. In case of tie: portion with highest average quality survives."
        )
    sliding_window_group.add_argument(
        "--slider-filter-flag", "-sf", action="store_true", default = False,
        help="[FLAG] Turn on sliding window quality trimming module. Default: off."
    )
    sliding_window_group.add_argument(
        "--slider-window", "-sw", type = int, default = 5, metavar="",
        help="[INT] Window size over which average quality is calculated. Default: 5."
    )
    sliding_window_group.add_argument(
        "--slider-quality", "-sq", type = int, default = 20, metavar="",
        help="[INT] Minimum average quality in sliding window. Default: 20."
    )
    sliding_window_group.add_argument(
        "--slider-step", "-ss", type = int, default = 1, metavar="",
        help="[INT] Sliding window step size. Default: 1."
    )

    homopolymer_nucleotide_trimming = parser.add_argument_group(
        "Homopolymer nucleotide trimming",
        "Illumina NovaSeq, NextSeq, and MiniSeq use a two-color chemistry, in which guanine bases are unlabeled. "
        "In event of short fragments, this can result in homopolymer G calls at the end of reads. "
        "This filter allows these bases to be cut."
        )
    homopolymer_nucleotide_trimming.add_argument(
        "--poly-filter-flag", "-pf", action="store_true", default = False,
        help="[FLAG] Turn on homopolymer read-end trimming module. Default: off."
    )
    homopolymer_nucleotide_trimming.add_argument(
        "--poly-bases-start", "-pbs", type = str, default = None, metavar = "",
        help='[STR] Base(s) to check for a homopolymer run at start of read. Comma-separated bases are checked independently. Default: none.'
    )
    homopolymer_nucleotide_trimming.add_argument(
        "--poly-bases-end", "-pbe", type = str, default = "G", metavar = "",
        help='[STR] Base(s) to check for a homopolymer run at end of read. Comma-separated bases are checked independently. Default: "G".'
    )
    homopolymer_nucleotide_trimming.add_argument(
        "--poly-bases-both", "-pbb", type = str, default = None, metavar = "",
        help='[STR] Base(s) to check for a homopolymer run at both read ends. Comma-separated bases are checked independently. Overridden by --poly-bases-start and --poly-bases-end. Default: none.'
    )
    homopolymer_nucleotide_trimming.add_argument(
        "--poly-length-start", "-pls", type = int, default = None, metavar = "",
        help="[INT] Minimum length of homopolymer run at start of read required to trigger trimming. Default: 0."
    )
    homopolymer_nucleotide_trimming.add_argument(
        "--poly-length-end", "-ple", type = int, default = 10, metavar = "",
        help="[INT] Minimum length of homopolymer run at end of read required to trigger trimming. Default: 10."
    )
    homopolymer_nucleotide_trimming.add_argument(
        "--poly-length-both", "-plb", type = int, default = None, metavar = "",
        help="[INT] Minimum length of homopolymer run at start and end of read required to trigger trimming. Overridden by --poly-length-start and --poly-length-end. Default: 0."
    )

    adapter_trimming = parser.add_argument_group(
        "Adapter trimming",
        "Trim reads for adapter sequences. "
        "Standard sequences included are TruSeq3, Nextera, and Illumina adapters. "
        "Adapter trimming is performed independent of quality. "
        "Can be combined with overlap trimming (--overlap-filter-flag)."
        )
    adapter_trimming.add_argument(
        "--adapter-filter-flag", "-af", action="store_true", default = False,
        help='[FLAG] Turn on adapter trimming module. Default: off.'
    )
    adapter_trimming.add_argument(
        "--adapter-seed", "-as", type = int, default = 6, metavar="", choices=range(1,100),
        help="[INT] Minimal 5'-end match length. Shorter seeds can result in more partial hits found at end of reads. Default: 6."
    )
    adapter_trimming.add_argument(
        "--adapter-group", "-ag", nargs = "+", default = ["TruSeq"], choices = ["all"] + [name for name, _ in DEFAULT_ADAPTERS], metavar="",
        help='[STR] Specify the group(s) of adapters to be used. Ignored if --adapter-excl is set. Choices: ' + ", ".join(["all"] + [name for name, _ in DEFAULT_ADAPTERS]) + '. For matched paired reads, forward reads will only be trimmed for "Read 1", and reverse reads for "Read 2" adapters (when applicable). Default: TruSeq.'
        )
    adapter_trimming.add_argument(
        "--adapter-mismatch", "-am", type = int, default = 2, metavar="",
        help="[INT] Number of mismatches allowed in adapter finding. Default: 2."
    )
    adapter_trimming.add_argument(
        "--adapter-add", "-ad", type = str, default = None, metavar="",
        help="[STR] FASTA file or comma-separated sequences to add to selected built-in sequences."
    )
    adapter_trimming.add_argument(
        "--adapter-excl", "-ax", type = str, default = None, metavar="",
        help="[STR] FASTA file or comma-separated sequences to trim instead of the built-in sequences."
    )
    adapter_trimming.add_argument(
        "--adapter-read-1", "-a1", type = str, default = None, metavar="",
        help="[STR] FASTA file or comma-separated sequences of read 1 specific adapter. Specifying this will override other adapter specifications. Leaving this empty when --adapter-read-2 is set will disable read 1 adapter trimming. Default: none."
    )
    adapter_trimming.add_argument(
        "--adapter-read-2", "-a2", type = str, default = None, metavar="",
        help="[STR] FASTA file or comma-separated sequences of read 2 specific adapter. Specifying this will override other adapter specifications. Leaving this empty when --adapter-read-1 is set will disable read 2 adapter trimming. Default: none."
    )
    
    overlap_trimming = parser.add_argument_group(
        "Overlap trimming",
        "Perform overlap analysis to trim adapter sequences. "
        "Works independently of adapter sequence. "
        "Only available for paired and interleaved reads. Can be combined with adapter trimming (--adapter-filter-flag)."
        )
    overlap_trimming.add_argument(
        "--overlap-filter-flag", "-of", action="store_true", default = False,
        help="[FLAG] Turn on overlap trimming module. Default: off."
    )
    overlap_trimming.add_argument(
        "--overlap-mismatches", "-om", type = int, default = 5, metavar="",
        help='[INT] Maximum mismatches allowed in the overlapping region. Default: 5.'
    )
    overlap_trimming.add_argument(
        "--overlap-min-length", "-ol", type = int, default = 50, metavar="",
        help="[INT] Minimum overlap length between the paired reads. Default: 50."
    )
    overlap_trimming.add_argument(
        "--overlap-portion-mismatch", "-op", type = int, default = 4, metavar="",
        help="[INT] Maximum percentage of mismatched bases allowed in the overlapping region. Default: 4."
    )

    low_complexity_group = parser.add_argument_group(
        "Low-complexity filtering",
        "Detect complexity of reads using k-mer-based nucleotide frequencies. "
        "Low-complexity reads are discarded entirely."
        )
    low_complexity_group.add_argument(
        "--kmer-filter-flag", "-kf", action="store_true", default=False,
        help="[FLAG] Turn on the k-mer-based complexity filtering module. Default: off."
    )
    low_complexity_group.add_argument(
        "--kmer-size", "-ks", default = 4, metavar="",
        help="[INT] K-mer length for k-mer-based complexity filtering. Comma-separated values are checked independently. Maximum value: 21. Default: 4."
    )
    low_complexity_group.add_argument(
        "--kmer-cutoff", "-kc", type = int, default = 50, metavar="",
        help="[INT] Minimum percentage of unique k-mers (relative to the maximum possible) required to pass the complexity filter. Higher values are stricter. Default: 50."
    )
    
    mgi_convert_group = parser.add_argument_group(
        "MGI header conversion",
        "Convert read header from MGI (BGI) format to Illumina format. "
        "Original header will be discarded."
        )
    mgi_convert_group.add_argument(
        "--mgi-convert-flag", "-mf", action="store_true", default=False,
        help="[FLAG] Turn on the MGI-to-Illumina header conversion module. Default: off."
    )
    mgi_convert_group.add_argument(
        "--mgi-bc5", "-m5", type = str, default = "PLACEHOLDERi5", metavar="",
        help="[STR] Input an i5 barcode for Illumina header conversion. Default: 'PLACEHOLDERi5'."
    )
    mgi_convert_group.add_argument(
        "--mgi-bc7", "-m7", type = str, default = "PLACEHOLDERi7", metavar="",
        help="[STR] Input an i7 barcode for Illumina header conversion. Default: 'PLACEHOLDERi7'."
    )
    mgi_convert_group.add_argument(
        "--mgi-instrument", "-mi", type = str, default = "PLACEHOLDERinstrument", metavar="",
        help="[STR] Instrument name for Illumina header conversion. Default: 'PLACEHOLDERinstrument'."
    )
    mgi_convert_group.add_argument(
        "--mgi-run", "-mr", type = str, default = "PLACEHOLDERrun", metavar="",
        help="[STR] Run ID for Illumina header conversion. Default: 'PLACEHOLDERrun'."
    )

    advanced_group = parser.add_argument_group(
        "Advanced options",
        "Further options that can be specified to alter the behaviour of Readzor.")
    advanced_group.add_argument(
        "--threads", "-t", type = int, default = None, metavar="",
        help="[INT] Number of threads to use. Default: platform-dependent through auto-detection: detection of assigned CPUs on HPC clusters, all-1 otherwise. Fallback: 1."
    )
    advanced_group.add_argument(
        "--reads-for-phred-offset", type = int, default = 500, metavar="",
        help="[INT] Number of reads to sample per file for detection of Phred quality encoding offset. Default: 500."
    )
    advanced_group.add_argument(
        "--chunk-size", type=int, default = 1500, metavar="",
        help="[INT] Number of reads per chunk sent to each worker. Changing can alter processing speed. Default: 1500."
    )
    advanced_group.add_argument(
        "--phred-offset", type = int, choices=[33, 64], default = None, metavar="",
        help="[INT] Define phred offset for all FASTQ files. When set, per-file auto-detection will not be performed. Possible values: 33, 64. Default: off (auto-detection per file)."
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
        help="[INT] Convert Phred encoding from 33 to 64, and vice versa. Possible values: 33, 64. Default: off."
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
            
    def priority_arguments(specific, both):
        if specific is not None:
            return specific
        else:
            return both

    parameters = {}
    parameters["testrun"] = args.testrun
    parameters["full_auto"] = args.full_auto
    parameters["unspecified_files"] = args.input_files
    parameters["unpaired_files"] = args.input_unpaired
    parameters["paired_files"] = group_paired_input_into_pairs(files = args.input_paired, parser = parser)
    parameters["interleaved_files"] = args.input_interleaved
    parameters["minimum_average_qual_post"] = args.min_average_qual_post
    parameters["minimum_average_qual_pre"] = args.min_average_qual_pre    
    parameters["slider_window"] = args.slider_window
    parameters["slider_step"] = args.slider_step
    parameters["slider_quality"] = args.slider_quality
    parameters["gzip_output"] = args.gzip
    parameters["gzip_level"] = args.gzip_level
    parameters["reads_for_phred_offset"] = args.reads_for_phred_offset
    parameters["adapter_filter_flag"] = args.adapter_filter_flag
    parameters["adapter_add"] = args.adapter_add
    parameters["adapter_excl"] = args.adapter_excl
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
    parameters["stats"] = args.stats
    parameters["adapter_group"] = (
        [name for name, _ in DEFAULT_ADAPTERS] if "all" in args.adapter_group
        else list(dict.fromkeys(args.adapter_group))
    )
    parameters["adapter_read_1"] = args.adapter_read_1
    parameters["adapter_read_2"] = args.adapter_read_2
    parameters["discard_singles"] = args.discard_singles
    parameters["overlap_filter_flag"] = args.overlap_filter_flag
    parameters["overlap_mismatches"] = args.overlap_mismatches
    parameters["overlap_min_length"] = args.overlap_min_length
    parameters["overlap_portion_mismatch"] = args.overlap_portion_mismatch
    parameters["cut_start"] = priority_arguments(args.cut_start, args.cut_both)
    parameters["cut_end"] = priority_arguments(args.cut_end, args.cut_both)
    parameters["endqual_min_start"] = priority_arguments(args.endqual_min_start, args.endqual_min_both)
    parameters["endqual_min_end"] = priority_arguments(args.endqual_min_end, args.endqual_min_both)
    parameters["poly_length_start"] = priority_arguments(args.poly_length_start, args.poly_length_both)
    parameters["poly_length_end"] = priority_arguments(args.poly_length_end, args.poly_length_both)
    parameters["poly_bases_start"] = priority_arguments(args.poly_bases_start, args.poly_bases_both)
    parameters["poly_bases_end"] = priority_arguments(args.poly_bases_end, args.poly_bases_both)
    parameters["adapter_per_read"] = None
    
    if parameters["gzip_output"]:
        parameters["stdout"] = False
        
    if parameters["stdout"]:
        parameters["verbose"] = False
        parameters["show_progress"] = False
        parameters["write_rejected"] = False
        parameters["discard_singles"] = True
        parameters["interleaved_out"] = True

    if parameters.get("adapter_filter_flag"):
        read1 = parameters.get("adapter_read_1")
        read2 = parameters.get("adapter_read_2")
        if read1 or read2:
            r1 = [seq.encode("utf-8") for _, seq in load_adapters_from_arg(read1)] if read1 else []
            r2 = [seq.encode("utf-8") for _, seq in load_adapters_from_arg(read2)] if read2 else []
            parameters["adapter_sequences"] = r1 + r2
            parameters["adapter_per_read"] = {"read_1": r1, "read_2": r2}
        else:
            if parameters.get("adapter_excl"):
                raw_adapters = load_adapters_from_arg(parameters["adapter_excl"])
                if raw_adapters == []:
                    raise ValueError(f"No sequences in file '{parameters['adapter_excl']}' detected.")
            else:
                selected_groups = set(parameters["adapter_group"])
                flat_default_adapters = [entry for name, entries in DEFAULT_ADAPTERS if name in selected_groups for entry in entries]
                if parameters.get("adapter_add"):
                    raw_adapters = flat_default_adapters + load_adapters_from_arg(parameters["adapter_add"])
                else:
                    raw_adapters = flat_default_adapters
            parameters["adapter_sequences"] = list(dict.fromkeys(seq.encode("utf-8") for _, seq in raw_adapters))
    else:
        parameters["adapter_sequences"] = []

    parameters["threads"] = worker_determination(parameters["threads"])
    
    return parameters

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

def main():
    """
    Readzor entry point.

    Selects the multiprocessing start method (fork if available, otherwise
    forkserver, otherwise spawn) and parses the command line. With
    --testrun, hands over to test_run(). Otherwise creates the timestamped
    results folder, starts file logging, logs all parameters, processes
    all inputs (input_handler()), writes results_summary.txt (with detailed
    statistics if --stats is set) and prints the completion message. On Ctrl+C, exits with status 130. Temporary stdin
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
        summary_results, collector = input_handler(unspecified_files = parameters["unspecified_files"], unpaired_files = parameters["unpaired_files"], paired_files = parameters["paired_files"], interleaved_files = parameters["interleaved_files"], output_dir = created_output_dir, threads = parameters["threads"], chunk_size = parameters["chunk_size"], show_progress = parameters["show_progress"], stdout = parameters["stdout"], interleaved_out = parameters["interleaved_out"], discard_singles = parameters["discard_singles"], write_rejected = parameters["write_rejected"], parameters = parameters)
        write_summary_and_statistics(summary_results, parameters, output_dir = created_output_dir, collector = collector)
        logger.info("Analysis successfully completed!")
        print_final_message(stdout = parameters["stdout"])
    except KeyboardInterrupt:
        logger.info("Readzor interrupted by user — shutting down.")
        if not parameters["verbose"]:
            print("Readzor interrupted by user — shutting down.", file=sys.stderr)
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
            summary_results, collector = input_handler(unspecified_files = parameters["unspecified_files"], unpaired_files = parameters["unpaired_files"], paired_files = parameters["paired_files"], interleaved_files = parameters["interleaved_files"], output_dir = created_output_dir, threads = parameters["threads"], chunk_size = parameters["chunk_size"], show_progress = parameters["show_progress"], stdout = parameters["stdout"], interleaved_out = parameters["interleaved_out"],discard_singles = parameters["discard_singles"], write_rejected = parameters["write_rejected"], parameters = parameters)
            write_summary_and_statistics(summary_results, parameters, output_dir = created_output_dir, collector = collector)
    finally:
        cleanup_stdin_temp_files()
        logger.info("All files deleted.")
        logger.info("Testrun completed!")
