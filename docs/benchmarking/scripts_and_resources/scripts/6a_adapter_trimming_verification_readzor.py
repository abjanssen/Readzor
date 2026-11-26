#!/usr/bin/env python3
"""
Run Readzor on paired-end FASTQ files with a grid of adapter- and overlap-trimming
settings, judge every output against the precomputed NEAT truth (precompute_truth.py)
with judge_trimmed_fast.py and append the rows to one TSV. Only that TSV is kept.

Settings (adapter / overlap grids, modes) are edited at the top of this script.
The truth directory must belong to the input FASTQs; with several pairs, run the
script once per pair.

Output: DIR/verify_results.tsv, one row per configuration x sample x mate [R1, R2, ALL].
        Rows are appended; the header is only written when the file is new.
        adapter_bases_removed = adapter bases - adapter_bases_left;
        under_k_lt_seed / under_k_ge_seed = under-trimmed reads whose adapter remnant k
        is shorter than / at least the seed.

Usage:
  adapter_verification_readzor.py -T truth_pre/ [options] R1_a R2_a [R1_b R2_b ...]
"""
import argparse
import glob
import os
import re
import shutil
import subprocess
import sys

# SETTINGS
THREADS = 10
READZOR = "readzor"
OUT_DIR = "trim_verify_out"
JUDGE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "6b_adapter_trimming_verifier.py")

RESULTS_HEADER = ["adapter_seed", "adapter_mismatch", "overlap_min_len", "overlap_mismatches",
                  "overlap_portion_mismatch_pct", "mate", "removed_reads", "exact_reads",
                  "under_trimmed_reads", "over_trimmed_reads", "adapter_bases_removed",
                  "adapter_bases_left", "genomic_bases_lost", "genomic_bases_in_removed_reads",
                  "adapter_recall", "precision", "f1", "under_k_lt_seed", "under_k_ge_seed"]
FASTQ_EXTENSION = re.compile(r"\.(fastq|fq)(\.gz|\.gzip)?$", re.IGNORECASE)

# sequence search trimming
ADAPTER_GROUP = "TruSeq"
SEED_MIN, SEED_MAX = 1, 10
MISMATCH_MIN, MISMATCH_MAX = 0, 10

# overlap search
OVERLAP_MIN, OVERLAP_MAX = 50, 50
OVERLAP_MISMATCH_MIN, OVERLAP_MISMATCH_MAX = 5, 5
OVERLAP_PCT_MIN, OVERLAP_PCT_MAX = 4, 4

# Modes
MODES = "adapter"  # raw, adapter, overlap, combined


def adapter_settings():
    settings = []
    for seed in range(SEED_MIN, SEED_MAX + 1):
        for mismatches in range(MISMATCH_MIN, MISMATCH_MAX + 1):
            settings.append((f"seed{seed}_mm{mismatches}", seed, mismatches))
    return settings


def overlap_settings():
    settings = []
    for overlap in range(OVERLAP_MIN, OVERLAP_MAX + 1, 10):
        for mismatches in range(OVERLAP_MISMATCH_MIN, OVERLAP_MISMATCH_MAX + 1, 2):
            for pct in range(OVERLAP_PCT_MIN, OVERLAP_PCT_MAX + 1, 10):
                settings.append((f"ol{overlap}_mm{mismatches}_p{pct}", overlap, mismatches, pct))
    return settings


def config_list(modes, adapters, overlaps):
    configs = []
    if "adapter" in modes:
        for adapter in adapters:
            configs.append((f"adapter_{adapter[0]}", "adapter", adapter, None))
    if "overlap" in modes:
        for overlap in overlaps:
            configs.append((f"overlap_{overlap[0]}", "overlap", None, overlap))
    if "combined" in modes:
        for adapter in adapters:
            for overlap in overlaps:
                configs.append((f"combined_{adapter[0]}_{overlap[0]}", "combined", adapter, overlap))
    return configs


def judge_and_append(prefix, seed, files, truth, judge, results_tsv):
    """Judge files given in R1,R2,R1,R2 order and append rows (R1, R2, ALL) to results_tsv."""
    command = [sys.executable, judge, "--truth-dir", truth, *files, "--no-header", "--seed", str(seed)]
    summary = subprocess.run(command, stdout=subprocess.PIPE, text=True, check=True).stdout

    with open(results_tsv, "a") as output:
        for line in summary.splitlines():
            fields = line.split("\t")
            adapter_bases_removed = str(int(fields[9]) - int(fields[10]))
            kept = [fields[1], fields[5], fields[6], fields[7], fields[8], adapter_bases_removed,
                    fields[10], fields[12], fields[13], fields[14], fields[15], fields[16],
                    fields[18], fields[19]]
            output.write("\t".join(prefix + kept) + "\n")


def run_config(config, r1, r2, sample_dir, readzor_command, threads, truth, judge, results_tsv):
    config_id, mode, adapter, overlap = config
    readzor_dir = os.path.join(sample_dir, config_id)

    readzor_args = []
    if adapter is not None:
        _, adapter_seed, adapter_mismatches = adapter
        readzor_args += ["-af", "-as", str(adapter_seed), "-am", str(adapter_mismatches),
                         "-ag", *ADAPTER_GROUP.split()]
    else:
        adapter_seed = adapter_mismatches = "NA"
    if overlap is not None:
        _, overlap_min_length, overlap_mismatches, overlap_pct = overlap
        readzor_args += ["-of", "-ol", str(overlap_min_length), "-om", str(overlap_mismatches),
                         "-op", str(overlap_pct)]
    else:
        overlap_min_length = overlap_mismatches = overlap_pct = "NA"

    os.makedirs(readzor_dir, exist_ok=True)
    command = [*readzor_command, "-ip", r1, r2, "-o", readzor_dir, "--threads", str(threads), *readzor_args]
    subprocess.run(command, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)

    results_dir = glob.glob(os.path.join(readzor_dir, "Readzor_results_*"))[0]
    files = []
    for pattern in ("*_R1_paired_filtered.fastq", "*_R2_paired_filtered.fastq",
                    "*_R1_unpaired_filtered.fastq", "*_R2_unpaired_filtered.fastq"):
        files.append(glob.glob(os.path.join(results_dir, pattern))[0])

    prefix = [str(adapter_seed), str(adapter_mismatches),
              str(overlap_min_length), str(overlap_mismatches), str(overlap_pct)]
    judge_and_append(prefix, adapter_seed, files, truth, judge, results_tsv)

    shutil.rmtree(readzor_dir)


def main():
    arguments = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    arguments.add_argument("fastqs", nargs="+", metavar="FASTQ", help="R1 R2 pairs")
    arguments.add_argument("-T", dest="truth", required=True, help="precomputed truth (precompute_truth.py)")
    arguments.add_argument("-o", dest="out_dir", default=OUT_DIR, help=f"output directory (default: {OUT_DIR})")
    arguments.add_argument("-t", dest="threads", type=int, default=THREADS, help=f"Readzor threads (default: {THREADS})")
    arguments.add_argument("-r", dest="readzor", default=READZOR,
                           help=f'Readzor command (default: "{READZOR}"; e.g. "python3 /path/readzor.py")')
    arguments.add_argument("-j", dest="judge", default=JUDGE, help="judge_trimmed_fast.py (default: next to this script)")
    arguments.add_argument("-m", dest="modes", default=MODES, help=f"modes to run, comma-separated (default: {MODES})")
    arguments_parsed = arguments.parse_args()

    modes = set(arguments_parsed.modes.split(","))
    configs = config_list(modes, adapter_settings(), overlap_settings())
    readzor_command = arguments_parsed.readzor.split()
    fastqs = arguments_parsed.fastqs
    out_dir = arguments_parsed.out_dir
    results_tsv = os.path.join(out_dir, "verify_results.tsv")
    work_dir = os.path.join(out_dir, "work")
    os.makedirs(out_dir, exist_ok=True)

    if not os.path.exists(results_tsv):
        with open(results_tsv, "w") as file:
            file.write("\t".join(RESULTS_HEADER) + "\n")

    for pair_index, (r1, r2) in enumerate(zip(fastqs[0::2], fastqs[1::2]), 1):
        sample = FASTQ_EXTENSION.sub("", os.path.basename(r1))
        sample_dir = os.path.join(work_dir, f"p{pair_index}_{sample}")

        if "raw" in modes:
            judge_and_append(["NA"] * 5, "NA", [r1, r2], arguments_parsed.truth, arguments_parsed.judge, results_tsv)

        for config in configs:
            run_config(config, r1, r2, sample_dir, readzor_command, arguments_parsed.threads,
                       arguments_parsed.truth, arguments_parsed.judge, results_tsv)

    shutil.rmtree(work_dir, ignore_errors=True)


if __name__ == "__main__":
    main()