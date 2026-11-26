#!/usr/bin/env python3
"""
Fast version of judge_trimmed_files.py: judges trimmed FASTQ files against a truth
table that was precomputed once with precompute_truth.py.

Only read names and lengths are used, so no truth sequences are loaded and the
per-read comparison is vectorised with numpy. Summary columns, --by-adapter-length
and the meaning of exact / under-trimmed / over-trimmed / removed are identical to
judge_trimmed_files.py.

Differences from the original:
  * --truth-dir DIR (output of precompute_truth.py) replaces --truth.
  * the "trimmed read is a prefix of the original" check is dropped (needs sequences).
  * --per-read has no kept_seq / cut_seq / seq_matches_truth columns.
  * a read seen twice (e.g. duplicated across files) counts once (last length wins);
    the original counted it twice.
  * trimmed reads whose name is not in the truth are skipped.

Usage:
  judge_trimmed_fast.py --truth-dir truth_pre/ trimmed_R1 trimmed_R2 [trimmed_R1_unpaired trimmed_R2_unpaired ...]
Trimmed files are given in the order R1, R2, R1, R2, ... (an odd last file is treated as R1).
"""
import argparse
import gzip
import os
import sys

import numpy as np

try:
    from isal import igzip as gzip_module
except ImportError:
    gzip_module = gzip

NOT_SEEN = 65535  # trimmed length meaning "read not seen in any trimmed file"
MATES = ("R1", "R2")
SUMMARY_COLUMNS = ["label", "mate", "truth_reads", "reads_with_adapter", "retained_reads", "removed_reads",
                   "exact_reads", "under_trimmed_reads", "over_trimmed_reads", "adapter_bases", "adapter_bases_left",
                   "genomic_bases", "genomic_bases_lost", "genomic_bases_in_removed_reads",
                   "adapter_recall", "precision", "f1", "pct_exact_of_retained"]
STATUS_LABELS = ["correct", "under_trimmed", "over_trimmed", "removed"]


def open_file(path):
    with open(path, "rb") as probe:
        magic_bytes = probe.read(2)
    if magic_bytes == b"\x1f\x8b":
        return gzip_module.open(path, "rb")
    return open(path, "rb", buffering=1 << 20)


def read_names(headers):
    names = []
    for header in headers:
        name = header.split(None, 1)[0].lstrip(b"@")
        if name[-2:] in (b"/1", b"/2"):
            name = name[:-2]
        names.append(name)
    return np.array(names, dtype="S")


def parse_chunk(lines):
    names = read_names(lines[0::4])
    read_lengths = np.fromiter(map(len, lines[1::4]), dtype=np.int64, count=len(lines) // 4)
    return names, read_lengths


def fastq_chunks(path, chunk_size=1 << 26):
    """Yield (names, read_lengths) per ~64 MB of FASTQ."""
    with open_file(path) as file:
        leftover = b""
        while True:
            data = file.read(chunk_size)
            if data and b"\r" in data:
                data = data.replace(b"\r\n", b"\n")
            lines = (leftover + data).split(b"\n")
            if not data:
                if lines and lines[-1] == b"":
                    lines.pop()
                complete_lines = len(lines) // 4 * 4
                if complete_lines:
                    yield parse_chunk(lines[:complete_lines])
                return
            # the last line may be cut off, so it always goes to the next chunk
            complete_lines = (len(lines) - 1) // 4 * 4
            if complete_lines:
                yield parse_chunk(lines[:complete_lines])
            leftover = b"\n".join(lines[complete_lines:])


def ratio(numerator, denominator):
    return float("nan") if denominator == 0 else numerator / denominator


def recall_precision_f1(adapter_bases_removed, adapter_bases, genomic_bases_lost):
    recall = ratio(adapter_bases_removed, adapter_bases)
    precision = ratio(adapter_bases_removed, adapter_bases_removed + genomic_bases_lost)
    if adapter_bases_removed:
        f1 = ratio(2 * precision * recall, precision + recall)
    elif adapter_bases:
        f1 = 0.0
    else:
        f1 = float("nan")
    return recall, precision, f1


def format_value(value):
    if isinstance(value, float):
        return "nan" if value != value else f"{value:.6f}"
    return str(value)


def stats_per_adapter_length(trimmed_lengths, genomic_lengths, adapter_lengths):
    """Per adapter length k: counts of retained / exact / under / over and residual / lost bases."""
    seen = trimmed_lengths != NOT_SEEN
    trimmed = trimmed_lengths.astype(np.int32)
    genomic = genomic_lengths.astype(np.int32)
    difference = trimmed - genomic
    number_of_lengths = int(adapter_lengths.max()) + 1

    def count_per_length(mask, weights=None):
        if weights is not None:
            weights = weights[mask]
        return np.bincount(adapter_lengths[mask], weights=weights, minlength=number_of_lengths)

    under = seen & (difference > 0)
    over = seen & (difference < 0)
    return {
        "truth_reads": np.bincount(adapter_lengths, minlength=number_of_lengths).astype(np.int64),
        "retained": count_per_length(seen).astype(np.int64),
        "exact": count_per_length(seen & (difference == 0)).astype(np.int64),
        "under": count_per_length(under).astype(np.int64),
        "over": count_per_length(over).astype(np.int64),
        "adapter_left": np.rint(count_per_length(under, difference)).astype(np.int64),
        "genomic_lost": np.rint(count_per_length(over, -difference)).astype(np.int64),
        "genomic_total": int(genomic.sum(dtype=np.int64)),
        "genomic_retained": int(genomic[seen].sum(dtype=np.int64)),
    }


def summarise(label, mate, stats):
    adapter_length = np.arange(len(stats["truth_reads"]), dtype=np.int64)
    truth_reads = int(stats["truth_reads"].sum())
    reads_with_adapter = int(stats["truth_reads"][1:].sum())
    retained = int(stats["retained"].sum())
    adapter_bases = int((adapter_length * stats["retained"]).sum())
    adapter_bases_left = int(stats["adapter_left"].sum())
    genomic_bases_lost = int(stats["genomic_lost"].sum())
    recall, precision, f1 = recall_precision_f1(adapter_bases - adapter_bases_left, adapter_bases,
                                                genomic_bases_lost)
    exact = int(stats["exact"].sum())
    return [label, mate, truth_reads, reads_with_adapter, retained, truth_reads - retained, exact,
            int(stats["under"].sum()), int(stats["over"].sum()), adapter_bases, adapter_bases_left,
            stats["genomic_retained"], genomic_bases_lost, stats["genomic_total"] - stats["genomic_retained"],
            recall, precision, f1, ratio(exact, retained)]


def write_per_read(path, names, genomic_lengths, adapter_lengths, trimmed_lengths, errors_only):
    with open(path, "w") as output:
        output.write("\t".join(["mate", "name", "had_adapter", "adapter_k", "expected_len", "trimmed_len",
                                "status", "adapter_left", "genomic_lost"]) + "\n")
        for mate in MATES:
            trimmed = trimmed_lengths[mate].astype(np.int32)
            genomic = genomic_lengths[mate].astype(np.int32)
            seen = trimmed_lengths[mate] != NOT_SEEN
            # status: 0 correct, 1 under-trimmed, 2 over-trimmed, 3 removed
            status = np.where(~seen, 3, np.where(trimmed == genomic, 0, np.where(trimmed > genomic, 1, 2)))
            if errors_only:
                selected = np.flatnonzero(status != 0)
            else:
                selected = np.arange(len(status))
            for index in selected:
                read_status = int(status[index])
                trimmed_length = int(trimmed[index]) if seen[index] else 0
                adapter_left = trimmed_length - int(genomic[index]) if read_status == 1 else 0
                if read_status == 2:
                    genomic_lost = int(genomic[index]) - trimmed_length
                elif read_status == 3:
                    genomic_lost = "NA"
                else:
                    genomic_lost = 0
                adapter_k = int(adapter_lengths[mate][index])
                output.write("\t".join(map(str, [
                    mate, names[index].decode(), "yes" if adapter_k else "no", adapter_k,
                    int(genomic[index]), trimmed_length, STATUS_LABELS[read_status],
                    adapter_left, genomic_lost])) + "\n")


def main():
    arguments = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    arguments.add_argument("--truth-dir", "--truth", dest="truth_dir", required=True, metavar="DIR",
                           help="directory written by precompute_truth.py")
    arguments.add_argument("trimmed", nargs="+", help="trimmed FASTQ files in order R1, R2, R1, R2, ...")
    arguments.add_argument("-o", "--out", help="summary TSV (default: stdout)")
    arguments.add_argument("--by-adapter-length", metavar="FILE", help="write per-k table to FILE")
    arguments.add_argument("--per-read", metavar="FILE", help="write per-read correct/incorrect table to FILE")
    arguments.add_argument("--errors-only", action="store_true", help="with --per-read: only incorrect reads")
    arguments.add_argument("--seed", metavar="N",
                           help="add columns under_k_lt_seed / under_k_ge_seed: under-trimmed reads whose adapter "
                                "length k is < N / >= N (use NA to fill both columns with NA)")
    arguments.add_argument("--label", help="value for the label column (default: name of first trimmed file)")
    arguments.add_argument("--no-header", action="store_true", help="omit the TSV header line")
    arguments_parsed = arguments.parse_args()
    truth_dir = arguments_parsed.truth_dir
    if not os.path.isdir(truth_dir):
        sys.exit(f"ERROR: {truth_dir} is not a directory; run precompute_truth.py on the truth TSV first.")

    names = np.load(os.path.join(truth_dir, "names.npy"), mmap_mode="r")
    genomic_lengths = {}
    adapter_lengths = {}
    for mate in MATES:
        genomic_lengths[mate] = np.load(os.path.join(truth_dir, f"{mate.lower()}_genomic.npy"))
        adapter_lengths[mate] = np.load(os.path.join(truth_dir, f"{mate.lower()}_adapter.npy")).astype(np.int64)
    truth_count = len(names)
    name_width = names.dtype.itemsize
    trimmed_lengths = {}
    for mate in MATES:
        trimmed_lengths[mate] = np.full(truth_count, NOT_SEEN, dtype=np.uint16)

    # Look up every trimmed read in the sorted truth names and store its length.
    for file_index, path in enumerate(arguments_parsed.trimmed):
        mate = MATES[file_index % 2]
        for chunk_names, chunk_lengths in fastq_chunks(path):
            # names longer than the truth names cannot be in the truth
            if chunk_names.dtype.itemsize > name_width:
                fits = np.char.str_len(chunk_names) <= name_width
            else:
                fits = np.ones(len(chunk_names), dtype=bool)
            padded_names = chunk_names.astype(f"S{name_width}")
            indices = np.searchsorted(names, padded_names)
            indices[indices >= truth_count] = truth_count - 1
            found = fits & (names[indices] == padded_names)
            indices = indices[found]
            trimmed_lengths[mate][indices] = np.minimum(chunk_lengths[found], NOT_SEEN - 1)

    label = arguments_parsed.label or os.path.basename(arguments_parsed.trimmed[0])
    stats = {}
    rows = []
    for mate in MATES:
        stats[mate] = stats_per_adapter_length(trimmed_lengths[mate], genomic_lengths[mate], adapter_lengths[mate])
        rows.append(summarise(label, mate, stats[mate]))

    # ALL row: sum the count columns of R1 and R2, then recompute the ratios.
    r1_row, r2_row = rows
    summed = []
    for r1_value, r2_value in zip(r1_row[2:14], r2_row[2:14]):
        summed.append(r1_value + r2_value)
    # summed: 0 truth_reads, 1 reads_with_adapter, 2 retained, 3 removed, 4 exact, 5 under, 6 over,
    #         7 adapter_bases, 8 adapter_bases_left, 9 genomic_bases, 10 genomic_bases_lost, 11 genomic_in_removed
    recall, precision, f1 = recall_precision_f1(summed[7] - summed[8], summed[7], summed[10])
    rows.append([label, "ALL"] + summed + [recall, precision, f1, ratio(summed[4], summed[2])])

    columns = list(SUMMARY_COLUMNS)
    if arguments_parsed.seed is not None:
        columns += ["under_k_lt_seed", "under_k_ge_seed"]
        if arguments_parsed.seed.upper() == "NA":
            for row in rows:
                row += ["NA", "NA"]
        else:
            seed = int(arguments_parsed.seed)
            under_by_seed = []
            for mate in MATES:
                under = stats[mate]["under"]
                under_by_seed.append((int(under[:seed].sum()), int(under[seed:].sum())))
            under_by_seed.append((under_by_seed[0][0] + under_by_seed[1][0],
                                  under_by_seed[0][1] + under_by_seed[1][1]))
            for row, (below_seed, at_least_seed) in zip(rows, under_by_seed):
                row += [below_seed, at_least_seed]

    output = open(arguments_parsed.out, "w") if arguments_parsed.out else sys.stdout
    if not arguments_parsed.no_header:
        output.write("\t".join(columns) + "\n")
    for row in rows:
        output.write("\t".join(format_value(value) for value in row) + "\n")
    if arguments_parsed.out:
        output.close()

    if arguments_parsed.by_adapter_length:
        with open(arguments_parsed.by_adapter_length, "w") as file:
            file.write("\t".join(["label", "mate", "adapter_bases_k", "truth_reads", "retained_reads",
                                  "removed_reads", "exact_reads", "under_trimmed_reads", "over_trimmed_reads",
                                  "mean_adapter_bases_left", "mean_genomic_bases_lost",
                                  "pct_exact_of_retained"]) + "\n")
            for mate in MATES:
                mate_stats = stats[mate]
                for adapter_length in np.flatnonzero((mate_stats["truth_reads"] > 0) | (mate_stats["retained"] > 0)):
                    truth_reads = int(mate_stats["truth_reads"][adapter_length])
                    retained = int(mate_stats["retained"][adapter_length])
                    file.write("\t".join(format_value(value) for value in [
                        label, mate, int(adapter_length), truth_reads, retained, truth_reads - retained,
                        int(mate_stats["exact"][adapter_length]), int(mate_stats["under"][adapter_length]),
                        int(mate_stats["over"][adapter_length]),
                        ratio(int(mate_stats["adapter_left"][adapter_length]), retained),
                        ratio(int(mate_stats["genomic_lost"][adapter_length]), retained),
                        ratio(int(mate_stats["exact"][adapter_length]), retained)]) + "\n")

    if arguments_parsed.per_read:
        write_per_read(arguments_parsed.per_read, names, genomic_lengths, adapter_lengths, trimmed_lengths,
                       arguments_parsed.errors_only)


if __name__ == "__main__":
    main()