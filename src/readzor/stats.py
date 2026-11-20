"""Per-file read statistics and the results_summary.txt grid.

Workers collect small, fixed-size numeric aggregates per chunk (sums, counts,
histograms) and return them with their results. The main process merges them
per input and mate. write_summary() then writes results_summary.txt as one
tab-separated table: one row per input and mate, with read counts, length
and quality before/after trimming, one group of columns per trimming module
and one per adapter sequence. Cells that do not apply are "NA".
"""

import csv

import numpy as np

QUAL_BINS_PER_UNIT = 10
FIELDS = ("invalid", "reads_in", "bases_in", "qual_sum_in",
          "reads_out", "bases_out", "qual_sum_out", "rejected_after_trim",
          "read_qual_sum_in", "read_qual_sum_out")
HISTOGRAMS = ("len_hist_in", "len_hist_out", "rq_hist_in", "rq_hist_out")
MODULE_COUNTERS = ("reads", "bases", "decisive_reads")
NA = "NA"
BASE_COLUMNS = [
    "input", "type", "mate", "input_files", "output_files",
    "reads_pre", "invalid_reads", "rejected_after_trimming", "reads_lost",
    "reads_post", "reads_kept_percent", "kept_in_pairs", "kept_as_singletons",
    "bases_pre", "bases_post", "bases_kept_percent",
    "mean_length_pre", "mean_length_post", "median_length_pre", "median_length_post",
    "mean_base_qual_pre", "mean_base_qual_post",
    "mean_read_qual_pre", "mean_read_qual_post", "median_read_qual_pre", "median_read_qual_post",
]
MODULE_FIELDS = ["reads_trimmed", "reads_trimmed_percent", "bases_trimmed", "bases_trimmed_percent", "decisive_reads"]

def new_read_stats():
    """Empty statistics for one mate of an input (or an unpaired input)."""
    stats = {field: 0 for field in FIELDS}
    for hist_name in HISTOGRAMS:
        stats[hist_name] = np.zeros(0, dtype=np.int64)
    stats["modules"] = {}
    stats["adapter_hits"] = {}
    return stats

def add_histogram(histogram_1, histogram_2):
    """
    Add two histograms position by position. They may differ in length;
    the result has the length of the longer one. Neither input is changed.
    """
    if len(histogram_1) < len(histogram_2):
        histogram_1, histogram_2 = histogram_2, histogram_1
    combined = histogram_1.copy()
    combined[:len(histogram_2)] += histogram_2
    return combined

def read_qual_histogram(read_mean_quals):
    """Histogram of per-read mean qualities, in bins of 1 / QUAL_BINS_PER_UNIT Phred."""
    bins = np.floor(read_mean_quals * QUAL_BINS_PER_UNIT + 1e-9).astype(np.int64)
    return np.bincount(bins)

def add_invalid(stats, n_invalid):
    """Count records rejected by validate_fastq()."""
    stats["invalid"] += int(n_invalid)

def collect_input(stats, quality_arr, raw_lengths):
    """Length/quality of the valid input reads. Padding positions are 0 in
    quality_arr (see qual_to_array()), so plain row sums are exact."""
    qual_sums_per_read = quality_arr.sum(axis=1, dtype=np.int32)
    read_mean_quals = qual_sums_per_read / np.maximum(raw_lengths, 1)
    stats["reads_in"] += int(len(raw_lengths))
    stats["bases_in"] += int(raw_lengths.sum())
    stats["qual_sum_in"] += int(qual_sums_per_read.sum(dtype=np.int64))
    stats["read_qual_sum_in"] += float(read_mean_quals.sum())
    stats["len_hist_in"] = add_histogram(stats["len_hist_in"], np.bincount(raw_lengths))
    stats["rq_hist_in"] = add_histogram(stats["rq_hist_in"], read_qual_histogram(read_mean_quals))

def collect_modules(stats, module_names, module_lefts, module_rights, real_lengths, final_lefts, final_rights):
    """
    Per-module trimming statistics.

    Modules are merged by max(left)/min(right), so their effects overlap.
    Two numbers are kept per module:
      - reads/bases: what the module would remove on its own;
      - decisive_reads: reads for which this module set a final boundary
        (several modules can be decisive for the same read).

    Args:
        stats (dict): Statistics of one mate, updated in place.
        module_names (list[str]): Name of each module, in pipeline order.
        module_lefts (list[numpy.ndarray]): Per module, its left boundary per read.
        module_rights (list[numpy.ndarray]): Per module, its right boundary per read.
        real_lengths (numpy.ndarray): Untrimmed length of each read.
        final_lefts (numpy.ndarray): Merged left boundary per read.
        final_rights (numpy.ndarray): Merged right boundary per read.
    """
    real_lengths = np.asarray(real_lengths, dtype=np.int64)
    final_lefts = np.clip(np.asarray(final_lefts, dtype=np.int64), 0, real_lengths)
    final_rights = np.clip(np.asarray(final_rights, dtype=np.int64), 0, real_lengths)
    for module_name, lefts, rights in zip(module_names, module_lefts, module_rights):
        lefts = np.clip(np.asarray(lefts, dtype=np.int64), 0, real_lengths)
        rights = np.clip(np.asarray(rights, dtype=np.int64), 0, real_lengths)
        bases_removed = np.minimum(lefts + (real_lengths - rights), real_lengths)
        trimmed = bases_removed > 0
        sets_final_cut = (((lefts > 0) & (lefts == final_lefts))
                          | ((rights < real_lengths) & (rights == final_rights)))
        module_counts = stats["modules"].setdefault(module_name, dict.fromkeys(MODULE_COUNTERS, 0))
        module_counts["reads"] += int(trimmed.sum())
        module_counts["bases"] += int(bases_removed.sum())
        module_counts["decisive_reads"] += int((sets_final_cut & trimmed).sum())

def collect_adapter_hits(stats, adapter_sequences, hit_index):
    """
    Count, per adapter sequence, the reads it was cut at.

    Args:
        stats (dict): Statistics of one mate, updated in place.
        adapter_sequences (list[bytes]): The adapters searched for in these
            reads (after per-mate selection). Each is recorded, with 0 hits
            if not found, so that the summary can tell "not found" (0) from
            "not searched for" (NA).
        hit_index (numpy.ndarray): Per read, index into adapter_sequences of
            the adapter that set the cut, or -1 (from adapter_trimming()).
    """
    hit_adapter_indices = hit_index[hit_index >= 0]
    if hit_adapter_indices.size:
        reads_per_adapter = np.bincount(hit_adapter_indices, minlength=len(adapter_sequences))
    else:
        reads_per_adapter = np.zeros(len(adapter_sequences), dtype=np.int64)
    for adapter_seq, read_count in zip(adapter_sequences, reads_per_adapter):
        stats["adapter_hits"][adapter_seq] = stats["adapter_hits"].get(adapter_seq, 0) + int(read_count)

def collect_output(stats, quality_arr, final_lefts, final_rights, keep_mask):
    """Length/quality of the kept reads, within their final [left, right) window."""
    keep_mask = np.asarray(keep_mask, dtype=bool)
    kept_rows = np.flatnonzero(keep_mask)
    stats["rejected_after_trim"] += int(keep_mask.size - kept_rows.size)
    if kept_rows.size == 0:
        return
    kept_quals = quality_arr[kept_rows]
    width = kept_quals.shape[1]
    index_dtype = np.int16 if width < 32000 else np.int32
    kept_lefts = np.asarray(final_lefts)[kept_rows]
    kept_rights = np.asarray(final_rights)[kept_rows]
    positions = np.arange(width, dtype=index_dtype)
    in_window = positions >= kept_lefts.astype(index_dtype)[:, None]
    in_window &= positions < kept_rights.astype(index_dtype)[:, None]
    qual_sums_per_read = np.where(in_window, kept_quals, 0).sum(axis=1, dtype=np.int32)
    kept_lengths = kept_rights.astype(np.int64) - kept_lefts.astype(np.int64)
    read_mean_quals = qual_sums_per_read / kept_lengths
    stats["reads_out"] += int(kept_rows.size)
    stats["bases_out"] += int(kept_lengths.sum())
    stats["qual_sum_out"] += int(qual_sums_per_read.sum(dtype=np.int64))
    stats["read_qual_sum_out"] += float(read_mean_quals.sum())
    stats["len_hist_out"] = add_histogram(stats["len_hist_out"], np.bincount(kept_lengths))
    stats["rq_hist_out"] = add_histogram(stats["rq_hist_out"], read_qual_histogram(read_mean_quals))

def merge_read_stats(totals, part):
    """Add the statistics ``part`` (one chunk) into ``totals``, in place. Returns ``totals``."""
    for field in FIELDS:
        totals[field] += part[field]
    for hist_name in HISTOGRAMS:
        totals[hist_name] = add_histogram(totals[hist_name], part[hist_name])
    for module_name, part_counts in part["modules"].items():
        total_counts = totals["modules"].setdefault(module_name, dict.fromkeys(MODULE_COUNTERS, 0))
        for counter in MODULE_COUNTERS:
            total_counts[counter] += part_counts[counter]
    for adapter_seq, read_count in part["adapter_hits"].items():
        totals["adapter_hits"][adapter_seq] = totals["adapter_hits"].get(adapter_seq, 0) + read_count
    return totals

def new_chunk_report(mates):
    """
    Per-chunk report sent from a worker back to the main process. Use in
    unified_worker():

        report = new_chunk_report(("R1", "R2"))
        result = process_paired_task_flat(task, parameters, report["mates"])
        return result + (report,)
    """
    return {"mates": {mate: new_read_stats() for mate in mates}}

class StatsCollector:
    """Main-process side: merges chunk reports per input."""

    def __init__(self):
        self.inputs = {}

    def add(self, label, report, input_files=None):
        entry = self.inputs.setdefault(label, {"input_files": input_files or label, "mates": {}})
        for mate, mate_stats in report["mates"].items():
            if mate in entry["mates"]:
                merge_read_stats(entry["mates"][mate], mate_stats)
            else:
                entry["mates"][mate] = mate_stats

def divide_NA(numerator, denominator, digits=2):
    """numerator / denominator rounded to ``digits``, or NA if the denominator is 0."""
    return round(numerator / denominator, digits) if denominator else NA

def histogram_median(histogram, bins_per_unit=1):
    """Median of the values counted in ``histogram`` (index = value * bins_per_unit), or NA if empty."""
    n_values = int(histogram.sum())
    if n_values == 0:
        return NA
    median_bin = int(np.searchsorted(np.cumsum(histogram), (n_values + 1) / 2))
    return median_bin if bins_per_unit == 1 else round(median_bin / bins_per_unit, 1)

def resolve_input_type(counts):
    if "type" in counts:
        return counts["type"]
    if "kept" in counts:
        return "unpaired"
    return "paired" if ", " in counts.get("input", "") else "interleaved"

def build_summary_grid(file_stats, collector=None, display_name=None):
    """
    One row per input and mate, combining input_handler()'s kept/rejected
    counts (``file_stats``) with the read statistics (``collector``, keyed
    the same way as file_stats). Module and adapter columns are added for
    every module/adapter seen in the run; cells that do not apply are "NA".

    Args:
        file_stats (dict): Output of input_handler(): unpaired inputs keyed
            by filepath, paired/interleaved inputs keyed by name prefix.
            An optional "type" entry ("unpaired", "paired", "interleaved")
            is used if present.
        collector (StatsCollector | None): Merged statistics, or None to
            fill all statistics columns with NA.
        display_name (callable | None): Maps an unpaired file_stats key
            (a filepath) to the name shown in the "input" column.

    Returns:
        tuple[list[str], list[dict]]: Column names and rows.
    """
    stats_by_input = collector.inputs if collector is not None else {}
    module_names = []
    adapter_seqs = []
    for entry in stats_by_input.values():
        for mate_stats in entry["mates"].values():
            module_names += [name for name in mate_stats["modules"] if name not in module_names]
            adapter_seqs += [seq for seq in mate_stats["adapter_hits"] if seq not in adapter_seqs]
    adapter_names = [seq.decode() if isinstance(seq, bytes) else seq for seq in adapter_seqs]

    columns = list(BASE_COLUMNS)
    columns += [f"{module_name}_{field}" for module_name in module_names for field in MODULE_FIELDS]
    columns += [f"adapter_{adapter_name}_reads" for adapter_name in adapter_names]

    rows = []
    for input_key, counts in file_stats.items():
        input_type = resolve_input_type(counts)
        mates = ("all",) if input_type == "unpaired" else ("R1", "R2")
        label = display_name(input_key) if (input_type == "unpaired" and display_name) else input_key
        stats_per_mate = stats_by_input.get(input_key, {}).get("mates", {})
        for mate in mates:
            row = dict.fromkeys(columns, NA)
            row.update(input=label, type=input_type, mate=mate,
                       input_files=counts.get("input", NA), output_files=counts.get("output", NA))
            if input_type == "unpaired":
                kept, rejected = counts["kept"], counts["rejected"]
            else:
                singletons = counts[f"kept_{mate}_singletons"]
                kept = counts["kept_pairs"] + singletons
                rejected = counts[f"rejected_{mate}"]
                row.update(kept_in_pairs=counts["kept_pairs"], kept_as_singletons=singletons)
            row.update(reads_pre=kept + rejected, reads_lost=rejected, reads_post=kept,
                       reads_kept_percent=divide_NA(100 * kept, kept + rejected))

            stats = stats_per_mate.get(mate)
            if stats is not None:
                row.update(
                    invalid_reads=stats["invalid"],
                    rejected_after_trimming=stats["rejected_after_trim"],
                    bases_pre=stats["bases_in"],
                    bases_post=stats["bases_out"],
                    bases_kept_percent=divide_NA(100 * stats["bases_out"], stats["bases_in"], 4),
                    mean_length_pre=divide_NA(stats["bases_in"], stats["reads_in"]),
                    mean_length_post=divide_NA(stats["bases_out"], stats["reads_out"]),
                    median_length_pre=histogram_median(stats["len_hist_in"]),
                    median_length_post=histogram_median(stats["len_hist_out"]),
                    mean_base_qual_pre=divide_NA(stats["qual_sum_in"], stats["bases_in"]),
                    mean_base_qual_post=divide_NA(stats["qual_sum_out"], stats["bases_out"]),
                    mean_read_qual_pre=divide_NA(stats["read_qual_sum_in"], stats["reads_in"]),
                    mean_read_qual_post=divide_NA(stats["read_qual_sum_out"], stats["reads_out"]),
                    median_read_qual_pre=histogram_median(stats["rq_hist_in"], QUAL_BINS_PER_UNIT),
                    median_read_qual_post=histogram_median(stats["rq_hist_out"], QUAL_BINS_PER_UNIT),
                )
                for module_name, module_counts in stats["modules"].items():
                    row[f"{module_name}_reads_trimmed"] = module_counts["reads"]
                    row[f"{module_name}_reads_trimmed_percent"] = divide_NA(100 * module_counts["reads"], stats["reads_in"])
                    row[f"{module_name}_bases_trimmed"] = module_counts["bases"]
                    row[f"{module_name}_bases_trimmed_percent"] = divide_NA(100 * module_counts["bases"], stats["bases_in"])
                    row[f"{module_name}_decisive_reads"] = module_counts["decisive_reads"]
                for adapter_name, adapter_seq in zip(adapter_names, adapter_seqs):
                    row[f"adapter_{adapter_name}_reads"] = stats["adapter_hits"].get(adapter_seq, NA)
            rows.append(row)
    return columns, rows

def write_summary(path, file_stats, collector=None, display_name=None):
    """Write the summary grid to ``path`` as one tab-separated table."""
    columns, rows = build_summary_grid(file_stats, collector, display_name)
    with open(path, "w", encoding="utf-8", newline="") as summary_file:
        writer = csv.DictWriter(summary_file, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)