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

QUAL_BINS_PER_UNIT = 10  # read-mean-quality histogram resolution: 0.1 Phred

_SCALARS = ("invalid", "reads_in", "bases_in", "qual_sum_in",
            "reads_out", "bases_out", "qual_sum_out", "rejected_after_trim",
            "read_qual_sum_in", "read_qual_sum_out")
_HISTS = ("len_hist_in", "len_hist_out", "rq_hist_in", "rq_hist_out")


# ---------------------------------------------------------------------------
# Per-mate read statistics (built in the workers)
# ---------------------------------------------------------------------------

def new_read_stats():
    """Empty statistics for one mate (or an unpaired input)."""
    stats = {k: 0 for k in _SCALARS}
    for h in _HISTS:
        stats[h] = np.zeros(0, dtype=np.int64)
    stats["modules"] = {}       # name -> {"reads", "bases", "decisive_reads"}
    stats["adapter_hits"] = {}  # adapter sequence (bytes) -> reads
    return stats


def _add_hist(a, b):
    if len(a) < len(b):
        a, b = b, a
    a = a.copy()
    a[:len(b)] += b
    return a


def _read_qual_hist(read_means):
    return np.bincount(np.floor(read_means * QUAL_BINS_PER_UNIT + 1e-9).astype(np.int64))


def add_invalid(stats, n_invalid):
    """Count records rejected by validate_fastq()."""
    stats["invalid"] += int(n_invalid)


def collect_input(stats, quality_arr, raw_lengths):
    """Length/quality of the valid input reads. Padding positions are 0 in
    quality_arr (see qual_to_array()), so plain row sums are exact."""
    row_sums = quality_arr.sum(axis=1, dtype=np.int32)
    read_means = row_sums / np.maximum(raw_lengths, 1)
    stats["reads_in"] += int(len(raw_lengths))
    stats["bases_in"] += int(raw_lengths.sum())
    stats["qual_sum_in"] += int(row_sums.sum(dtype=np.int64))
    stats["read_qual_sum_in"] += float(read_means.sum())
    stats["len_hist_in"] = _add_hist(stats["len_hist_in"], np.bincount(raw_lengths))
    stats["rq_hist_in"] = _add_hist(stats["rq_hist_in"], _read_qual_hist(read_means))


def collect_modules(stats, names, lefts, rights, real_lengths, final_left, final_right):
    """
    Per-module trimming statistics.

    Modules are merged by max(left)/min(right), so their effects overlap.
    Two numbers are kept per module:
      - reads/bases: what the module would remove on its own;
      - decisive_reads: reads for which this module set a final boundary
        (several modules can be decisive for the same read).
    """
    real = np.asarray(real_lengths, dtype=np.int64)
    fl = np.clip(np.asarray(final_left, dtype=np.int64), 0, real)
    fr = np.clip(np.asarray(final_right, dtype=np.int64), 0, real)
    for name, left, right in zip(names, lefts, rights):
        l = np.clip(np.asarray(left, dtype=np.int64), 0, real)
        r = np.clip(np.asarray(right, dtype=np.int64), 0, real)
        removed = np.minimum(l + (real - r), real)
        affected = removed > 0
        decisive = ((l > 0) & (l == fl)) | ((r < real) & (r == fr))
        m = stats["modules"].setdefault(name, {"reads": 0, "bases": 0, "decisive_reads": 0})
        m["reads"] += int(affected.sum())
        m["bases"] += int(removed.sum())
        m["decisive_reads"] += int((decisive & affected).sum())


def collect_adapter_hits(stats, adapter_sequences, hit_index):
    """
    Count, per adapter sequence, the reads it was cut at.

    Args:
        adapter_sequences (list[bytes]): The adapters searched for in these
            reads (after per-mate selection). Each is recorded, with 0 hits
            if not found, so that the summary can tell "not found" (0) from
            "not searched for" (NA).
        hit_index (numpy.ndarray): Per read, index into adapter_sequences of
            the adapter that set the cut, or -1 (from adapter_trimming()).
    """
    hits = hit_index[hit_index >= 0]
    counts = np.bincount(hits, minlength=len(adapter_sequences)) if hits.size else np.zeros(len(adapter_sequences), dtype=np.int64)
    for seq, c in zip(adapter_sequences, counts):
        stats["adapter_hits"][seq] = stats["adapter_hits"].get(seq, 0) + int(c)


def collect_output(stats, quality_arr, lefts, rights, keep_mask):
    """Length/quality of the kept reads, within their final window."""
    keep_mask = np.asarray(keep_mask, dtype=bool)
    keep = np.flatnonzero(keep_mask)
    stats["rejected_after_trim"] += int(keep_mask.size - keep.size)
    if keep.size == 0:
        return
    q = quality_arr[keep]
    width = q.shape[1]
    idx_dtype = np.int16 if width < 32000 else np.int32
    l = np.asarray(lefts)[keep]
    r = np.asarray(rights)[keep]
    cols = np.arange(width, dtype=idx_dtype)
    window = cols >= l.astype(idx_dtype)[:, None]
    window &= cols < r.astype(idx_dtype)[:, None]
    row_sums = np.where(window, q, 0).sum(axis=1, dtype=np.int32)
    lengths = r.astype(np.int64) - l.astype(np.int64)
    read_means = row_sums / lengths
    stats["reads_out"] += int(keep.size)
    stats["bases_out"] += int(lengths.sum())
    stats["qual_sum_out"] += int(row_sums.sum(dtype=np.int64))
    stats["read_qual_sum_out"] += float(read_means.sum())
    stats["len_hist_out"] = _add_hist(stats["len_hist_out"], np.bincount(lengths))
    stats["rq_hist_out"] = _add_hist(stats["rq_hist_out"], _read_qual_hist(read_means))


def merge_read_stats(total, part):
    for k in _SCALARS:
        total[k] += part[k]
    for h in _HISTS:
        total[h] = _add_hist(total[h], part[h])
    for name, m in part["modules"].items():
        t = total["modules"].setdefault(name, {"reads": 0, "bases": 0, "decisive_reads": 0})
        for k in t:
            t[k] += m[k]
    for seq, c in part["adapter_hits"].items():
        total["adapter_hits"][seq] = total["adapter_hits"].get(seq, 0) + c
    return total


# ---------------------------------------------------------------------------
# Chunk reports (worker -> main)
# ---------------------------------------------------------------------------

def new_chunk_report(mates):
    """
    Per-chunk report sent from a worker back to the main process. Use in
    unified_worker():

        report = new_chunk_report(("R1", "R2"))
        result = process_paired_task_flat(task, parameters, report["mates"])
        return result + (report,)
    """
    return {"mates": {m: new_read_stats() for m in mates}}


class StatsCollector:
    """Main-process side: merges chunk reports per input."""

    def __init__(self):
        self.inputs = {}   # label -> {"input_files", "mates"}

    def add(self, label, report, input_files=None):
        entry = self.inputs.setdefault(label, {"input_files": input_files or label, "mates": {}})
        for mate, st in report["mates"].items():
            if mate in entry["mates"]:
                merge_read_stats(entry["mates"][mate], st)
            else:
                entry["mates"][mate] = st


# ---------------------------------------------------------------------------
# Summary grid
# ---------------------------------------------------------------------------

NA = "NA"

BASE_COLUMNS = [
    "input", "type", "mate", "input_files", "output_files",
    "reads_pre", "invalid_reads", "rejected_after_trimming", "reads_rejected",
    "reads_post", "reads_kept_percent", "kept_in_pairs", "kept_as_singletons",
    "bases_pre", "bases_post", "bases_kept_percent",
    "mean_length_pre", "mean_length_post", "median_length_pre", "median_length_post",
    "mean_base_qual_pre", "mean_base_qual_post",
    "mean_read_qual_pre", "mean_read_qual_post", "median_read_qual_pre", "median_read_qual_post",
]
MODULE_FIELDS = ["reads_trimmed", "reads_trimmed_percent", "bases_trimmed", "bases_trimmed_percent", "decisive_reads"]


def _div(a, b, digits=2):
    return round(a / b, digits) if b else NA


def _hist_median(h, scale=1):
    total = int(h.sum())
    if total == 0:
        return NA
    idx = int(np.searchsorted(np.cumsum(h), (total + 1) / 2))
    return idx if scale == 1 else round(idx / scale, 1)


def _input_type(counts):
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
    stats_by_key = collector.inputs if collector is not None else {}
    modules, adapters = [], []
    for entry in stats_by_key.values():
        for st in entry["mates"].values():
            modules += [m for m in st["modules"] if m not in modules]
            adapters += [a for a in st["adapter_hits"] if a not in adapters]
    adapter_names = [a.decode() if isinstance(a, bytes) else a for a in adapters]

    columns = list(BASE_COLUMNS)
    columns += [f"{m}_{f}" for m in modules for f in MODULE_FIELDS]
    columns += [f"adapter_{name}_reads" for name in adapter_names]

    rows = []
    for key, counts in file_stats.items():
        kind = _input_type(counts)
        mates = ("all",) if kind == "unpaired" else ("R1", "R2")
        label = display_name(key) if (kind == "unpaired" and display_name) else key
        mate_stats = stats_by_key.get(key, {}).get("mates", {})
        for mate in mates:
            row = dict.fromkeys(columns, NA)
            row.update(input=label, type=kind, mate=mate,
                       input_files=counts.get("input", NA), output_files=counts.get("output", NA))
            if kind == "unpaired":
                kept, rejected = counts["kept"], counts["rejected"]
            else:
                singles = counts[f"kept_{mate}_singletons"]
                kept = counts["kept_pairs"] + singles
                rejected = counts[f"rejected_{mate}"]
                row.update(kept_in_pairs=counts["kept_pairs"], kept_as_singletons=singles)
            row.update(reads_pre=kept + rejected, reads_rejected=rejected, reads_post=kept,
                       reads_kept_percent=_div(100 * kept, kept + rejected))

            s = mate_stats.get(mate)
            if s is not None:
                row.update(
                    invalid_reads=s["invalid"],
                    rejected_after_trimming=s["rejected_after_trim"],
                    bases_pre=s["bases_in"],
                    bases_post=s["bases_out"],
                    bases_kept_percent=_div(100 * s["bases_out"], s["bases_in"], 4),
                    mean_length_pre=_div(s["bases_in"], s["reads_in"]),
                    mean_length_post=_div(s["bases_out"], s["reads_out"]),
                    median_length_pre=_hist_median(s["len_hist_in"]),
                    median_length_post=_hist_median(s["len_hist_out"]),
                    mean_base_qual_pre=_div(s["qual_sum_in"], s["bases_in"]),
                    mean_base_qual_post=_div(s["qual_sum_out"], s["bases_out"]),
                    mean_read_qual_pre=_div(s["read_qual_sum_in"], s["reads_in"]),
                    mean_read_qual_post=_div(s["read_qual_sum_out"], s["reads_out"]),
                    median_read_qual_pre=_hist_median(s["rq_hist_in"], QUAL_BINS_PER_UNIT),
                    median_read_qual_post=_hist_median(s["rq_hist_out"], QUAL_BINS_PER_UNIT),
                )
                for m, v in s["modules"].items():
                    row[f"{m}_reads_trimmed"] = v["reads"]
                    row[f"{m}_reads_trimmed_percent"] = _div(100 * v["reads"], s["reads_in"])
                    row[f"{m}_bases_trimmed"] = v["bases"]
                    row[f"{m}_bases_trimmed_percent"] = _div(100 * v["bases"], s["bases_in"])
                    row[f"{m}_decisive_reads"] = v["decisive_reads"]
                # 0 = searched for in this mate but not found; NA = not searched for
                for adapter_name, seq in zip(adapter_names, adapters):
                    row[f"adapter_{adapter_name}_reads"] = s["adapter_hits"].get(seq, NA)
            rows.append(row)
    return columns, rows


def write_summary(path, file_stats, collector=None, display_name=None):
    """Write the summary grid to ``path`` as one tab-separated table."""
    columns, rows = build_summary_grid(file_stats, collector, display_name)
    with open(path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
