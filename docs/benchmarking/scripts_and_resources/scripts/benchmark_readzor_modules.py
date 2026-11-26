#!/usr/bin/env python3
"""
Benchmark the cost of each Readzor filter module.

Every filter in FILTER_FLAGS is switched on alone on top of the same base
command (gz -> gz), plus a baseline without filters and "all filters" runs
(one per adapter mismatch setting). Each run is repeated REPEATS times per
thread count. Per run, wall-clock time, peak memory and output size are recorded.

Requires: pip install psutil
"""
import argparse
import csv
import os
import shutil
import subprocess
import time

import psutil

# ============================================================
# CONFIGURATION
# ============================================================
THREAD_COUNTS = [4]
REPEATS = 3
SAMPLE_INTERVAL = 0.1
RESULTS_CSV = "./results.csv"

R1 = "data/r1.fastq.gz"
R2 = "data/r2.fastq.gz"
GZIP_LEVEL = 1
READZOR_BASE = f"readzor --input-paired {R1} {R2} --gzip --gzip-level {GZIP_LEVEL} --threads {{threads}}"

FILTER_FLAGS = {
    "cut":         "--cut-flag --cut-both 5",
    "endqual":     "--endqual-filter-flag",
    "nfilter":     "--n-filter",
    "slider":      "--slider-filter-flag",
    "poly":        "--poly-filter-flag",
    "adapter_mm0": "--adapter-filter-flag --adapter-mismatch 0",
    "adapter_mm1": "--adapter-filter-flag --adapter-mismatch 1",
    "adapter_mm2": "--adapter-filter-flag --adapter-mismatch 2",
    "adapter_mm3": "--adapter-filter-flag --adapter-mismatch 3",
    "adapter_mm4": "--adapter-filter-flag --adapter-mismatch 4",
    "adapter_mm5": "--adapter-filter-flag --adapter-mismatch 5",
    "overlap":     "--overlap-filter-flag",
    "kmer":        "--kmer-filter-flag",
}

RUN_BASELINE = True
RUN_ALL_FILTERS = True

CSV_FIELDS = ["package", "filter", "experiment", "input_type", "threads", "repeat", "wall_time_s", "peak_rss_mb", "output_size_mb", "exit_code", "command"]

def build_runs(selected):
    runs = {}
    if RUN_BASELINE:
        runs["nofilter"] = ""
    for name in selected:
        runs[name] = FILTER_FLAGS[name]
    if RUN_ALL_FILTERS:
        others = [flags for name, flags in FILTER_FLAGS.items() if not name.startswith("adapter")]
        for name, adapter_flags in FILTER_FLAGS.items():
            if name.startswith("adapter"):
                runs["allfilter" + name.removeprefix("adapter")] = " ".join(others + [adapter_flags])
    return runs

def snapshot():
    paths = set()
    for directory in {os.path.abspath(path) for path in (".", os.path.dirname(R1) or ".", os.path.dirname(R2) or ".")}:
        for name in os.listdir(directory):
            paths.add(os.path.join(directory, name))
    return paths

def fastq_size_and_cleanup(before):
    total = 0
    for path in snapshot() - before:
        if path == os.path.abspath(RESULTS_CSV):
            continue
        if os.path.isdir(path):
            for root, _, files in os.walk(path):
                for name in files:
                    if name.endswith((".fastq", ".fastq.gz")):
                        total += os.path.getsize(os.path.join(root, name))
            shutil.rmtree(path)
        else:
            if path.endswith((".fastq", ".fastq.gz")):
                total += os.path.getsize(path)
            os.remove(path)
    return total

def tree_rss(pid):
    try:
        processes = [psutil.Process(pid)]
        processes += processes[0].children(recursive=True)
    except psutil.NoSuchProcess:
        return 0
    total = 0
    for process in processes:
        try:
            total += process.memory_info().rss
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass
    return total

def run_once(filter_name, command, threads, repeat, repeats):
    print(f">>> {filter_name} | threads={threads} | repeat={repeat}/{repeats}", flush=True)
    before = snapshot()
    peak = 0
    start = time.monotonic()
    process = subprocess.Popen(command, shell=True)
    try:
        while process.poll() is None:
            peak = max(peak, tree_rss(process.pid))
            time.sleep(SAMPLE_INTERVAL)
        elapsed = time.monotonic() - start
    finally:
        if process.poll() is None:
            process.terminate()
        output_bytes = fastq_size_and_cleanup(before)

    return {"package": f"readzor_{filter_name}", "filter": filter_name, "experiment": "all",
            "input_type": "gz", "threads": threads, "repeat": repeat,
            "wall_time_s": f"{elapsed:.3f}", "peak_rss_mb": f"{peak / 1024 ** 2:.2f}",
            "output_size_mb": f"{output_bytes / 1024 ** 2:.2f}",
            "exit_code": process.returncode, "command": command}

def main():
    parser = argparse.ArgumentParser(description="Benchmark Readzor filter flags across thread counts.")
    parser.add_argument("threads", nargs="*", type=int, help=f"thread count(s) (default: {THREAD_COUNTS})")
    parser.add_argument("-r", "--repeats", type=int, default=REPEATS, help=f"repeats (default: {REPEATS})")
    parser.add_argument("-f", "--flags", nargs="+", choices=list(FILTER_FLAGS), default=list(FILTER_FLAGS),
                        help="only test these filter flags (default: all)")
    arguments = parser.parse_args()

    os.makedirs(os.path.dirname(RESULTS_CSV), exist_ok=True)
    write_header = not os.path.exists(RESULTS_CSV) or os.path.getsize(RESULTS_CSV) == 0

    for filter_name, flags in build_runs(arguments.flags).items():
        for threads in arguments.threads or THREAD_COUNTS:
            command = READZOR_BASE.format(threads=threads) + (" " + flags if flags else "")
            for repeat in range(1, arguments.repeats + 1):
                row = run_once(filter_name, command, threads, repeat, arguments.repeats)
                with open(RESULTS_CSV, "a", newline="") as file:
                    writer = csv.DictWriter(file, fieldnames=CSV_FIELDS)
                    if write_header:
                        writer.writeheader()
                        write_header = False
                    writer.writerow(row)

    print(f"\nDone. Results written to {RESULTS_CSV}")


if __name__ == "__main__":
    main()