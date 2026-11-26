#!/usr/bin/env python3
"""
Benchmark Readzor's gzip settings.

Every command in COMMANDS is run on both input types (.fastq.gz and plain
.fastq), REPEATS times per thread count. Per run, wall-clock time, peak memory
and output size are appended to RESULTS_CSV.

  readzor_nogzip_*read  no --gzip (plain output); --gzip-method only selects how
                        .fastq.gz input is unpacked, so the difference between
                        the two is the cost of the decompression backend.
                        (On plain .fastq input --gzip-method does nothing.)
  readzor_isal0-3       --gzip output with isal, levels 0-3
  readzor_gzip1-9       --gzip output with stdlib gzip, levels 1-9

Peak memory is the RSS summed over the whole process tree (Readzor plus all its
worker processes), polled every SAMPLE_INTERVAL seconds. `/usr/bin/time -v`
only sees one process and under-reports multiprocessing tools. Summing RSS can
slightly overcount shared copy-on-write pages, but is fine for comparing trends.

The plain .fastq inputs are made once from the .gz files (before any timing)
and kept. Files and folders created during a run (in the working directory and
the input directories) are deleted after each repeat.

Usage:
    python3 benchmark_gziplevels.py                   # THREAD_COUNTS, both input types
    python3 benchmark_gziplevels.py 1 2 4 8           # these thread counts
    python3 benchmark_gziplevels.py 8 -r 1 -i gz      # 1 repeat, only .fastq.gz input

Requires: pip install psutil
"""
import argparse
import csv
import gzip
import os
import shutil
import subprocess
import time

import psutil

# ============================================================
# CONFIGURATION
# ============================================================
THREAD_COUNTS = [4]
REPEATS = 5
INPUT_TYPES = ["gz", "fastq"]
SAMPLE_INTERVAL = 0.1
RESULTS_CSV = "./results/gzip_bench/benchmark_results_gziplevels.csv"

R1 = "data/fastq/BHN418_syn_r1.fastq.gz"  # always the .gz paths; plain .fastq is derived
R2 = "data/fastq/BHN418_syn_r2.fastq.gz"
READZOR_BASE = "python ./src/readzor_isal_gzip_comparison.py --input-paired {r1} {r2} --threads {threads}"

COMMANDS = {
    "readzor_nogzip_isalread": "--gzip-method isal",
    "readzor_nogzip_gzipread": "--gzip-method gzip",
    "readzor_isal0": "--gzip --gzip-method isal --gzip-level 0",
    "readzor_isal1": "--gzip --gzip-method isal --gzip-level 1",
    "readzor_isal2": "--gzip --gzip-method isal --gzip-level 2",
    "readzor_isal3": "--gzip --gzip-method isal --gzip-level 3",
    "readzor_gzip1": "--gzip --gzip-method gzip --gzip-level 1",
    "readzor_gzip2": "--gzip --gzip-method gzip --gzip-level 2",
    "readzor_gzip3": "--gzip --gzip-method gzip --gzip-level 3",
    "readzor_gzip4": "--gzip --gzip-method gzip --gzip-level 4",
    "readzor_gzip5": "--gzip --gzip-method gzip --gzip-level 5",
    "readzor_gzip6": "--gzip --gzip-method gzip --gzip-level 6",
    "readzor_gzip7": "--gzip --gzip-method gzip --gzip-level 7",
    "readzor_gzip8": "--gzip --gzip-method gzip --gzip-level 8",
    "readzor_gzip9": "--gzip --gzip-method gzip --gzip-level 9",
}

CSV_FIELDS = ["package", "experiment", "input_type", "threads", "repeat", "wall_time_s", "peak_rss_mb", "output_size_mb", "exit_code", "command"]

def plain_fastq(gz_path):
    plain_path = gz_path[:-3]
    if not os.path.exists(plain_path):
        print(f"Decompressing {gz_path} -> {plain_path} (one-time) ...", flush=True)
        with gzip.open(gz_path, "rb") as source, open(plain_path + ".part", "wb") as target:
            shutil.copyfileobj(source, target, length=1024 * 1024)
        os.replace(plain_path + ".part", plain_path)
    return plain_path

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

def run_once(package, input_type, command, threads, repeat, repeats):
    print(f">>> {package} | input={input_type} | threads={threads} | repeat={repeat}/{repeats}", flush=True)
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

    return {"package": package, "experiment": "all", "input_type": input_type,
            "threads": threads, "repeat": repeat,
            "wall_time_s": f"{elapsed:.3f}", "peak_rss_mb": f"{peak / 1024 ** 2:.2f}",
            "output_size_mb": f"{output_bytes / 1024 ** 2:.2f}",
            "exit_code": process.returncode, "command": command}


def main():
    parser = argparse.ArgumentParser(description="Benchmark Readzor gzip settings across thread counts and input types.")
    parser.add_argument("threads", nargs="*", type=int, help=f"thread count(s) (default: {THREAD_COUNTS})")
    parser.add_argument("-r", "--repeats", type=int, default=REPEATS, help=f"repeats (default: {REPEATS})")
    parser.add_argument("-i", "--input-types", nargs="+", choices=["gz", "fastq"], default=INPUT_TYPES,
                        help=f"input file type(s) (default: {INPUT_TYPES})")
    arguments = parser.parse_args()

    # make the plain inputs before any timing starts
    inputs = {}
    for input_type in arguments.input_types:
        inputs[input_type] = (R1, R2) if input_type == "gz" else (plain_fastq(R1), plain_fastq(R2))

    os.makedirs(os.path.dirname(RESULTS_CSV), exist_ok=True)
    write_header = not os.path.exists(RESULTS_CSV) or os.path.getsize(RESULTS_CSV) == 0

    for package, flags in COMMANDS.items():
        for input_type in arguments.input_types:
            r1, r2 = inputs[input_type]
            for threads in arguments.threads or THREAD_COUNTS:
                command = READZOR_BASE.format(r1=r1, r2=r2, threads=threads) + " " + flags
                for repeat in range(1, arguments.repeats + 1):
                    row = run_once(package, input_type, command, threads, repeat, arguments.repeats)
                    with open(RESULTS_CSV, "a", newline="") as file:
                        writer = csv.DictWriter(file, fieldnames=CSV_FIELDS)
                        if write_header:
                            writer.writeheader()
                            write_header = False
                        writer.writerow(row)

    print(f"\nDone. Results written to {RESULTS_CSV}")

if __name__ == "__main__":
    main()