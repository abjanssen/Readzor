#!/usr/bin/env python3
"""
Benchmark fastp, Trimmomatic, Readzor and cutadapt on single-end data.

Every timed run covers ALL files in INPUT_DIR:
  fastp        ./src/fastp_parallel_fastq.py over the input directory
  readzor      one call, all files via --input-unpaired
  trimmomatic, cutadapt   a shell loop over the files

Per run, wall-clock time and peak memory are appended to the results CSV.
Peak memory is the RSS summed over the whole process tree, polled every
SAMPLE_INTERVAL seconds (`/usr/bin/time -v` only sees one process).

Usage:
    python3 benchmark_tools_3_se.py 1 2 4 8                    # these thread counts
    python3 benchmark_tools_3_se.py 8 -r 1                     # 1 repeat (quick test)
    python3 benchmark_tools_3_se.py 8 -o results/bench_t{threads}.csv
"""
import argparse
import csv
import glob
import os
import shutil
import subprocess
import time

import psutil

THREAD_COUNTS = [4]
REPEATS = 3
SAMPLE_INTERVAL = 0.1
RESULTS_CSV = "benchmark_results_se.csv"

INPUT_DIR = "./data/SE"

COMMANDS = {
    "fastp": (
        "mkdir -p {out}_dir && "
        "python ./src/fastp_parallel_fastq.py -i {input_dir}/ -o {out}_dir/ "
        "-p 1 -c ./src/fastp_wrapper.sh "
        '-a "--thread {threads} --compression 5 '
        "--cut_front --cut_front_window_size 1 --cut_front_mean_quality 20 "
        "--cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 20 "
        "--adapter_sequence AGATCGGAAGAGCACACGTCTGAACTCCAGTCA "
        "--json /dev/null --html /dev/null "
        "--length_required 50 --n_base_limit 50 --disable_quality_filtering "
        '--disable_trim_poly_g --dont_eval_duplication"'
    ),
    "trimmomatic": (
        "for input_file in {input_dir}/*; do "
        "micromamba run -p ./envs/mm_trimmomatic trimmomatic -Xmx10G SE "
        "-phred33 -threads {threads} -compressLevel 5 "
        '"$input_file" "{out}_$(basename "$input_file")" '
        "ILLUMINACLIP:./adapters_truseq_fwd.fna:2:30:10:8:true "
        "LEADING:20 TRAILING:20 MINLEN:50; "
        "done"
    ),
    "cutadapt": (
        "for input_file in {input_dir}/*; do "
        "micromamba run -p ./envs/mm_cutadapt cutadapt "
        "-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA --no-indels -e 0.1 --overlap 3 "
        "--cores {threads} --compression-level 5 "
        "--quality-cutoff 20,20 --minimum-length 50 "
        '-o "{out}_$(basename "$input_file")" "$input_file"; '
        "done"
    ),
    "readzor": (
        "micromamba run -p ./envs/mm_readzor readzor --input-unpaired {input_dir}/* "
        "--adapter-filter-flag --adapter-read-1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCA "
        "--adapter-mismatch 2 --adapter-seed 6 "
        "--endqual-filter-flag --endqual-min-both 20 "
        "--min-length-output 50 "
        "--gzip --gzip-level 1 --threads {threads} "
        "--output {out}"
    ),
}

CSV_FIELDS = ["package", "experiment", "threads", "repeat", "wall_time_s", "peak_rss_mb", "exit_code", "command"]


def tree_rss(pid):
    """RSS (bytes) of a process plus all its descendants."""
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


def run_once(package, command, threads, repeat, repeats):
    print(f">>> {package} | threads={threads} | repeat={repeat}/{repeats}", flush=True)
    peak = 0
    start = time.monotonic()
    process = subprocess.Popen(command, shell=True)
    try:
        while process.poll() is None:
            peak = max(peak, tree_rss(process.pid))
            time.sleep(SAMPLE_INTERVAL)
    except KeyboardInterrupt:
        process.terminate()
        raise
    return {"package": package, "experiment": "all", "threads": threads, "repeat": repeat,
            "wall_time_s": f"{time.monotonic() - start:.3f}", "peak_rss_mb": f"{peak / 1024 ** 2:.2f}",
            "exit_code": process.returncode, "command": command}


def remove_output(out_prefix):
    """Remove everything matching '<out_prefix>*'."""
    for path in glob.glob(glob.escape(out_prefix) + "*"):
        if os.path.isdir(path):
            shutil.rmtree(path)
        else:
            os.remove(path)


def main():
    parser = argparse.ArgumentParser(description="Benchmark trimmers on single-end data across thread counts.")
    parser.add_argument("threads", nargs="*", type=int, help=f"thread count(s) (default: {THREAD_COUNTS})")
    parser.add_argument("-r", "--repeats", type=int, default=REPEATS, help=f"repeats (default: {REPEATS})")
    parser.add_argument("-o", "--output", default=RESULTS_CSV,
                        help=f"results CSV (default: {RESULTS_CSV}); {{threads}} is replaced by the thread count(s)")
    parser.add_argument("--keep-output", action="store_true", help="keep each run's output (default: delete)")
    arguments = parser.parse_args()

    thread_counts = arguments.threads or THREAD_COUNTS
    results_csv = arguments.output.format(threads="-".join(str(threads) for threads in thread_counts))
    if os.path.dirname(results_csv):
        os.makedirs(os.path.dirname(results_csv), exist_ok=True)
    write_header = not os.path.exists(results_csv) or os.path.getsize(results_csv) == 0
    print(f"Threads: {thread_counts} | repeats: {arguments.repeats} | results: {results_csv}", flush=True)

    for package, template in COMMANDS.items():
        for threads in thread_counts:
            out_prefix = f"out_{package}_all_t{threads}"
            command = template.format(threads=threads, out=out_prefix, input_dir=INPUT_DIR)
            for repeat in range(1, arguments.repeats + 1):
                row = run_once(package, command, threads, repeat, arguments.repeats)
                with open(results_csv, "a", newline="") as file:
                    writer = csv.DictWriter(file, fieldnames=CSV_FIELDS)
                    if write_header:
                        writer.writeheader()
                        write_header = False
                    writer.writerow(row)
                if not arguments.keep_output:
                    remove_output(out_prefix)

    print(f"\nDone. Results written to {results_csv}")


if __name__ == "__main__":
    main()