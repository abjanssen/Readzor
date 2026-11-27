#!/usr/bin/env python3
"""
Benchmark fastp, Trimmomatic, Readzor and cutadapt on paired-end data.

Every timed run covers ALL read pairs in INPUT_DIR:
  fastp        ./src/fastp_parallel_fastq.py over the input directory
  readzor      one call, all pairs via --input-paired (r1 r2 r1 r2 ...)
  trimmomatic, cutadapt   a shell loop over the pairs

Per run, wall-clock time and peak memory are appended to the results CSV.
Peak memory is the RSS summed over the whole process tree, polled every
SAMPLE_INTERVAL seconds (`/usr/bin/time -v` only sees one process).

Usage:
    python3 benchmark_tools_3_pe.py 1 2 4 8                    # these thread counts
    python3 benchmark_tools_3_pe.py 8 -r 1                     # 1 repeat (quick test)
    python3 benchmark_tools_3_pe.py 8 -o results/bench_t{threads}.csv
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
REPEATS = 1
SAMPLE_INTERVAL = 0.1
RESULTS_CSV = "benchmark_results.csv"

INPUT_DIR = "./data/fastq/PE/wgs_pneumo/fastq"
READ_ONE_SUFFIX = "_R1.fastq.gz"
READ_TWO_SUFFIX = "_R2.fastq.gz"
ADAPTER_FILE = "./adapters_truseq.fna"
ADAPTER_FILE_OVERLAP_ADAPTER_TRIMMOMATIC = "./adapters_truseq_palindrome_simple_trimmomatic"
ADAPTER_FILE_OVERLAP_TRIMMOMATIC = "./adapters_truseq_palindrome_trimmomatic.fna"

# Shell loop over all pairs in INPUT_DIR; sets $read_one, $read_two, $sample_name.
# Shell braces are doubled because the commands go through str.format().
LOOP_HEAD = (
    "for read_one in {input_dir}/*{read_one_suffix}; do "
    'sample_prefix="${{read_one%{read_one_suffix}}}"; '
    'read_two="${{sample_prefix}}{read_two_suffix}"; '
    'sample_name=$(basename "$sample_prefix"); '
)

# name -> command. Placeholders: {threads} {out} {input_dir} {experiments} and the adapter files.
COMMANDS = {
    "fastp_overlap_adapter": (
        "mkdir -p {out}_dir && "
        "python ./src/fastp_parallel_fastq.py -i {input_dir}/ -o {out}_dir/ "
        "-p 1 -c ./src/fastp_wrapper.sh "
        '-a "--thread {threads} --compression 5 '
        "--json /dev/null --html /dev/null "
        "--cut_front --cut_front_window_size 1 --cut_front_mean_quality 20 "
        "--cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 20 "
        "--overlap_len_require 30 --overlap_diff_limit 5 --overlap_diff_percent_limit 20 "
        "--adapter_fasta {adapter_file} "
        "--length_required 50 --n_base_limit 50 --disable_quality_filtering "
        '--disable_trim_poly_g --dont_eval_duplication"'
    ),
    "fastp_overlap": (
        "mkdir -p {out}_dir && "
        "python ./src/fastp_parallel_fastq.py -i {input_dir}/ -o {out}_dir/ "
        "-p 1 -c ./src/fastp_wrapper.sh "
        '-a "--thread {threads} --compression 5 '
        "--json /dev/null --html /dev/null "
        "--cut_front --cut_front_window_size 1 --cut_front_mean_quality 20 "
        "--cut_tail --cut_tail_window_size 1 --cut_tail_mean_quality 20 "
        "--overlap_len_require 30 --overlap_diff_limit 5 --overlap_diff_percent_limit 20 "
        "--length_required 50 --n_base_limit 50 --disable_quality_filtering "
        '--disable_trim_poly_g --dont_eval_duplication"'
    ),
    "trimmomatic_overlap_adapter": (
        LOOP_HEAD
        + "micromamba run -p ./envs/mm_trimmomatic trimmomatic -Xmx10G PE "
        "-phred33 -threads {threads} -compressLevel 5 "
        '"$read_one" "$read_two" '
        '"{out}_${{sample_name}}_R1.fastq.gz" "{out}_${{sample_name}}_R1_unpaired.fastq.gz" '
        '"{out}_${{sample_name}}_R2.fastq.gz" "{out}_${{sample_name}}_R2_unpaired.fastq.gz" '
        "ILLUMINACLIP:{adapter_file_overlap_adapter_trimmomatic}:2:30:10:8:true "
        "LEADING:20 TRAILING:20 MINLEN:50; "
        "done"
    ),
    "trimmomatic_overlap": (
        LOOP_HEAD
        + "micromamba run -p ./envs/mm_trimmomatic trimmomatic -Xmx10G PE "
        "-phred33 -threads {threads} -compressLevel 5 "
        '"$read_one" "$read_two" '
        '"{out}_${{sample_name}}_R1.fastq.gz" "{out}_${{sample_name}}_R1_unpaired.fastq.gz" '
        '"{out}_${{sample_name}}_R2.fastq.gz" "{out}_${{sample_name}}_R2_unpaired.fastq.gz" '
        "ILLUMINACLIP:{adapter_file_overlap_trimmomatic}:2:30:10:8:true "
        "LEADING:20 TRAILING:20 MINLEN:50; "
        "done"
    ),
    "trimmomatic_adapter": (
        LOOP_HEAD
        + "micromamba run -p ./envs/mm_trimmomatic trimmomatic -Xmx10G PE "
        "-phred33 -threads {threads} -compressLevel 5 "
        '"$read_one" "$read_two" '
        '"{out}_${{sample_name}}_R1.fastq.gz" "{out}_${{sample_name}}_R1_unpaired.fastq.gz" '
        '"{out}_${{sample_name}}_R2.fastq.gz" "{out}_${{sample_name}}_R2_unpaired.fastq.gz" '
        "ILLUMINACLIP:{adapter_file}:2:30:10:8:true "
        "LEADING:20 TRAILING:20 MINLEN:50; "
        "done"
    ),
    "readzor_adapter_mm0_overlap": (
        "micromamba run -p ./envs/mm_readzor readzor "
        "--input-paired {experiments} "
        "--adapter-filter-flag --adapter-excl {adapter_file} --adapter-mismatch 0 --adapter-seed 6 "
        "--overlap-filter-flag --overlap-min-length 50 --overlap-mismatches 5 --overlap-portion-mismatch 4 "
        "--endqual-filter-flag --endqual-min-both 20 "
        "--min-length-output 50 "
        "--gzip --gzip-level 1 --threads {threads} "
        "--output {out}"
    ),
    "readzor_adapter_mm2_overlap": (
        "micromamba run -p ./envs/mm_readzor readzor "
        "--input-paired {experiments} "
        "--adapter-filter-flag --adapter-excl {adapter_file} --adapter-mismatch 2 --adapter-seed 6 "
        "--overlap-filter-flag --overlap-min-length 50 --overlap-mismatches 5 --overlap-portion-mismatch 4 "
        "--endqual-filter-flag --endqual-min-both 20 "
        "--min-length-output 50 "
        "--gzip --gzip-level 1 --threads {threads} "
        "--output {out}"
    ),
    "readzor_overlap": (
        "micromamba run -p ./envs/mm_readzor readzor "
        "--input-paired {experiments} "
        "--overlap-filter-flag --overlap-min-length 50 --overlap-mismatches 5 --overlap-portion-mismatch 4 "
        "--endqual-filter-flag --endqual-min-both 20 "
        "--min-length-output 50 "
        "--gzip --gzip-level 1 --threads {threads} "
        "--output {out}"
    ),
    "readzor_adapter_mm2": (
        "micromamba run -p ./envs/mm_readzor readzor "
        "--input-paired {experiments} "
        "--adapter-filter-flag --adapter-excl {adapter_file} --adapter-mismatch 2 --adapter-seed 6 "
        "--endqual-filter-flag --endqual-min-both 20 "
        "--min-length-output 50 "
        "--gzip --gzip-level 1 --threads {threads} "
        "--output {out}"
    ),
    "readzor_adapter_mm0": (
        "micromamba run -p ./envs/mm_readzor readzor "
        "--input-paired {experiments} "
        "--adapter-filter-flag --adapter-excl {adapter_file} --adapter-mismatch 0 --adapter-seed 6 "
        "--endqual-filter-flag --endqual-min-both 20 "
        "--min-length-output 50 "
        "--gzip --gzip-level 1 --threads {threads} "
        "--output {out}"
    ),
    "cutadapt_adapter": (
        LOOP_HEAD
        + "micromamba run -p ./envs/mm_cutadapt cutadapt "
        "-a file:{adapter_file} -A file:{adapter_file} --no-indels -e 0.1 --overlap 3 "
        "--cores {threads} --compression-level 5 "
        "--quality-cutoff 20,20 --minimum-length 50 "
        '-o "{out}_${{sample_name}}_R1.fastq.gz" -p "{out}_${{sample_name}}_R2.fastq.gz" '
        '"$read_one" "$read_two"; '
        "done"
    ),
}

CSV_FIELDS = ["package", "experiment", "threads", "repeat", "wall_time_s", "peak_rss_mb", "exit_code", "command"]


def all_pairs():
    """All pairs in INPUT_DIR as 'r1 r2 r1 r2 ...' (Readzor --input-paired)."""
    paths = []
    for read_one in sorted(glob.glob(os.path.join(INPUT_DIR, "*" + READ_ONE_SUFFIX))):
        paths += [read_one, read_one[:-len(READ_ONE_SUFFIX)] + READ_TWO_SUFFIX]
    return " ".join(paths)


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
    parser = argparse.ArgumentParser(description="Benchmark trimmers on paired-end data across thread counts.")
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

    experiments = all_pairs()
    for package, template in COMMANDS.items():
        for threads in thread_counts:
            out_prefix = f"out_{package}_all_t{threads}"
            command = template.format(
                threads=threads, out=out_prefix, experiments=experiments, input_dir=INPUT_DIR,
                read_one_suffix=READ_ONE_SUFFIX, read_two_suffix=READ_TWO_SUFFIX,
                adapter_file=ADAPTER_FILE,
                adapter_file_overlap_adapter_trimmomatic=ADAPTER_FILE_OVERLAP_ADAPTER_TRIMMOMATIC,
                adapter_file_overlap_trimmomatic=ADAPTER_FILE_OVERLAP_TRIMMOMATIC)
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