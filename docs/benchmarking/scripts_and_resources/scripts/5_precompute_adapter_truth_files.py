#!/usr/bin/env python3
"""
Precompute the truth for judging (run once per truth TSV).

Reads the truth TSV from the adapter parser and stores the read names and the
genomic / adapter length of every read as .npy files in OUT_DIR:

  names.npy              read names, sorted (without the leading '@')
  r1_genomic.npy  r1_adapter.npy     genomic / adapter length of R1 (same order as names)
  r2_genomic.npy  r2_adapter.npy     same for R2

Usage:
  precompute_truth.py parsed_adapters.tsv truth_pre/
"""
import argparse
import csv
import os

import numpy as np

COLUMNS = {"r1_genomic_seq": "r1_genomic", "r1_adapter_seq": "r1_adapter", "r2_genomic_seq": "r2_genomic", "r2_adapter_seq": "r2_adapter"}

def main():
    arguments = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    arguments.add_argument("truth", help="truth TSV from the adapter parser")
    arguments.add_argument("outdir", help="output directory (created)")
    arguments.add_argument("--chunk", type=int, default=10_000, help="rows converted to numpy at a time")
    arguments_parsed = arguments.parse_args()

    buffer_names = []
    buffer_lengths = {column: [] for column in COLUMNS}
    name_chunks = []
    length_chunks = {column: [] for column in COLUMNS}

    def flush_buffer():
        name_chunks.append(np.array(buffer_names, dtype="S"))
        buffer_names.clear()
        for column in COLUMNS:
            length_chunks[column].append(np.array(buffer_lengths[column], dtype=np.uint32))
            buffer_lengths[column].clear()

    with open(arguments_parsed.truth, newline="") as file:
        for row in csv.DictReader(file, delimiter="\t"):
            buffer_names.append(row["read_pair_id"].lstrip("@"))
            for column in COLUMNS:
                buffer_lengths[column].append(len(row[column]))
            if len(buffer_names) >= arguments_parsed.chunk:
                flush_buffer()
    flush_buffer()

    names = np.concatenate(name_chunks)
    sort_order = np.argsort(names)

    os.makedirs(arguments_parsed.outdir, exist_ok=True)
    np.save(os.path.join(arguments_parsed.outdir, "names.npy"), names[sort_order])
    for column, file_name in COLUMNS.items():
        column_lengths = np.concatenate(length_chunks[column])[sort_order]
        np.save(os.path.join(arguments_parsed.outdir, file_name + ".npy"), column_lengths)

if __name__ == "__main__":
    main()