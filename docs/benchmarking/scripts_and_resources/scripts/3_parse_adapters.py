#!/usr/bin/env python3
"""
Build the adapter ground truth from a (renamed) NEAT golden SAM: split every read
into its genomic part and its adapter part, one TSV row per read pair.

The split uses the soft clip at the read's 3' end in the CIGAR (the end of the
CIGAR, or the start for reverse-strand reads). Reverse-strand reads are turned
back to FASTQ orientation. insert_length (END - START from the read name) is only
reported, not used for the split.

Unmapped, secondary and supplementary records and names without _START_END
coordinates are skipped; a missing mate is written as NA.

Output columns: read_pair_id, insert_length, has_adapter, r1/r2_adapter_start_pos
(0-based = number of genomic bases; 0 = no adapter),
r1/r2_genomic_seq, r1/r2_adapter_seq.

Usage:
  neat_adapter_parse.py --sam renamed_golden.sam -o adapters.tsv
"""

import argparse
import re

CIGAR_REGEX = re.compile(r"(\d+)([MIDNSHP=X])")
COORDINATE_REGEX = re.compile(r"_(\d+)_(\d+)$")
COMPLEMENTATION = str.maketrans("ACGTNacgtn", "TGCANtgcan")
FILE_HEADER = ["read_pair_id", "insert_length", "has_adapter", "r1_adapter_start_pos", "r2_adapter_start_pos", "r1_genomic_seq", "r1_adapter_seq", "r2_genomic_seq", "r2_adapter_seq"]

def revcomp(seq):
    rev_comp = seq.translate(COMPLEMENTATION)[::-1]
    return rev_comp

def analyse_read(flag, cigar, seq):
    reverse = bool(flag & 0x10)
    leading_clip = trailing_clip = 0
    first = True
    for length, operation in CIGAR_REGEX.findall(cigar):
        if operation == "H":
            continue
        if operation == "S":
            if first:
                leading_clip = int(length)
            trailing_clip = int(length)
        else:
            trailing_clip = 0
        first = False
    adapter_length = leading_clip if reverse else trailing_clip
    if reverse:
        seq = revcomp(seq)
    genomic_length = len(seq) - adapter_length
    return seq[:genomic_length], seq[genomic_length:]

def main():
    arguments = argparse.ArgumentParser()
    arguments.add_argument("--sam", required=True, help="SAM file")
    arguments.add_argument("-o", "--out", required=True, help="output TSV")
    arguments_parsed = arguments.parse_args()
    pending_pairs = {}

    with open(arguments_parsed.sam) as file_in, open(arguments_parsed.out, "w") as file_out:
        file_out.write("\t".join(FILE_HEADER) + "\n")

        def write_pair(name, pair):
            r1_genomic, r1_adapter = pair.get(1, ("NA", ""))
            r2_genomic, r2_adapter = pair.get(2, ("NA", ""))
            has_adapter = bool(r1_adapter or r2_adapter)

            def adapter_start(genomic_seq, adapter_seq):
                return len(genomic_seq) if adapter_seq else 0

            file_out.write("\t".join(map(str, [
                "@" + name, pair["insert_length"], "YES" if has_adapter else "NO",
                adapter_start(r1_genomic, r1_adapter), adapter_start(r2_genomic, r2_adapter),
                r1_genomic, r1_adapter, r2_genomic, r2_adapter])) + "\n")

        for line in file_in:
            if line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            read_name, flag, cigar, seq = fields[0], int(fields[1]), fields[5], fields[9]
            if flag & 0x904 or cigar == "*":
                continue

            name = re.sub(r"/[12]$", "", read_name)
            match = COORDINATE_REGEX.search(name)
            if not match:
                continue
            insert_length = int(match.group(2)) - int(match.group(1))
            mate = 2 if flag & 0x80 else 1

            pair = pending_pairs.setdefault(name, {"insert_length": insert_length})
            pair[mate] = analyse_read(flag, cigar, seq)

            if 1 in pair and 2 in pair:
                write_pair(name, pending_pairs.pop(name))

        for name, pair in pending_pairs.items():
            write_pair(name, pair)

if __name__ == "__main__":
    main()