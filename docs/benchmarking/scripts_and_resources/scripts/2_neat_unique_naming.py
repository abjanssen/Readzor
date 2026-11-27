#!/usr/bin/env python3
"""
NEAT (< v4.8.1) produces read headers that contain the contig name, chunk number, and start and stop position.
It is therefore possible for two reads to have the same header, if they come from the same position.
This can cause conflicts in read processing.
 
This script gives every read pair a unique name; in both fastq files, and the golden SAM.
The 4th field of the name (the chunk number) is replaced by an incrementing pair counter.

For example:
    NEAT output:
    NEAT_generated_BHN418_1_0000000014_0000000189
    
    Corrected output:
    NEAT_generated_BHN418_0000123456_0000000014_0000000189

Counters follow the FASTQ order (pair 1 = first record of R1/R2).

Usage:
  neat_rename_reads.py --r1 BHN418_syn_r1.fastq --r2 BHN418_syn_r2.fastq --sam BHN418_syn_golden.sam --out-prefix BHN418_syn_renamed
      
Writes <prefix>_r1.fastq, <prefix>_r2.fastq and <prefix>_golden.sam.

"""
import argparse
import gzip
import re
import sys
from collections import Counter, defaultdict, deque

NAME_REGEX = re.compile(r"^(.*_)(\d+)(_\d+_\d+)$")
COMPLEMENTATION = str.maketrans("ACGTNacgtn", "TGCANtgcan")

def revcomp(seq):
    rev_comp = seq.translate(COMPLEMENTATION)[::-1]
    return rev_comp

def open_close_file(path, mode="rt"):
    if path.endswith(".gz"):
        if "w" in mode:
            return gzip.open(path, mode, compresslevel=4)
        else:
            return gzip.open(path, mode)
    return open(path, mode)

def fastq(path):
    with open_close_file(path) as file:
        while True:
            header = file.readline().rstrip("\n")
            if not header:
                return
            seq = file.readline().rstrip("\n")
            plus = file.readline().rstrip("\n")
            quality = file.readline().rstrip("\n")
            yield header, seq, plus, quality

def split_header(header):
    name = header.lstrip("@").partition(" ")[0]
    suffix = ""
    if name.endswith("/1") or name.endswith("/2"):
        name, suffix = name[:-2], name[-2:]
    return name, suffix

def new_header(old_header, counter, width):
    prefix, _, coordinates = NAME_REGEX.match(old_header).groups()
    new_header = prefix + str(counter).zfill(width) + coordinates
    return new_header

def main():
    arguments = argparse.ArgumentParser()
    arguments.add_argument("--r1", required=True)
    arguments.add_argument("--r2", required=True)
    arguments.add_argument("--sam", required=True)
    arguments.add_argument("--out-prefix", required=True)
    arguments.add_argument("--width", type=int, default=10, help="zero-padding of the counter (default 10)")
    arguments.add_argument("--gzip", action="store_true", help="gzip the FASTQ outputs")
    arguments_parsed = arguments.parse_args()

    counts = Counter()
    for header, _, _, _ in fastq(arguments_parsed.r1):
        counts[split_header(header)[0]] += 1
    dup_names = set()
    for name, count in counts.items():
        if count > 1:
            dup_names.add(name)

    if arguments_parsed.gzip:
        ext = ".fastq.gz"
    else:
        ext = ".fastq"

    unique_map = {} 
    dup_map = defaultdict(deque)

    prefix = arguments_parsed.out_prefix
    with open_close_file(f"{prefix}_r1{ext}", "wt") as output_1, open_close_file(f"{prefix}_r2{ext}", "wt") as output_2:
        i = 0
        for (header1, seq1, plus1, qual1), (header2, seq2, plus2, qual2) in zip(fastq(arguments_parsed.r1), fastq(arguments_parsed.r2), strict=True):
            i += 1
            name1, suffix1 = split_header(header1)
            name2, suffix2 = split_header(header2)
            if name1 != name2:
                sys.exit(f"ERROR: R1/R2 out of sync at pair {i}: {name1} vs {name2}")
            new_name = new_header(name1, i, arguments_parsed.width)
            if name1 in dup_names:
                dup_map[(name1, 1, seq1)].append(new_name)
                dup_map[(name1, 2, seq2)].append(new_name)
            else:
                unique_map[name1] = new_name
            output_1.write(f"@{new_name}{suffix1}\n{seq1}\n{plus1}\n{qual1}\n")
            output_2.write(f"@{new_name}{suffix2}\n{seq2}\n{plus2}\n{qual2}\n")

    with open_close_file(arguments_parsed.sam) as file_in, open(f"{arguments_parsed.out_prefix}_golden.sam", "w") as file_out:
        for line in file_in:
            if line.startswith("@"):
                file_out.write(line)
                continue
            fields = line.rstrip("\n").split("\t")
            name, suffix = split_header(fields[0])
            flag = int(fields[1])
            if name in unique_map:
                new_name = unique_map[name]
            else:
                mate = 2 if flag & 0x80 else 1
                seq = revcomp(fields[9]) if flag & 0x10 else fields[9]
                candidates = dup_map.get((name, mate, seq)) or dup_map.get((name, mate, fields[9]))
                if not candidates:
                    sys.exit(f"ERROR: no FASTQ match for SAM record {fields[0]} (flag {flag})")
                new_name = candidates.popleft() if not flag & 0x900 else candidates[0]
            fields[0] = new_name + suffix
            file_out.write("\t".join(fields) + "\n")

if __name__ == "__main__":
    main()