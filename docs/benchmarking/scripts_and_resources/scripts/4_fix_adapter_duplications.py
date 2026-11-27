#!/usr/bin/env python3
"""
Replace NEAT's repeated adapter copies with a realistic TruSeq read-through.

When the insert is much shorter than the read, NEAT fills the rest of the read
with the adapter pasted again and again. A real read instead continues through
the full adapter, the index, the flow-cell (P7/P5) sequence and then into
poly-G (2-colour: NovaSeq/NextSeq/NovaSeq X) or random low-quality bases
(4-colour: MiSeq/HiSeq).

What the script does, per read with an adapter (from the truth TSV):
  * genomic part             -> unchanged
  * first NEAT adapter copy  -> unchanged (keeps NEAT's sequencing errors)
  * everything after it      -> replaced by the realistic continuation
                                (+ substitutions at --error-rate)
Read lengths, quality strings and adapter start positions do not change, so the
SAM CIGARs stay valid. The truth TSV is rewritten with the new adapter_seq.

The split is purely length-based: the adapter starts where the truth TSV says,
the first len(--adapter-rN) bases after that are NEAT's adapter copy (kept), and
everything beyond is replaced. No sequence matching is done.

NEAT adapters (defaults):   R1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCA   (33 bp)
                            R2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT   (33 bp)
Continuation after them:    R1 C <i7> ATCTCGTATGCCGTCTTCTGCTTG  then tail
                            R2   <i5> GTGTAGATCTCGGTGGTCGCCGTATCATT then tail
(R2 contains the i5 as sequenced in forward-strand workflows; set --i5 for yours.)

Usage:
  neat_fix_adapter_readthrough.py --r1 X_r1.fastq --r2 X_r2.fastq --truth parsed_adapters.tsv \
      --out-prefix X_fixed [--tail polyG|random] [--i7 ATCACG] [--i5 TATAGCCT]
Writes <prefix>_r1.fastq, <prefix>_r2.fastq, <prefix>_parsed_adapters.tsv
"""
import argparse
import csv
import gzip
import random

R1_READ_THROUGH = "C{index}ATCTCGTATGCCGTCTTCTGCTTG"
R2_READ_THROUGH = "{index}GTGTAGATCTCGGTGGTCGCCGTATCATT"

def open_file(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)

def fastq(path):
    with open_file(path) as file:
        while True:
            header = file.readline()
            if not header:
                return
            seq, plus, quality = file.readline(), file.readline(), file.readline()
            yield header.rstrip("\n"), seq.rstrip("\n"), plus.rstrip("\n"), quality.rstrip("\n")

def read_name(header):
    name = header.split()[0].lstrip("@")
    return name[:-2] if name[-2:] in ("/1", "/2") else name

def mutate(seq, error_rate, random_generator):
    if error_rate <= 0:
        return seq
    bases = list(seq)
    for position, base in enumerate(bases):
        if random_generator.random() < error_rate:
            bases[position] = random_generator.choice([other for other in "ACGT" if other != base])
    return "".join(bases)

def continuation(read_through, length, tail, error_rate, random_generator):
    bases = read_through[:length]
    missing = length - len(bases)
    if missing > 0:
        if tail == "polyG":
            bases += "G" * missing
        else:
            bases += "".join(random_generator.choice("ACGT") for _ in range(missing))
    return mutate(bases, error_rate, random_generator)

def main():
    arguments = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    arguments.add_argument("--r1", required=True)
    arguments.add_argument("--r2", required=True)
    arguments.add_argument("--truth", required=True, help="parsed adapter file belonging to these FASTQs")
    arguments.add_argument("--out-prefix", required=True)
    arguments.add_argument("--i7", default="ATCACG", help="i7 index in R1 read-through (default ATCACG)")
    arguments.add_argument("--i5", default="TATAGCCT", help="i5 index in R2 read-through (default TATAGCCT)")
    arguments.add_argument("--tail", choices=["polyG", "random"], default="random", help="what follows the flow-cell sequence")
    arguments.add_argument("--error-rate", type=float, default=0.002, help="substitution rate on the replaced bases (default 0.002)")
    arguments.add_argument("--adapter-r1", default="AGATCGGAAGAGCACACGTCTGAACTCCAGTCA", help="NEAT's R1 adapter (only its length is used)")
    arguments.add_argument("--adapter-r2", default="AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT", help="NEAT's R2 adapter (only its length is used)")
    arguments.add_argument("--seed", type=int, default=1)
    arguments_parsed = arguments.parse_args()

    read_through = {1: R1_READ_THROUGH.format(index=arguments_parsed.i7),
                    2: R2_READ_THROUGH.format(index=arguments_parsed.i5)}
    
    adapter_length = {1: len(arguments_parsed.adapter_r1), 
                      2: len(arguments_parsed.adapter_r2)}

    with open(arguments_parsed.truth, newline="") as file:
        reader = csv.DictReader(file, delimiter="\t")
        fieldnames = reader.fieldnames
        rows = list(reader)
    truth_by_name = {}
    for row in rows:
        truth_by_name[read_name(row["read_pair_id"])] = row

    for mate, input_path in ((1, arguments_parsed.r1), (2, arguments_parsed.r2)):
        random_generator = random.Random(arguments_parsed.seed * 10 + mate)
        output_path = f"{arguments_parsed.out_prefix}_r{mate}.fastq"
        if input_path.endswith(".gz"):
            output_path += ".gz"
        with open_file(output_path, "wt") as output:
            for header, seq, plus, quality in fastq(input_path):
                row = truth_by_name.get(read_name(header))
                if row is not None:
                    genomic_seq = row[f"r{mate}_genomic_seq"]
                    adapter_seq = row[f"r{mate}_adapter_seq"]
                    if adapter_seq and len(adapter_seq) > adapter_length[mate]:
                        first_adapter_copy = adapter_seq[:adapter_length[mate]]
                        new_adapter_seq = first_adapter_copy + continuation(
                            read_through[mate], len(adapter_seq) - adapter_length[mate],
                            arguments_parsed.tail, arguments_parsed.error_rate, random_generator)
                        seq = genomic_seq + new_adapter_seq
                        row[f"r{mate}_adapter_seq"] = new_adapter_seq
                output.write(f"{header}\n{seq}\n{plus}\n{quality}\n")

    with open(f"{arguments_parsed.out_prefix}_parsed_adapters.tsv", "w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

if __name__ == "__main__":
    main()