#!/usr/bin/env bash

printf '%s\n' \
  "reference: BHN418.fna" \
  "read_len: 151" \
  "ploidy: 1" \
  "coverage: 4000" \
  "threads: 8" \
  "paired_ended: True" \
  "fragment_mean: 280" \
  "fragment_st_dev: 50" \
  "produce_bam: True" \
  "produce_vcf: False" \
  "produce_fastq: True" \
  "adapters: True" \
  "adapter_preset: truseq" \
  "rng_seed: 6443347027567828" \
  > config.yml

neat read-simulator --config config.yml --output ./output --prefix BHN418_synthetic
samtools view -h output/BHN418_synthetic_golden.bam > output/BHN418_synthetic_golden.sam
