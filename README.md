[![GitHub release](https://img.shields.io/github/v/release/abjanssen/Readzor?label=Latest%20release)](https://github.com/abjanssen/Readzor/releases)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Python3](https://img.shields.io/badge/Language-Python3-steelblue.svg)](https://www.python.org/)
[![DOI](https://zenodo.org/badge/1333429627.svg)](https://doi.org/10.5281/zenodo.22336649)
[![PyPI Downloads](https://img.shields.io/pepy/dt/readzor?label=PyPI%20downloads)](https://pypi.org/project/Readzor/)
[![Bioconda Downloads](https://img.shields.io/conda/dn/bioconda/Readzor?label=Bioconda%20downloads&color=%2352b924)](https://anaconda.org/channels/bioconda/packages/readzor/overview)
[![Docker Downloads](https://img.shields.io/badge/dynamic/json?url=https://ghcr-badge.elias.eu.org/api/abjanssen/Readzor/readzor&query=downloadCount&label=Docker%20downloads&color=%2352b924)](https://github.com/abjanssen/Readzor/pkgs/container/readzor)

# Welcome to Readzor

<img src="docs/other/logo.png" alt="Readzor logo" width="256" height="256" align="left">

**Readzor** is a fast, modular, and fully-featured FASTQ quality trimming and filtering pipeline, written in Python.
Readzor can trim and filter FASTQ reads for quality, adapters, homopolymers, Ns and low-complexity. 
Designed with the user in mind, Readzor will detect read-pairing, Phred offsets, file compression formats, and integrates with HPC clusters.

By default, all processing modules are **off**, giving you full control over your workflow.

<br clear="left">

## Contents
- [How it works](#how-it-works)
- [Features](#features)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Command line options](#command-line-options)
- [Benchmarking](#benchmarking)
- [Dependencies](#dependencies)
- [Issues and bug reports](#issues-and-bug-reports)
- [License](#license)
- [Author](#author)
- [Reference](#reference)
- [Oligonucleotide sequences](#oligonucleotide-sequences)

## How it works
Each enabled module assesses every raw read independently. The results are then merged, and the most stringent trim per read end is applied. For example, if quality trimming would remove 5 bases from the 3' end and adapter trimming would remove 12, Readzor removes 12.
 
<img src="docs/other/summary.png" alt="Readzor summary" width="100%">

## Features
**Modular trimming approach**:
By design, each filter will independently assess each raw read. Outcomes are merged, and the most stringent trimming result per end is applied.\
For example, if the quality trimmer determines 5 bases should be removed from the 3' end, but the adapter trimmer identifies 12 bases to remove from that same end, Readzor will merge these outcomes and trim the most restrictive amount (12 bases) to ensure high-quality output.
*   Quality-dependent end trimming.
*   Sliding window quality trimming.
*   Homopolymer trimming, from both read ends.
*   Adapter trimming with built-in Illumina adapter groups (TruSeq, TruSeq small RNA, Illumina miRNA, Nextera, AmpliSeq, TruSeq DNA methylation, TruSeq Ribo-Profile, and Illumina RNA), with optional mismatches and detection of partial adapters at the 3' end.
*   Overlap-based adapter trimming for paired and interleaved reads, independent of the adapter sequence.
*   N-base trimming from both read ends, and overall N-base filtering.
*   K-mer based low complexity filtering.
*   Set-length read end trimming.
*   Average input and output quality thresholds.
*   Input and output length thresholds (absolute, or as a percentage of the input read length).
    
**Auto-detection**:
Using auto-detection methods, Readzor provides an easy-to-use platform for novice users, enabling high quality read trimming with minimal inputs.
*   File Pairing: paired vs. unpaired vs. interleaved FASTQ files detection using internal header information, independent of file names.
*   Read Numbering: read numbers detection when needed, directly from headers.
*   Compression: Gzip compression detection via magic bytes, regardless of the file extension.
*   Phred Offsets: per-file Phred quality offsets detection (i.e., Phred33 vs. Phred64) based on quality string symbols.

**Extra features**:
Equipped with a number of other features, like built-in validation and seamless HPC cluster detection methods, Readzor delivers a robust execution environment, ensuring reproducible outputs at any scale.
*   Fully automatic mode: invoking Readzor in fully automatic mode (optional file specification possible), will use predetermined settings to optimize read processing. For more info, invoke Readzor with flags --full-auto and --help.
*   MGI/BGI header conversion: Built-in module to convert MGI/BGI fastq headers into standard Illumina format for downstream compatibility (e.g., with SAMtools).
*   Safe and traceable outputs: Automatically generates isolated, timestamped output directories to prevent accidental file overwrites. Each run includes a log file recording the exact command and all parameters used, and a summary file with the kept and rejected read counts per input, for full reproducibility.
*   HPC integration: Automatically reads HPC cluster environment variables to seamlessly scale threads and optimize performance on HPC clusters.
*   Live progress monitoring: Features a dynamic progress tracker to accurately provide real-time feedback on processing speed, completion, and estimated time remaining.
*   Real-time record validation: Readzor detects malformed or corrupted FASTQ records on the fly, ensuring high-quality output.
*   Gzip compression: Compress your output files and control the compression depth to optimize storage size.
*   Test run: Perform a short test run according to specified settings, printing verbose output to screen, to catch config errors before committing to a full run. No output files are kept.
*   Piping: Read FASTQ data from stdin and stream the filtered reads to stdout, to use Readzor inside a pipeline.
*   Rejected reads: Optionally write all rejected reads to separate files for inspection.

**Input notes**: gzipped and plain files can be mixed freely. In `--full-auto` mode without input files, Readzor picks up files in the working directory ending in `.fastq` or `.fq` (optionally followed by `.gz` or `.gzip`). Extra adapter sequences are supplied as a FASTA file.

## Installation
Installation is made easy through several repositories. Check for correct installation by calling `--version`. 
 
**Bioconda**
```bash
conda install readzor
readzor --version
```
 
**Pip**
```bash
pip install readzor
readzor --version
```
 
**Docker**
```bash
docker pull ghcr.io/abjanssen/readzor:latest
docker run --rm ghcr.io/abjanssen/readzor:latest --version
```
 
**From source**
```bash
git clone https://github.com/abjanssen/Readzor.git
pip install numpy==2.5.2 isal==1.8.0
python Readzor/src/readzor/readzor.py --version
```

## Quick start
 
```bash
# Fully automatic: preset settings, files and pairing detected in the current directory
readzor -GO
 
# Automatic, but with your own files (pairing auto-detected)
readzor -GO --input-files sample*.fastq.gz
 
# Declare pairing explicitly (can be combined)
readzor -GO --input-paired s1_R1.fq.gz s1_R2.fq.gz --input-unpaired s2.fq.gz
 
# Choose modules yourself, with their default settings
readzor --input-files *.fq.gz --adapter-filter-flag --endqual-filter-flag
 
# Customize module settings
readzor --input-files *.fq.gz --endqual-filter-flag --endqual-min-both 25
readzor --input-files *.fq.gz --adapter-filter-flag --adapter-group Nextera --adapter-mismatch 1
readzor --input-files *.fq.gz --adapter-filter-flag --adapter-fasta-add my_adapters.fasta
 
# Stream through a pipeline
cat reads.fastq | readzor --endqual-filter-flag --stdout > filtered.fastq
```

## Command line options 
All modules are off by default. To use a module, specify a module flag. Further specification of module settings is possible using the appropiate options.

### General settings
General options for help, version information, automatic mode, and console output.\
`--help, -h` [FLAG]: Show this help message and exit. Combine with --full-auto/-GO for more information on fully automatic mode.\
`--version, -v` [FLAG]: Show Readzor version and exit.\
`--list-adapters` [FLAG]: Show all built-in adapter sequences and exit.\
`--full-auto, -GO` [FLAG]: Run Readzor in fully automatic mode. Combine with --help/-h for more information on fully automatic mode.\
`--progress` [FLAG]: Show a live progress bar and estimated time remaining during processing, based on estimated read counts. Default: off.\
`--verbose` [FLAG]: Write verbose output to terminal, in addition to the log file. Default: off.
    
### Input options
Specify input FASTQ files using any combination of --input-files, --input-paired, --input-unpaired, and --input-interleaved. Lists with any combination of regular (fastq/fq) and gzipped (fastq.gz/fq.gz) files accepted. If no input flags are provided, data will be read directly from standard input (stdin). With --full-auto and no input flags, FASTQ files in the current working directory are used instead.\
`--input-files, -i`: FASTQ files of unspecified pairing. Paired, interleaved, and unpaired files will be auto-detected.\
`--input-interleaved, -ii`: Interleaved FASTQ files. Note: these files will be split in forward and reverse reads.\
`--input-paired, -ip`: Paired-end FASTQ files, given as one or more R1/R2 pairs, e.g. --input-paired sample1_R1 sample1_R2 sample2_R1 sample2_R2.\
`--input-unpaired, -iu`: Unpaired FASTQ files.

### Output options
Filtered sequences can be either written to file, or output to stdout. Setting --stdout overrides other settings.\
`--gzip` [FLAG]: Compress filtered FASTQ files in gzip format using isal. Default: off.\
`--gzip-level`: Set gzip compression level. Higher compression decreases processing speed. Possible values: 0-3. Default: 1.\
`--interleaved-out` [FLAG]: Interleave surviving read pairs of paired and interleaved input files into one output file. Single surviving reads are still written to separate files. Default: off.\
`--output, -o`: Path to directory in which the timestamped results folder will be created. Default: current working directory.\
`--stdout` [FLAG]: Stream resulting FASTQ reads to stdout. Forces --interleaved-out and --discard-singles for paired and interleaved files. Overrides --verbose, --progress, and --write-rejected to off. Overridden to 'off' by --gzip. Note: all files will be streamed one after another, without any separators. Default: off.\
`--write-rejected` [FLAG]: Write rejected reads to file. Either one (for unpaired and when --interleaved-out is set), or two (for forward and reverse reads) are produced. Overridden to 'off' when --stdout is set. Default: off.\
`--discard-singles` [FLAG]: Discard single leftover reads. For paired and interleaved reads, single surviving reads will be discarded instead of written to a separate file. No effect on unpaired read filtering. Default: off.\
`--stats` [FLAG]: Write detailed summary read trimming statistics. Default: off.
             
### General quality filters
Filters raw input, or trimmed output sequences for overall parameters. When a read does not pass these filters, it is discarded entirely.\
`--min-average-qual-pre`: Minimum average quality of input read. Possible values: 0-127. Default: 0.\
`--min-average-qual-post`: Minimum average quality of output read. Possible values: 0-127. Default: 0.\
`--min-length-input`: Minimum length for input read. Default: off.\
`--max-length-input`: Maximum length for input read. Default: off.\
`--min-length-output`: Minimum length of output read in absolute number of nucleotides. Overrides --min-length-output-perc when both set. Default: off.\
`--min-length-output-perc`: Minimum length of output read as percentage of input read. Overridden by --min-length-output when both set. Default: off.\
`--max-length-output`: Maximum length of output read in absolute number of nucleotides. Overrides --max-length-output-perc when both set. Default: off.\
`--max-length-output-perc`: Maximum length of output read as percentage of input read. Overridden by --max-length-output when both set. Default: off.\
`--n-filter` [FLAG]: Reject reads containing `N` bases anywhere in read. Turns off N nucleotide end trimming. Default: off.

### Set-Length end trimming
Trim a set number of bases of the ends of each read, independent of sequence or quality.\
`--cut-flag, -cf` [FLAG]: Turn on the set-length end trimming module. Default: off.\
`--cut-start, -cs`: Number of bases to trim from the start of the read. Default: 0.\
`--cut-end, -ce`: Number of bases to trim from the end of the read. Default: 0.\
`--cut-both, -cb`: Number of bases to trim from both ends of the read. Overwritten by --cut-start and --cut-end. Default: 0.

### Quality-dependent end trimming
Trim the ends of each read for quality, up to the first position that fulfills quality requirements.\
`--endqual-filter-flag, -ef` [FLAG]: Turn on quality-dependent end trimming. Default: off.\
`--endqual-min-start, -ems`: Specific phred score threshold for the start of the read. Default: 20.\
`--endqual-min-end, -eme`: Specific phred score threshold for the end of the read. Default: 20.\
`--endqual-min-both, -emb`: Phred score threshold for the quality trimming of read ends. Overwritten by --endqual-min-start and --endqual-min-end. Default: 20

### Sliding window quality trimming
Trim the reads for quality based on a sliding window of size X, moved with stepsize Y. In case of mid-read quality dropoff, the longest portion survives. In case of a tie, the portion with the highest average quality survives.\
`--slider-filter-flag, -sf` [FLAG]: Turn on sliding window quality trimming module. Default: off.\
`--slider-window, -sw`: Window size over which average quality is calculated. Default: 5.\
`--slider-quality, -sq`: Minimum average quality in sliding window. Default: 20.\
`--slider-step, -ss`: Sliding window step size. Default: 1.

### Homopolymer nucleotide trimming
Illumina NovaSeq, NextSeq, and MiniSeq use a two-color chemistry, in which guanine bases are unlabeled. In event of short fragments, this can result in homopolymer G calls at the end of reads. This filter allows for G, and/or other bases, to be cut.\
`--poly-filter-flag, -pf` [FLAG]: Turn on homopolymer read-end trimming module. Default: off.\
`--poly-bases-start, -pbs`: Base(s) to check for a homopolymer run at start of read. Comma-separated bases are checked independently. Default: none.\
`--poly-bases-end, -pbe`: Base(s) to check for a homopolymer run at end of read. Comma-separated bases are checked independently. Default: "G".\
`--poly-bases-both, -pbb`: Base(s) to check for a homopolymer run at both read ends. Comma-separated bases are checked independently. Overwritten by --poly-bases-start and --poly-bases-end. Default: none.\
`--poly-length-start, -pls`: Minimum length of homopolymer run at start of read required to trigger trimming. A value of 0 disables trimming at the start. Default: off.\
`--poly-length-end, -ple`: Minimum length of homopolymer run at end of read required to trigger trimming. A value of 0 disables trimming at the end. Default: 10. \
`--poly-length-both, -plb`: Minimum length of homopolymer run at start and end of read required to trigger trimming. Used for an end whose own length is not given. Default: off.

### Adapter trimming
Trim reads for adapter sequences. Each read is cut at the earliest adapter occurrence; adapters running off the 3' end are detected down to --adapter-seed bases. For paired reads, only the adapter belonging to each mate (Read 1 or Read 2) of a built-in group is searched for. Use --list-adapters to show all built-in sequences. Adapter trimming is performed independent of quality. Can be combined with overlap trimming (--overlap-filter-flag).\
`--adapter-filter-flag, -af` [FLAG]: Turn on adapter trimming module. Default: off.\
`--adapter-group, -ag`: Specify the group(s) of adapters to be used. Ignored if --adapter-fasta-excl is set. Choices: all, TruSeq, TruSeq_small_RNA, Illumina_miRNA, Nextera, AmpliSeq, TruSeq_DNA_methylation, TruSeq_Ribo_profile, Illumina_RNA. Default: TruSeq.\
`--adapter-seed, -as`: Minimal 5'-end match length. Shorter seeds can result in more partial hits found at end of reads. Default: 6.\
`--adapter-mismatch, -am`: Number of mismatches allowed in adapter finding (substitutions only). For partial adapters at the 3' end, the allowed number is scaled down to the matched length. Default: 2.\
`--adapter-fasta-add, -ad`: Fasta file with adapter sequences to trim for, in addition to defined sequences from --adapter-group.\
`--adapter-fasta-excl, -ax`: Fasta file with adapter sequences to trim for, excluding predefined and additional sequences specified.

### Overlap trimming
Perform overlap analysis to trim adapter sequences. Works independently of adapter sequence. Only available for paired and interleaved reads. Can be combined with --adapter-filter-flag.\
`--overlap-filter-flag, -of` [FLAG]: Turn on overlap trimming module. Default: off.\
`--overlap-mismatches, -om`: Maximum mismatches allowed in the overlapping region. Default: 5.\
`--overlap-min-length, -ol`: Minimum overlap length between the paired reads. Default: 50.\
`--overlap-portion-mismatch, -op`: Maximum percentage of mismatched bases allowed in the overlapping region. Default: 4
    
### Low complexity filtering
Detect complexity of reads using k-mer-based nucleotide frequencies. Low complexity reads are discarded entirely.\
`--kmer-filter-flag, -kf` [FLAG]: Turn on the k-mer-based complexity filtering module. Default: off.\
`--kmer-size, -ks`: K-mer length for k-mer-based complexity filtering. Comma-separated values are checked independently. Maximum value: 21. Default: 4.\
`--kmer-cutoff, -kc`: Minimum percentage of unique k-mers (relative to the maximum possible for the read) required to pass the complexity filter. Higher values are stricter. Default: 50.

### MGI header conversion
Convert read header from MGI (BGI) format to Illumina format. Original header will be discarded.\
`--mgi-convert-flag, -mf` [FLAG]: Turn on the MGI-to-Illumina header conversion module. Default: off.\
`--mgi-bc5, -m5`: Input an i5 barcode for Illumina header conversion. Default: 'PLACEHOLDERi5'.\
`--mgi-bc7, -m7`: Input an i7 barcode for Illumina header conversion. Default: 'PLACEHOLDERi7'.\
`--mgi-instrument, -mi`: Instrument name for Illumina header conversion. Default: 'PLACEHOLDERinstrument'.\
`--mgi-run, -mr`: Run ID for Illumina header conversion. Default: 'PLACEHOLDERrun'.

### Advanced options
Further options that can be specified to alter the behavior of Readzor.\
`--threads, -t`: Number of threads to use. Default: platform-dependent through auto-detection (assigned CPUs on HPC cluster systems, all-1 otherwise. Fallback: 1).\
`--reads-for-phred-offset`: Number of reads to sample per file for detection of Phred quality encoding offset. Default: 500.\
`--chunk-size`: Number of reads per chunk sent to each worker. Changing can alter processing speed. Default: 1500.\
`--phred-offset`: Define phred offset for all FASTQ files. When set, per-file auto-detection will not be performed. Possible values: 33, 64. Default: off (auto-detection per file).\
`--testrun` [FLAG]: Perform a limited test run according to specified settings: only the first chunks of each input are processed (one per thread), and all output, including the log and summary files, is deleted afterwards. Also sets --verbose. Default: off.\
`--ordered-output` [FLAG]: Force writing output reads in the same order as the input reads, with the exception of filtered reads. Default: off.\
`--phred-out`: Convert Phred encoding from 33 to 64, and vice versa. Possible values: 33, 64. Default: off.

## Benchmarking

### Production of synthetic short read sequences
To provide for a ground-truth for adapter removal, we produced a synthetic dataset from the _Streptococcus pneumoniae_ BHN418 genome using [NEAT](https://github.com/ncsa/NEAT) (NCBI Accession [CP155538.1](https://www.ncbi.nlm.nih.gov/nuccore/CP155538.1)).
We produced a 151 bp, paired-end dataset, with a mean fragment size of 280 bp and standard deviation of 50 bp, with TruSeq adapters ligated where necessary. 
This dataset totals 56,009,434 reads, with 267,114 (0.476%) reads containing adapter sequences.
The exact adapter content was parsed for use in adapter removal benchmarking using the information in the read headers and the golden BAM file.

### Adapter sequence trimming
Readzor can use two methods to find and trim adapter sequences: sequence search and overlap search.
To assess how these methods performed, we determined the following outcomes (amongst others):
* Recall (portion of adapters succesfully found);
* Precision (portion of removed bases that belonged to adapter sequences);
* F1 (harmonic mean combining precision and recall, $`F_1 = 2 \cdot \frac{\text{precision} \cdot \text{recall}}{\text{precision} + \text{recall}}`$);
* Genomic basese lost.

#### Adapter sequence search approach
Adapter sequence search relies on the specification of known adapter sequences.
Readzor scans each read for the first `--adapter-seed` bases of each adapter, adapter fragments shorter than the seed are therefore not found. 
The mismatches allowed are scaled proportionally according to remaining number of nucleotides in the read.

A grid search identified `--adapter-seed 6` and `--adapter-mismatch 2` as the best balance for precision and recall.

<h3 align="center">Adapter sequence search performance</h3>
<p align="center">
  <img src="docs/benchmarking/adapter_trimming/readzor/adapter_sequence_search_recall.png" width="49%">
  <img src="docs/benchmarking/adapter_trimming/readzor/adapter_sequence_search_precision.png" width="49%">
  <img src="docs/benchmarking/adapter_trimming/readzor/adapter_sequence_search_f1.png" width="49%">
  <img src="docs/benchmarking/adapter_trimming/readzor/adapter_sequence_search_genomic_bases_lost.png" width="49%">
</p>

#### Overlapping sequences approach
For paired-end data, when the insert is shorter than the read length, mates read into the adapters, which can be found by aligning the mates to each other. 
Adapter sequences can thus be identified without requiring specification of adapter sequences.
Inserts shorter than `--overlap-min-length` are not found.

By setting `--overlap-portion-mismatch` to `100`, we isolated the effect of `--overlap-mismatches` on the adapter sequences found.
Vice versa, by setting `--overlap-mismatches` to `999`, the behaviour of `--overlap-portion-mismatch`, was investigated.

<h3 align="center">Overlapping sequences search performance: absolute number of mismatches</h3>
<p align="center">
  <img src="docs/benchmarking/adapter_trimming/readzor/overlap_mismatch_recall.png" width="49%">
  <img src="docs/benchmarking/adapter_trimming/readzor/overlap_mismatch_precision.png" width="49%">
  <img src="docs/benchmarking/adapter_trimming/readzor/overlap_mismatch_f1.png" width="49%">
  <img src="docs/benchmarking/adapter_trimming/readzor/overlap_mismatch_genomic_bases_lost.png" width="49%">
</p>

<h3 align="center">Overlapping sequences search performance: relative number of mismatches</h3>
<p align="center">
  <img src="docs/benchmarking/adapter_trimming/readzor/overlap_portion_recall.png" width="100%">
  <img src="docs/benchmarking/adapter_trimming/readzor/overlap_portion_precision.png" width="100%">
  <img src="docs/benchmarking/adapter_trimming/readzor/overlap_portion_f1.png" width="100%">
  <img src="docs/benchmarking/adapter_trimming/readzor/overlap_portion_genomic_bases_lost.png" width="100%">
</p>

From both these datasets, we observe an inverse correlation between the precision and F1, and number of allowed mismathches, through the rapid increase in genomic bases lost.
A minimal `--overlap-min-length` of`50` seems to strike a balance between different outcomes.

To further characterize the behaviour of this trimming method, we set the `--overlap-min-length` at `50`, and varied the mismatch parameters.
As Readzor uses the lowest number of these two settings to set the allowed number of mismatches for an identified overlap, we can finetune the overlap search.
For shorter overlaps, `--overlap-portion-mismatch` will determine the number of mismatches allowed, but at longer overlaps, mismatches are capped by `--overlap-mismatches`, preventing widespread aspecific matching.
As the combination of these two methods determine the allowed number of mismatches, a grid-search only needs to perform a limited number of searches, as one of both parameters acts as _de facto_ cap.

<h3 align="center">Overlapping sequences search performance: fixed length 50</h3>
<p align="center">
  <img src="docs/benchmarking/adapter_trimming/readzor/overlap_length50_recall.png" width="100%">
  <img src="docs/benchmarking/adapter_trimming/readzor/overlap_length50_precision.png" width="100%">
  <img src="docs/benchmarking/adapter_trimming/readzor/overlap_length50_f1.png" width="100%">
  <img src="docs/benchmarking/adapter_trimming/readzor/overlap_length50_genomic_bases_lost.png" width="100%">
</p>

From these results, the default `--overlap-min-length` has been set to `50`, `--overlap-mismatches` to `5`, and `--overlap-portion-mismatch` to `4`.
These settings balance adapter recall, precision, and loss of genomic base pairs.

#### Combined adapter and overlap approach
When desired, both adapter sequence search and overlap search can be combined to further optimize adapter trimming.
For this synthetic dataset produced, the combination of these methods results in the following performance:

<div align="center">

| Adapter trimming method  | Recall | Precision |    F1 | Genomic bases lost |
|:-------------------------|-------:|----------:|------:|-------------------:|
| Sequence search          |  0.951 |     0.993 | 0.971 |             27,172 | 
| Overlap search           |  0.982 |     0.780 | 0.869 |          1,202,900 | 
| **Combined search**      |  0.999 |     0.779 | 0.875 |          1,230,072 |  

</div>

The preset parameters strike a balance between the portion of adapters found, precision, and genomic bases lost. 
When a more rigorous adapter trimming is desired, the parameters can be adapted. 
A more detailed overview of the effect of different settings on these parameters can be found in the [benchmarking folder](docs/benchmarking/adapter_trimming/readzor) in this repository. 

#### Comparison to other short-read trimming tools
We compared the performance of Readzor's adpater trimming capabilities, to other popular read-trimming tools.

As different tools can have different preset parameters, we ran these tools either with preset parameters, or by adapting their parameters to mimick each other (when possible).
The commands used are listed in the [benchmarking folder](docs/benchmarking/adapter_trimming/other) in this repository. 
We provided each tool with the same adapter sequences, to further standardize adapter trimming.

#### Summary
**Paired-end**

| Tool        | Configuration                         | Recall | Precision |    F1 | Genomic bases lost |
|:------------|:--------------------------------------|-------:|----------:|------:|-------------------:|
| Readzor     | Specified adapters + overlap          |  0.999 |     0.779 | 0.875 |          1,230,072 |
| fastp       | Readzor adapters (per read) + overlap |  0.998 |     0.622 | 0.766 |          2,645,431 |
| Cutadapt    | Readzor adapters (per read)           |  0.988 |     0.467 | 0.634 |          4,914,110 |
| Trimmomatic | Readzor adapters, simple + palindrome |  0.918 |     0.994 | 0.955 |             20,867 |

**Single-end**

| Tool        | Configuration                         | Recall | Precision |    F1 | Genomic bases lost |
|:------------|:--------------------------------------|-------:|----------:|------:|-------------------:|
| Cutadapt    | Readzor adapters (per read)           |  0.988 |     0.467 | 0.634 |          4,914,110 |
| fastp       | Readzor adapters (per read)           |  0.968 |     0.921 | 0.944 |            360,242 |
| Readzor     | Specified adapters                    |  0.951 |     0.993 | 0.971 |             27,172 |
| Trimmomatic | Readzor adapters, simple              |  0.710 |     1.000 | 0.830 |                  0 |

#### Full results
<details>
<summary><b>fastp results</b></summary>

| Mode | Adapters                             | Overlap | Recall | Precision |        F1 | Genomic bases lost |
|:-----|:-------------------------------------|:--------|-------:|----------:|----------:|-------------------:|
| PE   | —                                    |   Yes   |  0.994 |     0.654 |     0.789 |          2,285,459 |
| PE   | Auto-detection¹                      |   Yes   |  0.994 |     0.654 |     0.789 |          2,285,459 |
| PE   | Readzor (per read)                   |   Yes   |  0.998 |     0.622 |     0.766 |          2,645,431 |
| PE   | Readzor (overall)                    |   Yes   |  0.998 |     0.616 |     0.761 |          2,712,500 |
| PE   | Auto-detection¹                      |   No²   |  0.000 |         — |     0.000 |                  0 |
| PE   | Auto-detection¹ + Readzor (per read) |   No²   |  0.968 |     0.921 |     0.944 |            360,242 |
| PE   | Auto-detection¹ + Readzor (overall)  |   No²   |  0.968 |     0.908 |     0.937 |            426,429 |
| SE   | Auto-detection¹                      |   No    |  0.000 |         — |     0.000 |                  0 |
| SE   | Readzor (per read)                   |   No    |  0.968 |     0.921 |     0.944 |            360,242 |
| SE   | Readzor (overall)                    |   No    |  0.968 |     0.908 |     0.937 |            426,429 |

</details>

<details>
<summary><b>Cutadapt results</b></summary>

| Mode | Adapters           | Recall | Precision |        F1 | Genomic bases lost |
|:-----|:-------------------|-------:|----------:|----------:|-------------------:|
| PE   | Readzor (per read) |  0.988 |     0.467 |     0.634 |          4,914,110 |
| PE   | Readzor (overall)  |  0.988 |     0.466 |     0.634 |          4,914,305 |
| SE   | Readzor (per read) |  0.988 |     0.467 |     0.634 |          4,914,110 |
| SE   | Readzor (overall)  |  0.988 |     0.466 |     0.634 |          4,914,305 |

</details>

<details>
<summary><b>Trimmomatic results</b></summary>

| Mode | Adapters     | Clipping mode         | Recall | Precision |        F1 | Genomic bases lost |
|:-----|:-------------|:----------------------|-------:|----------:|----------:|-------------------:|
| PE   | Readzor      | Simple                |  0.710 |     1.000 |     0.830 |                  0 |
| PE   | Readzor      | Palindrome            |  0.907 |     0.994 |     0.949 |             20,867 |
| PE   | Readzor      | Simple + palindromic  |  0.918 |     0.994 |     0.955 |             20,867 |
| PE   | TruSeq3-SE   | Simple                |  0.710 |     0.998 |     0.830 |                387 |
| PE   | TruSeq3-PE   | Palindromic           |  0.907 |     0.995 |     0.949 |             20,867 |
| PE   | TruSeq3-PE-2 | Simple + palindromic  |  0.918 |     0.995 |     0.955 |             21,466 |
| SE   | Readzor      | Simple                |  0.710 |     1.000 |     0.830 |                  0 |
| SE   | Readzor      | Palindromic³          |  0.000 |         — |     0.000 |                  0 |
| SE   | Readzor      | Simple + palindromic³ |  0.710 |     1.000 |     0.830 |                  0 |
| SE   | TruSeq3-SE   | Simple                |  0.710 |     0.999 |     0.830 |                387 |
| SE   | TruSeq3-PE   | Palindromic³          |  0.000 |         — |     0.000 |                  0 |
| SE   | TruSeq3-PE-2 | simple + palindromic³ |  0.710 |     0.999 |     0.830 |                599 |

</details>

<details>
<summary><b>Readzor results</b></summary>

| Mode | Adapters             | Overlap | Recall | Precision |        F1 | Genomic bases lost |
|:-----|:---------------------|:--------|-------:|----------:|----------:|-------------------:|
| PE   | Specified per read   |   No    |  0.951 |     0.993 |     0.971 |             27,172 |
| PE   | Specified overall    |   No    |  0.951 |     0.993 |     0.971 |             27,172 |
| SE   | Specified per read   |   No    |  0.951 |     0.993 |     0.971 |             27,172 |
| SE   | Specified overall    |   Yes   |  0.951 |     0.993 |     0.971 |             27,172 |
| PE   | Specified per read   |   Yes   |  0.999 |     0.779 |     0.875 |          1,230,072 |
| PE   | Specified overall    |   Yes   |  0.999 |     0.779 |     0.875 |          1,230,072 |
| PE   | —                    |   Yes   |  0.982 |     0.780 |     0.869 |          1,202,900 |

</details>

¹ fastp's adapter auto-detection found no adapters in this dataset.
² Overlap trimming disabled with `--overlap_len_require 999`.
³ Palindrome (overlap) clipping needs paired reads, so it has no effect in SE mode.

### Resource usage across processing modules
The impact of individual Readzor modules on processing time and memory usage was determined, and runs without any module (baseline) and with all modules combined were included for reference. 
Without any module, Readzor's processing time and memory usage are spent on unpacking, reading, pickling, validating, and compressing, writing files. 
To compare the different strategies that adapter sequence search uses for mismatch-allowed search, we allowed for 0-5 mismatches.
Most modules have a low processing burden, adding none to limited extra time and memory for processing. 
Sliding window quality trimming, overlap trimming and k-mer low complexity filtering, added substantial time to the processing.
When all modules are active, processing takes about 2.5× longer than the baseline scenario; roughly the sum of the individual module costs. 

Memory usage is largely consistent across modules.

<p align="center">
  <img src="docs/benchmarking/resource_usage_modules/resource_usage_modules.png" width="100%">
</p>

### File compression
The difference in performance between gzip and isal was determined using a Readzor version in which both were available.
Isal drastically reduces processing time, both when reading compressed files and compressing results, although peak memory usage increases slightly.
Because of this speed-up, Readzor uses [`isal`](https://github.com/pycompression/python-isal) for its high-performance implementation of the DEFLATE algorithm.
The default compression level for when `--gzip` is set, has been set to `--gzip-level 1`, as no appreciable reduction in file size was observed when set higher.
When aiming for comparable output file size (isal level `1` versus gzip level `5`), processing time is reduced about 2.2x (compressed files in- and output), with comparable memory usage.

<p align="center">
  <img src="docs/benchmarking/compression_methods/compression_methods.png" width="100%">
</p>

### Processing time and memory usage (compared to other tools)

Finally, we characterized the processing speed and maximum memory usage of Readzor and other popular short read trimming tools. 
For this comparison, we used multiple real-world sequencing data files, to mimick the actual usage of Readzor.

Readzor has built-in functionality to process multiple files. For fastp, we used [this](https://github.com/OpenGene/parallel/blob/main/fastp.py) recently published script to mimick this behaviour.
Cutadapt and Trimmomatic do not have the same functionality, and were run in a loop stead.

For the paired read comparison, whole-genome sequencing reads were obtained from BioProject [PRJEB74511](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB74511/) (ERR12875316, ERR13030616-ERR13030639).
For the comparison for unpaired reads, we used the _S. pneumoniae_ D39V data available from BioProject [PRJNA1347041](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1347041/) (SRX30876273-SRX30876312)

#### Paired read comparison
<p align="center">
  <img src="docs/benchmarking/comparison_tools/paired/paired_comparison.png" width="100%">
</p>

#### Unpaired read comparison
<p align="center">
  <img src="docs/benchmarking/comparison_tools/unpaired/unpaired_comparison.png" width="100%">
</p>

### Software and facilities
Readzor has been benchmarked using the zen3 nodes of the [Curnagl](https://wiki.unil.ch/ci/books/high-performance-computing-hpc/page/curnagl) high-performance computering cluster of the University of Lausanne.
For benchmarking, the following software packages and versions were used, all read trimmer were installedin their own dedictated conda environment:
* [Trimmomatic](https://github.com/usadellab/trimmomatic) (version 0.41)
* [fastp](https://github.com/opengene/fastp) (version 1.3.6)
* [Cutadapt](https://github.com/marcelm/cutadapt/) (version 5.2)
* [NEAT](https://github.com/ncsa/NEAT) (version 4.7.1)
* [psutil](https://github.com/giampaolo/psutil) (version 7.2.2)
* [python-isal](https://github.com/pycompression/python-isal) (version 1.8.0)
* [NumPy](https://github.com/numpy/numpy) (version 2.5.2)

All commands used for benchmarking, and subsequent results, are available from the [benchmarking folder](docs/benchmarking) in this repository. 

## Dependencies
Readzor depends on the following software packages:
* NumPy (version >= 2.5.2)
* Python (no version requirement)
* Python-isal (version >= 1.8.0)

Although Readzor works on older versions of Python, it has been developed and tested for performance using version 3.14.7.

## Issues and bug reports
Please leave a message in [issues](https://github.com/abjanssen/Readzor/issues) or [discussions](https://github.com/abjanssen/Readzor/discussions) if you notice an issue, bug, or otherwise. 

## License
This project is provided under the GNU General Public License v3.0 (GPLv3).

## Author
Axel B. Janssen ([Google Scholar](https://scholar.google.com/citations?user=TWi-ysEAAAAJ&hl)/[GitHub](https://github.com/abjanssen/))

## Reference
Please use the DOI (https://doi.org/10.5281/zenodo.22336649) for your references. 

## Oligonucleotide sequences
The oligonucleotide sequences used for adapter trimming, included in this work, are copyrighted and protected by intellectual property, including issued or pending patents, copyright, and trade secrets.\
Oligonucleotide sequences © 2026 Illumina, Inc. All rights reserved.