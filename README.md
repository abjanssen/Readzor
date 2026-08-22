
[![GitHub release](https://img.shields.io/github/v/release/abjanssen/readzor)](https://github.com/abjanssen/readzor/releases)
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![Yes, it's in Python'](https://img.shields.io/badge/Language-Python3-steelblue.svg)
[![PyPI Downloads](https://img.shields.io/pypi/dm/readzor)](https://pypi.org/project/Readzor/)
[![Bioconda Downloads](https://img.shields.io/conda/dn/bioconda/prokka.svg)](https://bioconda.github.io/recipes/prokka/README.html)

# Welcome to ReadZor 

ReadZor is a fast, modular, and fully-featured FASTQ quality trimming and filtering pipeline. 
By default, **<ins>all processing modules are off</ins>**, leaving you the explicit ability to adapt the workflow according to your needs and wishes.

## Feature overview
**Modular trimming approach**:
By design, each filter will indepedently assess each raw read. Outcomes are merged, and the most stringent trimming result per end is applied.
For example, if the quality trimmer determines 5 bases should be removed from the 3' end, but the adapter trimmer identifies 12 bases to remove from that same end, ReadZor will merge these outcomes and trim the most restrictive amount (12 bases) to ensure high-quality output.
*   Quality-dependent end trimming
*   Sliding windowquality trimming
*   Homopolymer trimming from both read ends
*   Adapter trimming for TruSeq and Nextera adapters.
*   N-base trimming from read ends, or in overall.
*   Low complexity filtering based on K-mer counts.
*   Set-length read end trimming.
*   Overall quality threshold.
*   Overall length threshold    
    
**Auto-detection**:
Using auto-detection methods, ReadZor provides an easy-to-use platform for novice users, enabling high quality read trimming with minimal inputs.
*   File Pairing: paired vs. unpaired FASTQ files detection using internal header information, independent of file names.
*   Read Numbering: read numbers detection, directly from headers.
*   Compression: GZIP compression detection via magic bytes, regardless of the file extension.
*   Phred Offsets: per-file Phred quality offsets detection (i.e., Phred33 vs. Phred64) based on quality string symbols.

**Extra features**:
Equipped with a number of other features, like built-in validation and seamless Slurm detection methods, ReadZor delivers a robust execution environment, ensuring reproducible outputs at any scale.
*   Fully automatic mode: invoking ReadZor in fully automatic mode (optional file specification possible), will use predetermined settings to optimize read processing. For more info, invoke ReadZor with flags --full-auto and --help.
*   MGI/BGI header conversion: Built-in module to convert MGI/BGI fastq headers into standard Illumina format for downstream compatibility (e.g., with SAMtools).
*   Safe and traceable outputs: Automatically generates isolated, timestamped output directories to prevent accidental file overwrites. Each run includes a comprehensive summary file detailing the exact parameters used and the final trimming outcomes for full reproducibility.
*   HPC integration: Automatically reads Slurm environment variables to seamlessly scale threads and optimize performance on HPC clusters.
*   Live progress monitoring: Features a dynamic progress tracker to accurately provide real-time feedback on processing speed, completion, and estimated time remaining.
*   Real-time record validation: ReadZor detects malformed or corrupted FASTQ records on the fly, ensuring high-quality output.
*   Gzip compression: Compress your output files and control the compression depth  to optimize storage size.

## Installation and testing
Installation of Readzor is made easy through pip and conda, but you can also clone this repository:

### Pip
```python
pip install readzor
import readzor
print(readzor.VERSION)
```

### Repository cloning
```bash
git clone https://github.com/abjanssen/Readzor.git
pip install numpy
pip install isal
cd Readzor/src/reazor
chmod u+x readzor.py
python readzor.py --version
```

## Using Readzor 

### Beginner
```
# Use predetermined settings, and let Readzor detect the files (and their pairing) in your current working directory
% readzor -GO
```

### Beginner+
```
# Use predetermined settings, but specify your files

# Use on Readzor's paired-filed detection:
% readzor -GO --input-files /path/to/input/files

# Specify paired files:
% readzor -GO --input-paired /path/to/paired/files

# Specify unpaired files:
% readzor -GO --input-unpaired /path/to/unpaired/files

# Or a combination of both: 
% readzor -GO --input-unpaired /path/to/unpaired/files --input-paired /path/to/paired/files
```

### Novice
```
# Adapt the workflow according to your need by specifying specific trimming modules

# Turn on quality-dependent end trimming, use its default settings
% readzor --input-files /path/to/input/files --endqual-filter-flag

# Turn on adapter trimming, use its default settings
% readzor --input-files /path/to/input/files --adapter-trim-flag

# Or use a combination
% readzor --input-files /path/to/input/files --adapter-trim-flag --endqual-filter-flag
```

### Advanced
```
# Detail the module settings according to your needs by specifying trimming parameters

# Turn on quality-dependent end trimming, use custom settings
% readzor --input-files /path/to/input/files --endqual-filter-flag --endqual-min-start 25 --endqual-min-end 25

# Turn on adapter trimming, add custom sequences
% readzor --input-files /path/to/input/files --adapter-trim-flag --adapter-fasta /path/to/fasta/file
```

## Command line options 
All modules are off by default. To use a module, specify a module flag. Further specifications with module settings possible.

### General settings:
`--help/-h` [FLAG]: Show this help message and exit. Combine with --full-auto/-GO for more information on fully automatic mode.\
`--version/-v` [FLAG]: Show ReadZor version and exit.\
`--full-auto, -GO` [FLAG]: Run Readzor in fully automatic mode. Combine with --help/-h for more information on fully automatic mode.\
`--progress` [FLAG]: Show a live progress bar and estimated time remaining during processing, based on estimated read counts. Default: off.\

### Input options:
Specify input FASTQ files using any combination of --input-files, --input-paired, and --input-unpaired. Lists with any combination of regular (fastq/fq) and gzipped (fastq.gz/fq.gz) files accepted.\
`--input-files, -i`: FASTQ files of unspecified pairing. Paired and unpaired files will be auto-detected.\
`--input-paired, -ip `: Paired-end FASTQ files, given as one or more R1/R2 pairs, e.g. --input-paired sample1_R1 sample1_R2 sample2_R1 sample2_R2.\
`--input-unpaired, -iu`: Unpaired FASTQ files.\

###Output options:
`--output, -o`: Path to directory in which the timestamped results folder will be created. Default: current working directory.\
`--gzip` [FLAG]: Compress filtered FASTQ files in gzip format using isal. Default: off.\
`--gzip-level`: Set gzip compression level. Higher compression decreases processing speed. Possible levels: 0-3. Default: 2.\
             
### General quality filters
`--min-average-qual <int>`: Minimum average quality of output read. Default: 0.\
`--min-length <int>`: Minimum length of output read. Default: 0.\
`--max-length <int>`: Maximum length of output read. Default: off.\
`--nucl-filter` [FLAG]: Reject reads containing `N` bases anywhere in read. Default: off.\

### Set-Length end trimming
Trim a set number of bases of the ends of each read, independent of sequence or quality.\
`--cut-flag, -clf` [FLAG]: Turn on the set-length end trimming module. Default: off.\
`--cut-start, -cs <int>`: Number of bases to trim from the start of the read. Default: 0.\
`--cut-end, -ce <int>`: Number of bases to trim from the end of the read. Default: 0.\
`--cut-both, -cb  <int>`: Number of bases to trim from both ends of the read. Overwritten by --cut-start and --cut-end. Default: 0.\

### Quality-dependent end trimming
Trim the ends of each read, dependent on quality. Ends of reads will be trimmed up to first position that fulfills quality requirement.\
`--endqual-filter-flag, -eff`: [FLAG] Turn on quality-dependent end trimming. Default: off.\
`--endqual-min-start, -ems`: Specific phred score threshold for the start of the read. Default: 25.\
`--endqual-min-end, -eme`: Specific phred score threshold for the end of the read. Default: 25.\
`--endqual-min-both, -emb`: Phred score threshold for the quality trimming of read ends. Overwritten by --endqual-min-start and --endqual-min-end. Default: 25.\
                               
### N-Nucleotide end trimming
Trim the ends of each read for N bases. Redundant when --nucl-filter is set.\
`--n-end-trimming-flag, -ntf`: [FLAG]: Turn on the N nucleotide end trimming module. Default: off.\

### Sliding window quality trimming
Trim the reads for quality based on a sliding window of size X, moved with stepsize Y. Longest portion survives in case of mid-read quality dropoff.\
`--slider-filter-flag, -sf`: [FLAG] Turn on sliding window quality trimming module. Default: off.\
`--slider-window, -sw`: Window size over which average quality is calculated. Default: 5.\
`--slider-quality, -sq`: Minimum average quality in sliding window. Default: 20.\
`--slider-step, -ss`: Sliding window step size. Default: 1.\

### Homopolymer nucleotide trimming
Illumina NovaSeq, NextSeq, and MiniSeq use a two-color chemistry, in which guanine bases are unlabeled. In event of short fragments, this can result in homolopolymer G calls at the end of reads.\
`--poly-filter-flag, -pf`: [FLAG] Turn on homopolymer read-end trimming module. Default: off.\
`--poly-bases-start, -pbs`: Base(s) to check for a homopolymer run at start of read. Comma-separated bases are checked independently. Default: none.\
`--poly-bases-end, -pbe`: Base(s) to check for a homopolymer run at end of read. Comma-separated bases are checked independently. Default: "G".\
`--poly-bases-both, -pbb`: Base(s) to check for a homopolymer run at both read ends. Comma-separated bases are checked independently. Overwritten by poly_bases_start and poly_bases_end. Default: none.\
`--poly-length-start, -pls`: Minimum length of homopolymer run at start of read required to trigger trimming. Default: 10.\
`--poly-length-end, -ple`: Minimum length of homopolymer run at end of read required to trigger trimming. Default: 10.\
`--poly-length-both, -plb`: Minimum length of homopolymer run at start and end of read required to trigger trimming. Default: 0.\

### Adapter trimming
Trim reads for Illumina adapter sequences. Standard sequences included are TruSeq3 universal and index adapters, and Nextera adapters. Only perfectly matching sequences are trimmed. Independent of quality.\
`--adapter-trim-flag, -af`: [FLAG] Turn on adapter trimming module. Default: off.\
`--adapter-fasta, -ad`: FASTA file with additional adapter sequences to trim for.\

### Low complexity filtering
Detect complexity of reads using k-mer-based nucleotide frequencies. Low complexity reads are discarded entirely.\
`--kmer-filter-flag`: [FLAG] Turn on the k-mer-based complexity filtering module. Default: off.\
`--kmer-size`: K-mer length for k-mer-based complexity filtering. Comma-separated values are checked independently. Default: 4.\
`--kmer-cutoff`: Minimum percentage of unique k-mers (relative to the maximum possible for the read) required to pass the complexity filter. Higher values are stricter. Default: 50.\

### MGI header conversion
Convert read header from MGI (BGI) format to Illumina format. Original header will be stored in the placeholder line. Conversion is necessary for downstream analysis with tools such as SAMtools.\
`--mgi-convert-flag`: [FLAG] Turn on the MGI-to-Illumina header conversion module. Default: off.\
`--mgi-bc5`: Input an i5 barcode for Illumina header conversion. Default: 'PLACEHOLDERi5'.\
`--mgi-bc7`: Input an i7 barcode for Illumina header conversion. Default: 'PLACEHOLDERi7'.\
`--mgi-instrument`: Instrument name for Illumina header conversion. Default: 'PLACEHOLDERinstrument'.\
`--mgi-run`: Run ID for Illumina header conversion. Default: 'PLACEHOLDERrun'.\

### Advanced options
Further options that can be specified to alter the behavior of ReadZor.\
`--threads, -t`: Number of threads to use. Default: platform-dependent through auto-detection (assigned CPUs on Slurm-managed systems, all-1 otherwise. Fallback: 1).\
`--min-raw-read-length`: Minimum length a raw (untrimmed) read must have to be considered valid. Default: 0.\
`--reads-for-phred-offset`: Number of reads to sample per file for detection of Phred quality encoding offset. Default: 500.\
`--chunk-size`: Number of reads per chunk sent to each worker. Default: platform-dependent (20,000 for Slurm-managed systems, 1000 otherwise).\
`--phred-offset`: Define Phred offset for all FASTQ files. Default: off (auto-detection per file).\

## Software version
Although Readzor can probably work on older versions of its dependencies, it has been developed for best performance using the following versions:
* Python: >= 3.14.7
* NumPy: >= 2.5.2
* Python-isal: >= 1.8.0

## Issues and bug reports
Please leave a message in [issues](https://github.com/abjanssen/ReadZor/issues) or [discussions](https://github.com/abjanssen/ReadZor/discussions) if you notice an issue, bug, or otherwise. 

## License
This project is provded under the GNU General Public License v3.0 (GPLv3).

## Author
Axel B. Janssen [Google Scholar](https://scholar.google.com/citations?user=TWi-ysEAAAAJ&hl) [GitHub](https://github.com/abjanssen/)

## Reference
A reference for Readzor will be available soon.

## Legal
The oligonucleotide sequences used for adapter trimming, included in this work, are copyrighted and protected by intellectual property, including issued or pending patents, copyright, and trade secrets.\
Oligonucleotide sequences © 2026 Illumina, Inc. All rights reserved.
