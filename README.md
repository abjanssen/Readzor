
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
![Yes, it's in Python'](https://img.shields.io/badge/Language-Python3-steelblue.svg)

# Welcome to ReadZor 

ReadZor is a fast, modular, and fully-featured FASTQ quality trimming and filtering pipeline. 
By default, ReadZor takes a modular approach: <ins>all processing modules are off by default</ins>, leaving you the explicit ability to adapt the workflow according to your needs and wishes.

## Feature overview
*   **Modular trimming approach**: Enable only the filters you need:
    *   Quality dependent end trimming
    *   Sliding window quality trimming
    *   Homopolymer trimming from both read ends
    *   Adapter trimming for TruSeq and Nextera adapters.
    *   N-base trimming from read ends or in general
    *   Low complexity filtering based on K-mer counts.
    *   Set-length read end trimming.
    
By design, each filter will indepedently assess each read. Outcomes are merged per end, and the most stringent trimming limit is applied.

For example, if the quality trimmer determines 5 bases should be removed from the 3' end, but the adapter trimmer identifies 12 bases to remove from that same end, ReadZor will merge these outcomes and trim the most restrictive amount—12 bases—to ensure the highest quality output.

*   **Auto-detection**:
    *   File Pairing: Identifies paired vs. unpaired FASTQ files using internal header information, independent of file names.
    *   Read Numbering: Determines read numbers directly from headers.
    *   Compression: Detects GZIP compression via magic bytes, regardless of the file extension.
    *   Phred Offsets: Detects per-file Phred quality offsets (i.e., Phred33 vs. Phred64) based on quality string symbols.

Using auto-detection methods, ReadZor provides an easy-to-use platform for novice users, enabling high quality read trimming with minimal inputs.

*   **Extra features**:
    *   Fully automatic mode: invoking ReadZor in fully automatic mode (optional file specification possible), will use predetermined settings to optimize read processing. For more info, invoke ReadZor with flags --full-auto and --help.
    *   MGI/BGI header conversion: Built-in module to convert MGI/BGI fastq headers into standard Illumina format for downstream compatibility (e.g., with SAMtools).
    *   Safe and traceable outputs: Automatically generates isolated, timestamped output directories to prevent accidental file overwrites. Each run includes a comprehensive summary file detailing the exact parameters used and the final trimming outcomes for full reproducibility.
    *   HPC integration: Automatically reads Slurm environment variables to seamlessly scale threads and optimize performance on HPC clusters.
    *   Live progress monitoring: Features a dynamic progress bar that provides accurate real-time feedback on processing speed, percentage completed, and estimated time remaining.
    *   Real-time record validation: ReadZor detects malformed or corrupted FASTQ records on-the-fly, ensuring high-quality output.
    *   Gzip compression: Optionally compress your output files and control the compression depth (1-9) to optimize storage size.

Equipped with a number of other features, like built-in validation and seamless Slurm detection methods, ReadZor delivers a robust execution environment, ensuring reproducible outputs at any scale.

## Installation

### Test
* Type `readzor` and it should output the help information.
* Type `readzor -v` and you should see an output like `ReadZor version: X.X.X`

## Quick Start
```bash
python readzor.py -GO
```

```bash
python readzor.py -GO --input-files /path/to/input/files
```

## Command line options 

## Detailed usage
### Invoking modules
All modules are <ins>off</ins> by default. Use the respective `--[module]-flag` to turn them on using their default settings.
Further function parameters adapt behaviour as follow:


### General quality filters
Filter reads based on their final length or average quality *after* trimming.
*   `--min-average-qual <int>`: Minimum average quality of the output read (default: off).
*   `--min-length <int>`: Minimum length of the output read (default: off).
*   `--max-length <int>`: Maximum length of the output read (default: off).
*   `--nucl-filter` [FLAG]: Reject reads containing `N` bases anywhere (default: off).


### Sliding window quality trimming
Trims reads based on a sliding window. Retains the longest portion of the read that passes, and the highest average window quality in the event of an even read split.
*   `--slider-filter-flag` [FLAG]: Turn on the sliding window quality trimming module (default: off).
*   `--slider-window <int>`: Use a sliding window of N size (default: 5)
*   `--slider-quality <int>`: Average quality cutoff for a window to pass (default: 20)
*   `--slider-step <int>`: Step size between sliding window locations (default: 1)


### Homopolymer nucleotide trimming
Trims reads for homopolymer runs of defined nucleotides. Useful for Illumina two-color chemistry (NovaSeq/NextSeq) where lack of signal is read as a 'G'. Trims homopolymer runs from the ends. Multiple basese are analysed for independently (i.e. no heteropolymer runs are found).
*   `--poly-filter-flag` [FLAG]: Turn on the homopolymer nucleotide trimming module (default: off).
*   `--poly-bases-start <str>`: Bases to analyse for homopolymer run, at the start of thread (default: None). Overrides `--poly-bases-both` for the start of the read.
*   `--poly-bases-end <str>`: Bases to analyse for homopolymer run, at the end of thread (default: "G"). Overrides `--poly-bases-both` for the end of the read.
*   `--poly-bases-both <str>`: Bases to analyse for homopolymer run, at both ends of the read (default: None). Overridden by `--poly-bases-start` and `--poly-bases-end` for their respective ends.
*   `--poly-length-start <int>`: Threshold value for homopolymer run, at the start of thread (default: off). Overrides `--poly-bases-both` for the start of the read.
*   `--poly-length-end <int>`: Threshold value for homopolymer run, at the end of thread (default: 10). Overrides `--poly-bases-both` for the end of the read.
*   `--poly-length-both <int>`: Threshold value for homopolymer run, at both ends of the read (default: off).  Overridden by `--poly-bases-start` and `--poly-bases-end` for their respective ends.


### Set-Length end trimming
Trim a set number of bases from the ends of each read, regardless of quality or sequence.
*   `--cut-flag` [FLAG]: Turn on the set-length end trimming module (default: off).
*   `--cut-start <int>`: Number of bases to trim of start of the read (default: 0). Overrides `--cut-both` for the start of the read
*   `--cut-end <int>`: Number of bases to trim of end of the read (default: 0). Overrides `--cut-both` for the end of the read
*   `--cut-both <int>`: Number of bases to trim of end of the read (default: 0). Overridden by `--cut-start` and `--cut-end` for their respective ends.

### N-Nucleotide end trimming
Preconfigured wrapper for homopolymer nucleotide trimming that trims of any N nucleotides of the ends of the reads. This becomes redundant when `--nucl-filter` is set.
*   `--n-end-trimming-flag` [FLAG]: Turn on the N nucleotide nucleotide trimming module (default: off).

### Low complexity filtering ()
Filter reads that fall below a certain complexity threshold, calculated via unique k-mer frequencies.
*   `--kmer-filter-flag` [FLAG]: Turn on k-mer based low complexity filtering.
*   `--kmer-size <int>`: K-mer size to use in analysis (default: 4). Multiple comma-seperated values may be provided, which will be analysed indepedently.
*   `--kmer-cutoff <int>`: Minimum percentage of unique k-mers present in the read, relative to the maximum possible, required to pass the complexity filter. Lower values are stricter (default: 50).




## Citation

## Issues and bug reports
Please leave a message in [issues](https://github.com/abjanssen/ReadZor/issues) or [discussions](https://github.com/abjanssen/ReadZor/discussions) if you notice an issue, bug, or otherwise. 

## License
This project is provded under the GNU General Public License v3.0 (GPLv3).

## Author
[Axel B. Janssen](https://scholar.google.com/citations?user=TWi-ysEAAAAJ&hl)


















## Installation

ReadZor requires **Python 3** and **NumPy**.

```bash
# Clone the repository
git clone https://github.com/abjanssen/ReadZor.git
cd ReadZor

# Install dependencies
pip install numpy

# Make the script executable
chmod +x readzor.py
```

## Quick Start

### 1. Fully Automatic Mode (-GO / --full-auto)
The simplest way to run ReadZor. 
In this mode will either take input files specified from the user, or auto-detect all FASTQ files in your current directory.
All other parameters will be ignored.
ReadZor will determine paired or unpaired files, and process them accordingly.
Through the fully automatic mode, reads will only be processed through quality-dependent ends trimming, and adapter trimming. In addition read containing N bases will be filtered out.



```bash
python readzor.py -GO --input-files /path/to/files/
python readzor.py -GO --input-unpaired /path/to/unpaired/files --input-paired /path/to/paired/files
```

### 2. Manual Input Specification
You can manually specify unpaired or paired files.

```bash
# Unpaired reads
python readzor.py --input-unpaired sample1.fastq sample2.fastq.gz --endsquality-filter-flag

# Paired-end reads (Provide R1 and R2 consecutively)
python readzor.py --input-paired sample_R1.fastq sample_R2.fastq --adapter-trim-flag --poly-filter-flag
```

## Modules and Usage




### 2. Quality-Dependent End Trimming (`--endsquality-filter-flag`)
Trim the ends of reads inward until a base meets the quality threshold.
*   `--min-quality-start <int>`, `--min-quality-end <int>`, `--min-quality-both <int>`







### 6. Adapter Trimming (`--adapter-trim-flag`)
Trims Illumina standard adapters (TruSeq3, Nextera, small RNA). 
*   `--adapter-fasta <path>`: Provide a custom FASTA file with additional adapter sequences.


### 8. MGI Header Conversion (`--mgi-convert-flag`)
Converts MGI (BGI) style headers (`<flowcell>L<lane>...`) to standard Illumina headers. Original headers are preserved in the plus-line.
*   Options: `--mgi-bc5`, `--mgi-bc7`, `--mgi-instrument`, `--mgi-run`

## Output

ReadZor automatically creates an isolated, timestamped output directory to prevent accidental overwrites:
`ReadZor_results_<YYYY-MM-DD_HH-MM-SS>/`

To specify the parent directory for this folder, use `--output /path/to/dir`. 

By default, output files are uncompressed. Use the `--gzip` flag to compress output files, and `--gzip-level` (1-9) to control compression depth.
