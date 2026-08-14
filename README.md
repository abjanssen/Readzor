
[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Yes, it's in Python'](https://img.shields.io/badge/Language-Python3-steelblue.svg)

# Welcome to ReadZor 

ReadZor is a fast, modular, and fully-featured FASTQ quality trimming and filtering pipeline. 
By default, ReadZor takes a modular approach: **all processing modules are off by default**, leaving you the explicit ability to adapt the workflow according to your needs and wishes.

## Feature overview
*   **Modular trimming approach**: Enable only the filters you need:
    *   Quality dependent end trimming
    *   Sliding window quality trimming
    *   Homopolymer trimming from read ends
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
* Type `readzor` and it should output its help screen.
* Type `readzor -v` and you should see an output like `ReadZor vX.X.X`

## Quick Start







## Command line options 

## Detailed usage
### Invoking modules
All modules are **off** by default. Use the respective `--[module]-flag` to turn them on using their default settings.
Further function parameters adapt module behaviour.

### General quality filters
Filter reads based on their final length or average quality *after* trimming.
*   `--min-average-qual <int>`: Minimum average quality of the output read (Default: 0).
*   `--min-length <int>`: Minimum length of the output read (Default: 1/3 of raw read length).
*   `--max-length <int>`: Maximum length of the output read.
*   `--nucl-filter`: Reject reads containing `N` bases anywhere.

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
git clone https://github.com/yourusername/ReadZor.git
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
python readzor.py -GO
python readzor.py --full-auto
```

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



### 1. Set-Length End Trimming (`--cut-flag`)
Trim a set number of bases from the ends of each read, regardless of quality.
*   `--cut-start <int>`, `--cut-end <int>`, `--cut-both <int>`

### 2. Quality-Dependent End Trimming (`--endsquality-filter-flag`)
Trim the ends of reads inward until a base meets the quality threshold.
*   `--min-quality-start <int>`, `--min-quality-end <int>`, `--min-quality-both <int>`

### 3. N-Nucleotide End Trimming (`--n-end-trimming-flag`)
Removes leading and trailing `N` bases from reads.

### 4. Sliding Window Quality Trimming (`--slider-filter-flag`)
Trims reads based on a sliding window. Retains the longest portion of the read where average window quality stays above the threshold.
*   `--slider-window <int>` (Default: 5)
*   `--slider-quality <int>` (Default: 20)
*   `--slider-step <int>` (Default: 1)

### 5. Homopolymer Nucleotide Trimming (`--poly-filter-flag`)
Useful for Illumina two-color chemistry (NovaSeq/NextSeq) where lack of signal is read as a 'G'. Trims homopolymer runs from the ends.
*   `--poly-bases-start <str>`, `--poly-bases-end <str>` (Default: "G"), `--poly-bases-both <str>`
*   `--poly-length-start <int>`, `--poly-length-end <int>` (Default: 10), `--poly-length-both <int>`

### 6. Adapter Trimming (`--adapter-trim-flag`)
Trims Illumina standard adapters (TruSeq3, Nextera, small RNA). 
*   `--adapter-fasta <path>`: Provide a custom FASTA file with additional adapter sequences.

### 7. Low Complexity Filtering (`--kmer-filter-flag`)
Discards reads that fall below a certain complexity threshold, calculated via unique k-mer frequencies.
*   `--kmer-size <int>` (Default: 4)
*   `--kmer-cutoff <int>`: Minimum % of unique k-mers required (Default: 50).

### 8. MGI Header Conversion (`--mgi-convert-flag`)
Converts MGI (BGI) style headers (`<flowcell>L<lane>...`) to standard Illumina headers. Original headers are preserved in the plus-line.
*   Options: `--mgi-bc5`, `--mgi-bc7`, `--mgi-instrument`, `--mgi-run`

## Output

ReadZor automatically creates an isolated, timestamped output directory to prevent accidental overwrites:
`ReadZor_results_<YYYY-MM-DD_HH-MM-SS>/`

To specify the parent directory for this folder, use `--output /path/to/dir`. 

By default, output files are uncompressed. Use the `--gzip` flag to compress output files, and `--gzip-level` (1-9) to control compression depth.
