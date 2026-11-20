"""Readzor-wide constants and shared run-time state.

Constants are imported by name (from readzor.constants import DEFAULT_ADAPTERS).
Run-time state is always accessed through the module (constants.GZIP_DETECTION),
because ACTIVE_PROGRESS_TRACKER and WORKER_PARAMETERS are reassigned during a run
and a by-name import would keep the old value.
"""

VERSION = "0.7.2"
DEFAULT_ADAPTERS = [
    ["AmpliSeq", [
        ["AmpliSeq", "CTGTCTCTTATACACATCT"]  #https://support-docs.illumina.com/SHARE/AdapterSequences/Content/AmpliSeq-Sequences.htm
    ]],
    ["Illumina_miRNA", [
        ["Illumina_miRNA", "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"]     #https://support-docs.illumina.com/SHARE/AdapterSequences/Content/Illumina-miRNA-Indexes.htm
    ]],
    ["Illumina_RNA", [
        ["Illumina_RNA", "ACTGTCTCTTATACACATCT"]    #https://support-docs.illumina.com/SHARE/AdapterSequences/Content/Nextera_Illumina-Sequences.htm
    ]],
    ["Nextera", [
        ["Nextera", "CTGTCTCTTATACACATCT"]  #https://support-docs.illumina.com/SHARE/AdapterSequences/Content/Nextera_Illumina-Sequences.htm
    ]],
    ["TruSeq", [
        ["Read 1", "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"],  #https://support-docs.illumina.com/SHARE/AdapterSequences/Content/UDIndexes.htm
        ["Read 2", "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"]  #https://support-docs.illumina.com/SHARE/AdapterSequences/Content/UDIndexes.htm
    ]],
    ["TruSeq_DNA_methylation", [
        ["Read 1", "AGATCGGAAGAGCACACGTCTGAAC"],  #https://support-docs.illumina.com/SHARE/AdapterSequences/Content/TruSeq-DNAMethyl.htm
        ["Read 2", "AGATCGGAAGAGCGTCGTGTAGGGA"]  #https://support-docs.illumina.com/SHARE/AdapterSequences/Content/TruSeq-DNAMethyl.htm
    ]],
    ["TruSeq_Ribo_profile", [
        ["TruSeq_Ribo_profile", "AGATCGGAAGAGCACACGTCT"]  #https://support-docs.illumina.com/SHARE/AdapterSequences/Content/TruSeq-RiboProfile.htm
    ]],
    ["TruSeq_small_RNA", [
        ["TruSeq_small_RNA", "TGGAATTCTCGGGTGCCAAGG"]    #https://support-docs.illumina.com/SHARE/AdapterSequences/Content/TruSeq-SmallRNA.htm
    ]]
]
FULL_AUTO_PRESERVED_DESTS = {"input_files", "input_paired", "input_unpaired", "input_interleaved", "full_auto"}
FULL_AUTO_OVERRIDES = {
    "endqual_filter_flag": True,
    "endqual_min_both": 20,
    "adapter_filter_flag": True,
    "adapter_group": ["TruSeq"],
    "adapter_seed": 6,
    "adapter_mismatch": 2,
    "n_filter": True,
    "min_length_output_perc": 33,
    "gzip": True,
    "progress": True
}
PHRED_ALLOWED = bytes(range(33, 127))
NUCL_ATCG = b"ATCG"
NUCL_ATCGN = b"ATCGN"
COMP_TABLE = bytes.maketrans(NUCL_ATCGN, b"TAGCN")
PHRED64_TO_33 = bytes.maketrans(
    bytes(range(59, 127)),
    bytes(max(33, b - 31) for b in range(59, 127))
)
PHRED33_TO_64 = bytes.maketrans(
    bytes(range(33, 127)),
    bytes(min(126, b + 31) for b in range(33, 127))
)

GZIP_DETECTION = {}
STDIN_TEMP_FILES = []
ESTIMATED_READ_COUNTS = {}
ESTIMATED_BYTE_PER_READ = {}
ESTIMATED_ZIP_RATIO = {}
ACTIVE_PROGRESS_TRACKER = None
WORKER_PARAMETERS = None
