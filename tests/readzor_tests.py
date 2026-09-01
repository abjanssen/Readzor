import subprocess
import hashlib
from pathlib import Path
import gzip

DATA_DIR = Path(__file__).parent / "data"

EXPECTED_HASHES = {
    "results_summary.txt": "f6b13b02882a2e109bc90e7419d2066e53efc34e05e9a350e7b8678b9e7017e8",
    "test_paired_R1_paired_filtered.fastq.gz": "b9def62118ac494217c3ca3df30818f90788973c9722fe0fa688ff37039c46dc",
    "test_paired_R2_paired_filtered.fastq.gz": "76a904ecd486f981d5fa4cf25512c40961e8fd0a9f88a3567a5956d5b90bc007",
    "test_paired_unpaired_filtered.fastq.gz": "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855",
    "test_unpaired_testhash3_filtered.fastq.gz": "6f941cea0eade596e16a5c4474fda8cfd9e4f340c33ca6da07534355b18addf1",
}

def test_readzor(tmp_path):
    # Map the exact output filenames to your expected hashes

    result = subprocess.run(
        ["readzor", "-h"],
        cwd=tmp_path,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    assert result.returncode == 0, f"Help message failed: {result.stderr}"
    
    input_files = [str(f) for f in DATA_DIR.glob("*.fastq*")]
    assert len(input_files) > 0, f"No FASTQ files found in {DATA_DIR}"

    result = subprocess.run(
        ["readzor", "--input-files", *input_files,
         "--gzip", "--nucl-filter",
         "--cut-flag", "--endqual-filter-flag",
         "--n-trimming-flag", "--slider-filter-flag",
         "--poly-filter-flag", "--adapter-trim-flag",
         "--kmer-filter-flag"],
        cwd=tmp_path,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True
    )
    
    assert result.returncode == 0, f"Readzor failed to execute: {result.stderr}"

    subdirs = [d for d in tmp_path.iterdir() if d.is_dir()]
    assert len(subdirs) == 1, f"Expected exactly one output folder, found: {len(subdirs)}"

    output_dir = subdirs[0]
    files_in_output = list(output_dir.iterdir())
    assert len(files_in_output) == 6, f"Expected 6 files in output folder, found: {len(files_in_output)}"
    
    for filepath in files_in_output:
        if filepath.name == "Readzor_log.txt":
            continue
    
        open_fn = gzip.open if filepath.name.endswith(".gz") else open
    
        with open_fn(filepath, "rb") as file:
            file_hash = hashlib.file_digest(file, "sha256").hexdigest()
    
        expected = EXPECTED_HASHES.get(filepath.name)
    
        assert expected is not None, f"Unexpected file found: {filepath.name}"
        assert file_hash == expected, f"File {filepath.name} incorrectly processed. Got: {file_hash}"