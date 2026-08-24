import subprocess
import hashlib
from pathlib import Path

DATA_DIR = Path(__file__).parent / "data"

EXPECTED_HASHES = {
    "output_file_1.txt": "your_sha256_hash_here",
    "output_file_2.txt": "another_sha256_hash_here",
    "output_file_3.txt": "another_sha256_hash_here",
    "output_file_4.txt": "another_sha256_hash_here",
    "output_file_5.txt": "another_sha256_hash_here",
    "output_file_6.txt": "another_sha256_hash_here",
    "output_file_7.txt": "another_sha256_hash_here"
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
        ["readzor", "--input-files", *input_files, "--gzip", "--nucl-filter","--cut-flag","--endqual-filter-flag", "--n-end-trimming-flag","--slider-filter-flag","--poly-filter-flag", "--adapter-trim-flag", "--kmer-filter-flag", "--mgi-convert-flag"],
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
    assert len(files_in_output) == 7, f"Expected 7 files in output folder, found: {len(files_in_output)}"
    
    for filepath in files_in_output:
        with open(filepath, "rb") as file:
            file_hash = hashlib.file_digest(file, "sha256").hexdigest()
            
            expected = EXPECTED_HASHES.get(filepath.name)
            
            assert expected is not None, f"Unexpected file found: {filepath.name}"
            assert file_hash == expected, f"File {filepath.name} incorrectly processed. Got: {file_hash}"