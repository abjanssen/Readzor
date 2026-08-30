import subprocess
import hashlib
from pathlib import Path
import gzip

DATA_DIR = Path(__file__).parent / "data"

EXPECTED_HASHES = {
    "results_summary.txt": "f6b13b02882a2e109bc90e7419d2066e53efc34e05e9a350e7b8678b9e7017e8",
    "test_paired_R1_paired_filtered.fastq.gz": "d3f327fa2442d49c87215a9c5bef50c71ca8141c6fb885961c5b9045e53c0c49",
    "test_paired_R2_paired_filtered.fastq.gz": "2b678c0851044e66fe7a8c04546bc2af468cf2d7f6544c7e3e243b072228e875",
    "test_paired_unpaired_filtered.fastq.gz": "afc3ac1631bc1f2f53116cb869b1fa23805d2f38bb81e89a1b0e1fe1ae938b0a",
    "test_unpaired_testhash3_filtered.fastq.gz": "d4d127cab35a4570456950e25b5384d92f970c11b50ee0f904a71a9e65c009ca",
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