import subprocess
import hashlib
from pathlib import Path
import gzip

DATA_DIR = Path(__file__).parent / "data"

EXPECTED_HASHES = {
    "results_summary.txt": "b1a12db3dbffa3e7b1a8e6757a2ab6fe68a4dd387cfceb361d549c2c9880d3e3",
    "test_paired_R1_paired_filtered.fastq.gz": "acf7526124d4987fe15634fc43862360b69418840b02b10309947b8e47d37e8a",
    "test_paired_R2_paired_filtered.fastq.gz": "89d8385c3544ed1ba01823cd572c21937be72c3ce570d8832035d77de4d29976",
    "test_paired_unpaired_filtered.fastq.gz": "f32092d0e15ee34b1d25a97cf3b98fbd8921acf9986192da6ae5b0beac625786",
    "test_unpaired_testhash3_filtered.fastq.gz": "1d4b0706fd1821ac3ec78927f61b69d40c0dce2c097c78de0a80f9f19d966a51",
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
    assert len(files_in_output) == 7, f"Expected 7 files in output folder, found: {len(files_in_output)}"
    
    for filepath in files_in_output:
        if filepath.name == "parameters.txt":
            continue
    
        open_fn = gzip.open if filepath.name.endswith(".gz") else open
    
        with open_fn(filepath, "rb") as file:
            file_hash = hashlib.file_digest(file, "sha256").hexdigest()
    
        expected = EXPECTED_HASHES.get(filepath.name)
    
        assert expected is not None, f"Unexpected file found: {filepath.name}"
        assert file_hash == expected, f"File {filepath.name} incorrectly processed. Got: {file_hash}"