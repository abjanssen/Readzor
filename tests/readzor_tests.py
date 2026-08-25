import subprocess
import hashlib
from pathlib import Path

DATA_DIR = Path(__file__).parent / "data"

EXPECTED_HASHES = {
    "parameters.txt": "c8e7bbaebb082458f12a0a7cad659c9a13d1a0f4d9b6db2481c665ffce7aa810",
    "results_summary.txt": "b1a12db3dbffa3e7b1a8e6757a2ab6fe68a4dd387cfceb361d549c2c9880d3e3",
    "test_paired_R1_paired_filtered.fastq.gz": "6c35bb8979940e18cdac535582319caffa73b46e28ed336c8802878d293108f8",
    "test_paired_R2_paired_filtered.fastq.gz": "e05fccaefbe5df4f70479e1f452b2b69c7cca42baa85965ed3f0bdcd76c6c754",
    "test_paired_unpaired_filtered.fastq.gz": "4228bc5e7ab669bcbb050ea5601cb0908542679a215e4316058a19574437355e",
    "test_unpaired_testhash3_filtered.fastq.gz": "ecba86c6016fc30c0496accd795397e98b80bb3ca8b9445f7529b64c6bc3faca",
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
        with open(filepath, "rb") as file:
            file_hash = hashlib.file_digest(file, "sha256").hexdigest()
            
            expected = EXPECTED_HASHES.get(filepath.name)
            
            assert expected is not None, f"Unexpected file found: {filepath.name}"
            assert file_hash == expected, f"File {filepath.name} incorrectly processed. Got: {file_hash}"