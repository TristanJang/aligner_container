import subprocess
import json
import os
import pytest

# Directory paths (adjust if needed)
TEST_CASE_DIR = "/data/test_cases/"
REFERENCE = "/data/Test_reference/hg38_chr22.fasta"

# Acceptance thresholds (REQ 04)
MIN_MAPPED_PERCENT = 95.0
MAX_PERCENT_DUPLICATES = 10.0


@pytest.mark.parametrize("case_number", [1, 2, 3, 4])
def test_aligner_output_and_metrics(case_number):
    """
    Executes aligner on specified test case, then verifies:
    - Output BAM and JSON files are created.
    - % mapped reads exceeds threshold.
    - Computational stats are present.
    """

    fq1 = os.path.join(TEST_CASE_DIR, f"test_case_{case_number}_R1.fastq.gz")
    fq2 = os.path.join(TEST_CASE_DIR, f"test_case_{case_number}_R2.fastq.gz")
    output_bam = f"/data/aligned_test_case_{case_number}.bam"
        # Check if outputs already exist
    output_json = output_bam.replace(".bam", "_results.json")

    if not (os.path.isfile(output_bam) and os.path.isfile(output_json)):
        # Run aligner.py only if output files are missing
        subprocess.run([
            "python", "/app/aligner/aligner.py",
            "--fq1", fq1,
            "--fq2", fq2,
            "--ref", REFERENCE,
            "--output", output_bam
        ], check=True)
    else:
        pass

    # Check if aligner output files exist
    assert os.path.isfile(output_bam), f"BAM file missing for test case {case_number}"
    assert os.path.isfile(output_json), f"Result JSON missing for test case {case_number}"

    # Load JSON results
    with open(output_json, "r") as f:
        results = json.load(f)

    # Verify alignment statistics (REQ 04)
    stats = results["alignment_stats"]
    assert stats["percent_mapped"] >= MIN_MAPPED_PERCENT, \
        f"Mapped reads below threshold in test case {case_number}"

    assert stats["percent_duplicates"] <= MAX_PERCENT_DUPLICATES, \
        f"Too many duplicate reads in test case {case_number}"

    # Verify computational stats exist (REQ 05/06)
    comp_stats = results["computational_stats"]
    assert "total_runtime_seconds" in comp_stats, "Missing runtime metric"
    assert "peak_memory_MB" in comp_stats, "Missing memory metric"
    assert "cpu_usage_percent" in comp_stats, "Missing CPU usage metric"

    print(f"Test case {case_number} PASSED.")
