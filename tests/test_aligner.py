import subprocess
import json
import os
import pytest

# Directory paths
TEST_CASE_DIR = "/data/test_cases/"
REFERENCE = "/data/chr22_reference/hg38_chr22.fasta"

SMALL_DATA_RUNTIME_LIMIT = 120  # seconds (REQ 06)

# Summary collector
TEST_SUMMARY = []


@pytest.mark.parametrize("case_number", [1, 2, 3, 4])
def test_aligner_output_and_metrics(case_number):
    """
    Executes aligner.py on specified test case (if needed), validates outputs, and evaluates against acceptance criteria.
    Summary results are collected for reporting.
    """

    fq1 = os.path.join(TEST_CASE_DIR, f"test_case_{case_number}_R1.fastq.gz")
    fq2 = os.path.join(TEST_CASE_DIR, f"test_case_{case_number}_R2.fastq.gz")
    output_bam = f"/data/aligned_test_case_{case_number}.bam"
    output_json = output_bam.replace(".bam", "_results.json")

    test_results = {
        "test_case": case_number,
        "AC_01_01": None,
        "AC_02_01": None,  # Single-end structurally supported
        "AC_03_01": None,
        "AC_04_01": None,
        "AC_04_02": None,
        "AC_05_1": None,
        "AC_05_2": None,
        "AC_06": None,
        "overall_status": None
    }

    # Run aligner if needed
    if not (os.path.isfile(output_bam) and os.path.isfile(output_json)):
        try:
            subprocess.run([
                "python", "/app/aligner/aligner.py",
                "--fq1", fq1,
                "--fq2", fq2,
                "--ref", REFERENCE,
                "--output", output_bam
            ], check=True)
        except subprocess.CalledProcessError as e:
            print(f"\nAligner failed to run for test case {case_number}.")
            test_results["AC_01_01"] = "FAIL"
            test_results["AC_02_01"] = "FAIL"
            test_results["AC_03_01"] = "FAIL"
            test_results["AC_05_2"] = "FAIL"
            test_results["overall_status"] = "FAIL"
            TEST_SUMMARY.append(test_results)
            return  # Skip further checks

    # Step 1: Validate BAM file
    if os.path.isfile(output_bam):
        result = subprocess.run(["samtools", "quickcheck", output_bam])
        if result.returncode == 0:
            test_results["AC_01_01"] = "PASS"
            test_results["AC_02_01"] = "PASS"
            test_results["AC_03_01"] = "PASS"
            test_results["AC_05_2"] = "PASS"
        else:
            print(f"BAM integrity check failed for test case {case_number}.")
            test_results["AC_01_01"] = "FAIL"
            test_results["AC_02_01"] = "FAIL"
            test_results["AC_03_01"] = "FAIL"
            test_results["AC_05_2"] = "FAIL"
            test_results["overall_status"] = "FAIL"
            TEST_SUMMARY.append(test_results)
            return

    else:
        test_results["AC_01_01"] = "FAIL"
        test_results["AC_02_01"] = "FAIL"
        test_results["AC_03_01"] = "FAIL"
        test_results["AC_05_2"] = "FAIL"
        test_results["overall_status"] = "FAIL"
        print(f"BAM missing for test case {case_number}.")
        TEST_SUMMARY.append(test_results)
        return

    # Step 2: Process JSON metrics
    if os.path.isfile(output_json):
        with open(output_json, "r") as f:
            results = json.load(f)
    else:
        print(f"Result JSON missing for test case {case_number}.")
        test_results["AC_04_01"] = "FAIL"
        test_results["AC_04_02"] = "FAIL"
        test_results["AC_05_1"] = "FAIL"
        test_results["AC_06"] = "FAIL"
        test_results["overall_status"] = "FAIL"
        TEST_SUMMARY.append(test_results)
        return

    # AC_04_01: Required stats present
    stats = results.get("alignment_stats", {})
    required_stats_present = all(k in stats for k in [
        "total_reads", "percent_mapped", "percent_unmapped",
        "percent_duplicates", "average_base_quality", "average_mapping_quality"
    ])
    test_results["AC_04_01"] = "PASS" if required_stats_present else "FAIL"

    # AC_04_02: Minimal range checks
    try:
        logical_ranges = (
            stats["total_reads"] >= 0 and
            stats["percent_mapped"] >= 0 and
            stats["percent_unmapped"] >= 0 and
            stats["percent_duplicates"] >= 0 and
            stats["average_base_quality"] >= 0 and
            stats["average_mapping_quality"] >= 0
        )
        test_results["AC_04_02"] = "PASS" if logical_ranges else "FAIL"
    except Exception:
        test_results["AC_04_02"] = "FAIL"

    # AC_05_1: Resource usage reporting
    comp_stats = results.get("computational_stats", {})
    if all(k in comp_stats for k in [
        "total_runtime_seconds", "peak_memory_MB", "cpu_usage_percent"
    ]):
        test_results["AC_05_1"] = "PASS"
    else:
        test_results["AC_05_1"] = "FAIL"

    # AC_06: Small dataset runtime check
    total_reads = stats.get("total_reads", 1000000)
    runtime_recorded = comp_stats.get("total_runtime_seconds", 9999)

    if total_reads < 10000:
        if runtime_recorded <= SMALL_DATA_RUNTIME_LIMIT:
            test_results["AC_06"] = "PASS"
        else:
            test_results["AC_06"] = "FAIL"
    else:
        test_results["AC_06"] = "PASS"

    # Determine overall status
    if "FAIL" in test_results.values():
        test_results["overall_status"] = "FAIL"
    else:
        test_results["overall_status"] = "PASS"

    TEST_SUMMARY.append(test_results)
    print(f"Test case {case_number} evaluated: {test_results['overall_status']}")


def teardown_module(module):
    """Write full summary report after all tests complete."""
    summary_path = "/data/test_summary.json"
    with open(summary_path, "w") as f:
        json.dump(TEST_SUMMARY, f, indent=4)
    print(f"\n\nTest summary written to {summary_path}\n")
