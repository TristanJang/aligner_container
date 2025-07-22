import json
import datetime
import os

# File paths
TEST_SUMMARY_JSON = "/data/test_summary.json"
RESULT_FILES_DIR = "/data/"
REPORT_PATH = "/data/final_test_report.md"

# Load test summary
with open(TEST_SUMMARY_JSON, "r") as f:
    summary_data = json.load(f)

# Build Markdown report
lines = []

# Report header
lines.append("# Alignment Pipeline Test Report\n")
lines.append(f"**Date of Test:** {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
lines.append("---\n")

# Software Requirements Overview
lines.append("## Software Requirements\n")
lines.append("- **REQ 01:** Paired-end read alignment.")
lines.append("- **REQ 02:** Single-end read alignment.")
lines.append("- **REQ 03:** Multi-reference genome alignment.")
lines.append("- **REQ 04:** Reporting of key alignment statistics.")
lines.append("- **REQ 05:** Resource usage reporting (CPU, RAM).")
lines.append("- **REQ 06:** Fast runtime on small datasets (<10,000 reads).")
lines.append("\n---\n")

# Acceptance Criteria
lines.append("## Acceptance Criteria Overview\n")
lines.append("See acceptance_criteria.md for full details of each AC.\n")
lines.append("\n---\n")

# Description of Test Cases
lines.append("## Test Cases\n")
lines.append("- **Test Case 1–4:** Provided FASTQ files (R1/R2) tested against hg38 or chr22 references.")
lines.append("\n---\n")

# Results Table
lines.append("## Results Summary\n")
lines.append("| Test Case | AC_01_01 | AC_02_01 | AC_03_01 | AC_04_01 | AC_04_02 | AC_05_1 | AC_05_2 | AC_06 | Overall |\n")
lines.append("|-----------|---------|----------|----------|----------|----------|---------|---------|-------|---------|")

for case in summary_data:
    row = f"| {case['test_case']} | {case['AC_01_01']} | {case['AC_02_01']} | {case['AC_03_01']} | {case['AC_04_01']} | {case['AC_04_02']} | {case['AC_05_1']} | {case['AC_05_2']} | {case['AC_06']} | {case['overall_status']} |"
    lines.append(row)

lines.append("\n---\n")

# Per-Test Metrics Summary
lines.append("## Detailed Results per Test Case\n")

for case in summary_data:
    case_number = case["test_case"]
    lines.append(f"### Test Case {case_number}\n")
    result_file = f"{RESULT_FILES_DIR}aligned_test_case_{case_number}_results.json"

    if os.path.isfile(result_file):
        with open(result_file, "r") as f:
            metrics = json.load(f)

        lines.append("**Alignment Statistics:**")
        stats = metrics.get("alignment_stats", {})
        lines.append("```")
        for k, v in stats.items():
            lines.append(f"{k}: {v}")
        lines.append("```")

        lines.append("**Computational Statistics:**")
        comp = metrics.get("computational_stats", {})
        lines.append("```")
        for k, v in comp.items():
            lines.append(f"{k}: {v}")
        lines.append("```")
    else:
        lines.append("_No result JSON generated (likely alignment failure)._")

    lines.append("\n---\n")

# Write report to Markdown file
with open(REPORT_PATH, "w") as f:
    f.write("\n".join(lines))

print(f"Markdown report generated: {REPORT_PATH}")
