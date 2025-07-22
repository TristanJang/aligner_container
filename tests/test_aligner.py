import json
import glob
import pytest

# === Acceptance Thresholds ===
PERCENT_MAPPED_THRESHOLD = 10.0  # Example: >= 10% mapped reads
RUNTIME_THRESHOLD_SECONDS = 600  # Example: ≤ 10 minutes
PEAK_MEMORY_MB_THRESHOLD = 8000  # Example: ≤ 8 GB

# === Load JSON Output ===
def load_results(json_path):
    with open(json_path, "r") as f:
        return json.load(f)

# === Collect All Result Files Automatically ===
# Adjust this glob as needed based on your folder structure:
result_files = glob.glob("*.json")

@pytest.mark.parametrize("json_file", result_files)
def test_alignment_requirements(json_file):
    results = load_results(json_file)
    stats = results["alignment_stats"]
    comp_stats = results["computational_stats"]

    # === Alignment Quality Tests ===
    assert stats["percent_mapped"] >= PERCENT_MAPPED_THRESHOLD, \
        f"{json_file}: Percent mapped too low: {stats['percent_mapped']}%"

    # === Computational Performance Tests ===
    assert comp_stats["total_runtime_seconds"] <= RUNTIME_THRESHOLD_SECONDS, \
        f"{json_file}: Runtime too long: {comp_stats['total_runtime_seconds']} seconds"

    assert comp_stats["peak_memory_MB"] <= PEAK_MEMORY_MB_THRESHOLD, \
        f"{json_file}: Peak memory usage too high: {comp_stats['peak_memory_MB']} MB"
