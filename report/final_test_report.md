# Alignment Pipeline Test Report

**Date of Test:** 2025-07-22 10:50:12

---

## Software Requirements

- **REQ 01:** Paired-end read alignment.
- **REQ 02:** Single-end read alignment.
- **REQ 03:** Multi-reference genome alignment.
- **REQ 04:** Reporting of key alignment statistics.
- **REQ 05:** Resource usage reporting (CPU, RAM).
- **REQ 06:** Fast runtime on small datasets (<10,000 reads).

---

## Acceptance Criteria Overview

See acceptance_criteria.md for full details of each AC.


---

## Test Cases

- **Test Case 1–4:** Provided FASTQ files (R1/R2) tested against hg38 or chr22 references.

---

## Results Summary

| Test Case | AC_01_01 | AC_02_01 | AC_03_01 | AC_04_01 | AC_04_02 | AC_05_1 | AC_05_2 | AC_06 | Overall |

|-----------|---------|----------|----------|----------|----------|---------|---------|-------|---------|
| 1 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS |
| 2 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS |
| 3 | FAIL | FAIL | FAIL | None | None | None | FAIL | None | FAIL |
| 4 | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS |

---

## Detailed Results per Test Case

### Test Case 1

**Alignment Statistics:**
```
total_reads: 50372
mapped_reads: 9188
percent_mapped: 18.24
percent_unmapped: 81.76
percent_duplicates: 0.74
average_base_quality: 33.28
average_mapping_quality: 23.8
```
**Computational Statistics:**
```
total_runtime_seconds: 26.76
peak_memory_MB: 136.21484375
cpu_usage_percent: 77.0
```

---

### Test Case 2

**Alignment Statistics:**
```
total_reads: 10080
mapped_reads: 1978
percent_mapped: 19.62
percent_unmapped: 80.38
percent_duplicates: 0.79
average_base_quality: 33.36
average_mapping_quality: 23.89
```
**Computational Statistics:**
```
total_runtime_seconds: 7.31
peak_memory_MB: 97.98828125
cpu_usage_percent: 93.0
```

---

### Test Case 3

_No result JSON generated (likely alignment failure)._

---

### Test Case 4

**Alignment Statistics:**
```
total_reads: 1002747
mapped_reads: 149401
percent_mapped: 14.9
percent_unmapped: 85.1
percent_duplicates: 0.27
average_base_quality: 32.87
average_mapping_quality: 34.01
```
**Computational Statistics:**
```
total_runtime_seconds: 174.51
peak_memory_MB: 198.0546875
cpu_usage_percent: 100.0
```

---
