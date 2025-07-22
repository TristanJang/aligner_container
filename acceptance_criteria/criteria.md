# Acceptance Criteria for Software Requirements

This document defines clear and testable acceptance criteria for each software requirement specified in the assignment instructions.

---

## REQ 01: The aligner wrapper should be able to align reads in pair-ended mode.

- **AC 01.1**: Given two paired-end FASTQ files (`R1` and `R2`) and an indexed reference genome, the aligner wrapper must execute successfully without errors.
- **AC 01.2**: The output BAM file must pass basic integrity checks (`samtools quickcheck`).
- **AC 01.3**: The BAM file must contain at least one properly paired read, as reported by `samtools flagstat`.

---

## REQ 02: The aligner wrapper should be able to align reads in single-ended mode.

- **AC 02.1**: Given a single FASTQ file and an indexed reference genome, the aligner wrapper must execute successfully without errors.
- **AC 02.2**: The output BAM file must pass integrity checks and contain aligned reads.

---

## REQ 03: The aligner wrapper should be able to align reads using different reference genomes.

- **AC 03.1**: The aligner must accept any valid, indexed reference genome (e.g., hg38 or alternative references).
- **AC 03.2**: Successful execution and valid BAM output must occur regardless of the specific reference genome used.

---

## REQ 04: The aligner wrapper should correctly report alignment statistics.

- **AC 04.1**: Total number of input reads must be accurately reported.
- **AC 04.2**: Percentage of mapped reads must be reported.
- **AC 04.3**: Percentage of unmapped reads must be reported.
- **AC 04.4**: If possible, percentage of duplicated reads must be reported.
- **AC 04.5**: Average base quality score (Phred score) should be reported (planned for future iteration).
- **AC 04.6**: Average mapping quality (MAPQ) must be reported.

*Note: For the current version, basic reporting via `samtools flagstat` fulfills the minimum requirement.*

---

## REQ 05: The aligner wrapper should not exceed pre-defined computational resources.

- **AC 05.1**: CPU usage should not exceed 400% when using 4 CPUs (can be monitored via `psutil`).
- **AC 05.2**: Memory usage must not exceed 16 GB when processing test cases.
- **AC 05.3**: CPU and memory usage should be reported alongside alignment stats (planned for future iteration).

---

## REQ 06: The aligner wrapper should run in less than 2 minutes with small input data (<10,000 reads).

- **AC 06.1**: Using a small input dataset, the aligner should complete execution (alignment + BAM output) within 2 minutes.
- **AC 06.2**: Execution time should be printed or logged as part of the run summary.

---

# Notes:
- Reporting and monitoring components are being progressively expanded.
- All PASS/FAIL outcomes will be reported in a structured Markdown report via `reporter.py`.

---

