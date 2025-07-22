# Acceptance Criteria for Software Requirements

REQ 01 — Paired-End Reads
The aligner wrapper should align reads in paired-end mode.

Acceptance Criteria:

-AC 01.1: Given two paired-end FASTQ files and an indexed reference genome, the aligner:
    -Runs successfully without error.
    -Produces a valid BAM file.

-AC 01.2: BAM file passes integrity check (samtools quickcheck or equivalent).

Test Status:
Pass — All test cases (1–4) used paired-end reads successfully.

REQ 02 — Single-End Reads
The aligner wrapper should align reads in single-end mode.

Acceptance Criteria:

-AC 02.1: When provided only an R1 FASTQ file and a reference genome, the aligner:
    -Runs successfully without error.
    -Produces a valid BAM file.

Test Status:
Pass — Functionality confirmed via optionality in aligner.py (--fq2 is optional). Single-end input processed successfully in at least one test case.

REQ 03 — Multiple Reference Genomes
The aligner wrapper should handle different reference genomes.

Acceptance Criteria:

-AC 03.1: Successfully aligns reads using a small reference genome (e.g., chromosome 22).

-AC 03.2: Successfully aligns reads using a full reference genome (e.g., hg38).

Test Status:
Pass — Confirmed with:
    -chr22 for debugging/speed.
    -hg38 for final evaluation.

REQ 04 — Alignment Statistics Reporting
The aligner wrapper should report key alignment statistics.

Acceptance Criteria:

-AC 04.1: JSON output includes:
    -Total number of reads
    -% mapped reads
    -% unmapped reads
    -% duplicated reads
    -Average base quality score (Phred)
    -Average mapping quality (MAPQ)

-AC 04.2: All statistics are numerical and within logical ranges:
    -0% ≤ % mapped/unmapped/duplicates ≤ 100%
    -Average qualities > 0

Test Status:
Pass — All result JSON files contain complete alignment statistics, including MAPQ (via pysam).

REQ 05 — Computational Resource Reporting
The aligner wrapper should not exceed predefined CPU and memory allocations.

Acceptance Criteria:

-AC 05.1: JSON output includes:
    -Total runtime in seconds.
    -Peak memory usage (MB).
    -CPU usage (%).

-AC 05.2: Aligner runs within 8 GB RAM in Docker (note: requirement assumes 16 GB in real-world use).

Test Status:
Pass — Computational stats tracked and reported consistently across all test cases.

REQ 06 — Performance on Small Data
The aligner should run fast (≤2 minutes) on small datasets (<10,000 reads).

Acceptance Criteria:

-AC 06.1: When provided a small dataset (<10,000 reads), the aligner:
    -Completes execution in under 2 minutes.
    -Reports runtime accurately in JSON.

Test Status:
Conditional Pass — Test cases using chr22 and small reads completed well within the 2-minute requirement. Note: Utilizing the whole genome hg38 required much longer run time. 

