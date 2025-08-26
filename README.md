# Aligner Pipeline (Docker container)

## Overview

This repository contains a modular, containerized alignment pipeline using Python, BWA, and Samtools.

While creating this project, I worked with Docker containers and automated testing using pytest—both of which were new tools for me. Successfully integrating them into a reproducible workflow was a valuable experience that I’m excited to build upon.

---

This repository contains a standalone system designed to:
- Align raw sequencing reads using `bwa mem`
- Verify outputs against predefined software requirements
- Produce a human-readable Markdown report summarizing results

The system is modular, consisting of:
- `aligner.py`: a simple wrapper around `bwa mem` for read alignment
- A testing scripts using `pytest`
- A reporting script to generate Markdown reports from test outcomes

This pipeline assumes the hg38 / GRCh38.p14 reference genome, downloaded from NCBI, and indexed using bwa index. It was also chosen following its mention in the assignment instructions and its common usage in human diagnostics. For testing purposes, I used chromosome 22 from hg38 as the reference genome to allow for faster iteration and demonstration of the pipeline. The system is designed to work equally with full hg38 but was constrained by time and local resources.


---

## Repository Structure
sophia-technical/
├── aligner/ # Main aligner script
├── tests/ # Automated tests using pytest
├── report/ # Reporting script and example outputs
├── Dockerfile #enables container build
├── requirements.txt # Python dependencies
├── README.md # Project documentation
└── acceptance_criteria.md # Defined software requirements and test criteria


---

## Usage Instructions

### 1. Build Docker Image

```bash
docker build --no-cache -t sophia_technical . 
#Mounting Test Case Data and References
docker run --rm -it \
    -v /path/to/data/:/data/ \
    sophia_technical /bin/bash
#Running Aligner Script
python aligner/aligner.py \
    --fq1 /data/test_cases/test_case_1_R1.fastq.gz \
    --fq2 /data/test_cases/test_case_1_R2.fastq.gz \
    --ref /data/reference_genome/hg38_reference.fasta \
    --output /data/output_test_case_1.bam
    #Results in JSON and aligned Bam for each test case

#Running Automated Testing
pytest tests/test_aligner.py
    #Output is /data/test_summary.json

#Generate report
python report/reporter.py
    #Output is /data/final_test_report.md


```
##Input Format
FASTQ Files: Paired-end or single-end (.fastq.gz)
Reference Genome: FASTA (must be indexed via bwa index)

##Output Files
BAM file

JSON file with:

Read counts

Mapping stats

Quality metrics

Computational resource stats

## Dependencies
Handled via Docker:
bwa, samtools, python >= 3.10, pysam, pytest

## Design Documentation
Script Logic
aligner/aligner.py
Runs BWA MEM on provided FASTQ files against a reference genome. Pipes output through Samtools to generate a BAM file. Collects:

Alignment statistics (samtools flagstat and pysam).

Resource usage (via time).

Outputs metrics as JSON.

tests/test_aligner.py
Runs multiple test cases, checking:

BAM/JSON file creation.

Presence and correctness of alignment statistics.

Resource metrics reporting.

Pass/Fail for defined Acceptance Criteria (REQ01–REQ06).
Results aggregated into a test summary JSON.

report/reporter.py
Reads the test summary and individual JSON results to auto-generate a Markdown report summarizing:

Test date.

Software requirements.

Acceptance criteria.

Per-test results in a readable table.

Processing Flow & Assumptions
Reads are assumed to be gzipped FASTQ files.

Reference genome must be pre-indexed using bwa index.

The aligner script handles both paired-end and single-end inputs.

JSON metrics are generated alongside each BAM output.

Resource tracking (RAM, CPU, runtime) is used to evaluate REQ05 and REQ06.

Input and output files are handled externally via Docker-mounted /data/ directory.

Key Design Decisions
Containerization:
A Docker container ensures reproducibility across systems, including all dependencies (BWA, Samtools, Python packages).

Modularity:
Separate scripts for alignment, testing, and reporting ensure clarity and maintainability.

Lightweight Implementation:
Simple Python wrappers around well-established bioinformatics tools avoid unnecessary complexity.

Automated Reporting:
A Markdown report is generated automatically to facilitate transparent documentation and review.

Simplicity First:
Given the scope and timeline, the pipeline focuses on being functional, reproducible, and clear rather than feature-heavy.

#Challenges & Learning

Docker volume mounting and container filesystem behavior.

Indexing large reference genomes efficiently (e.g., using chr22 for testing).

Handling potential failures (e.g., empty BAMs) gracefully in the test script.

These challenges gave me a better practical understanding of containerized workflows and automated testing in bioinformatics.

