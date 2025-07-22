# SOPHiA GENETICS Technical Assignment

## Project Overview

This repository contains a standalone system designed to:
- Align raw sequencing reads using `bwa mem`
- Verify outputs against predefined software requirements
- Produce a human-readable Markdown report summarizing results

The system is modular, consisting of:
- `aligner.py`: a simple wrapper around `bwa mem` for read alignment
- Testing scripts using `pytest`
- A reporting script to generate Markdown reports from test outcomes

This pipeline assumes the hg38 / GRCh38.p14 reference genome, downloaded from NCBI, and indexed using bwa index. It was also chosen following its mention in the assignment instructions and its common usage in human diagnostics. The indexed reference should be mounted into the container at runtime for alignments.


---

## Repository Structure

