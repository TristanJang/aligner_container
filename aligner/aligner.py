#!/usr/bin/env python3

import argparse
import subprocess

#define alignment command 
def run_bwa_mem(fq1, fq2, reference, output_bam):
    """Basic bwa mem wrapper that outputs BAM."""

    # Build bwa mem command
    bwa_cmd = ["bwa", "mem", reference, fq1]
    if fq2:
        bwa_cmd.append(fq2)

    # Output to BAM using samtools
    cmd = f"{' '.join(bwa_cmd)} | samtools view -bS - > {output_bam}"

    print(f"Running command:\n{cmd}\n")

    # Run the command
    subprocess.run(cmd, shell=True, check=True)

    print(f"Alignment successful. BAM file saved as: {output_bam}")

def main():
    parser = argparse.ArgumentParser(description="Simple bwa mem wrapper.")
    parser.add_argument("--fq1", required=True, help="FASTQ file R1 (or single-end).")
    parser.add_argument("--fq2", help="FASTQ file R2 (optional for paired-end).")
    parser.add_argument("--ref", required=True, help="Reference genome.")
    parser.add_argument("--output", default="output.bam", help="Output BAM filename.")

    args = parser.parse_args()

    run_bwa_mem(args.fq1, args.fq2, args.ref, args.output)

if __name__ == "__main__":
    main()
