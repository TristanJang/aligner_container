#!/usr/bin/env python3

import argparse
import subprocess
import time
import json

#define alignment command
def run_bwa_mem(fq1, fq2, reference, output_bam):
    """Basic bwa mem wrapper that outputs BAM.""" 
    # Build bwa mem command
    bwa_cmd = ["bwa", "mem", reference, fq1]
    if fq2:
        bwa_cmd.append(fq2)
# Building the command directly as a shell pipeline including memory and CPU usage tracking for REQ 05
    cmd = f"/usr/bin/time -v {' '.join(bwa_cmd)} | samtools view -bS - > {output_bam}"
    print(f"Running command:\n{cmd}\n")

    # Run the command
    #result = subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

    if result.returncode != 0:
        print("\n--- PIPELINE ERROR ---")
        print(result.stderr)
        print("-----------------------")
        raise subprocess.CalledProcessError(result.returncode, cmd)

    print(f"Alignment successful. BAM file saved as: {output_bam}")
    # Parse /usr/bin/time resource usage from stderr
    peak_memory = None
    cpu_percent = None
    for line in result.stderr.split('\n'):
        if "Maximum resident set size" in line:
            peak_memory = int(line.strip().split()[-1]) / 1024  # Convert KB to MB
        if "Percent of CPU this job got" in line:
            cpu_percent = float(line.strip().split()[-1].replace('%', ''))

    return peak_memory, cpu_percent

def alignment_stats(bam_file):
    """Collect basic alignment statistics using samtools flagstat."""
    cmd = ["samtools", "flagstat", bam_file]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    lines = result.stdout.strip().split('\n')
        # NOTE: Using hardcoded index based on current samtools flagstat output format
    total_reads = int(lines[0].split()[0])
    mapped_reads = int(lines[6].split()[0])
    percent_mapped = float(lines[6].split('(')[1].split('%')[0])
    percent_unmapped = 100.0 - percent_mapped

    duplicate_reads = int(lines[3].split()[0])
    percent_duplicates = (duplicate_reads / total_reads) * 100
    return {
        "total_reads": total_reads,
        "mapped_reads": mapped_reads,
        "percent_mapped": percent_mapped,
        "percent_unmapped": percent_unmapped,
        "percent_duplicates": round(percent_duplicates, 2)
    }

def quality_stats(bam_file):
    """Collect average base quality and average MAPQ using samtools stats."""
    cmd = ["samtools", "stats", bam_file]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    lines = result.stdout.strip().split('\n')

    avg_base_quality = None
    avg_mapq = None

    for line in lines:
        if "SN" in line and "average quality" in line:
            avg_base_quality = float(line.strip().split()[-1])
        if "SN" in line and "average mapping quality" in line:
            avg_mapq = float(line.strip().split()[-1])

    return {
        "average_base_quality": avg_base_quality,
        "average_mapping_quality": avg_mapq
    }

def main():
    parser = argparse.ArgumentParser(description="Simple bwa mem wrapper.")
    parser.add_argument("--fq1", required=True, help="FASTQ file R1 (or single-end).")
    parser.add_argument("--fq2", help="FASTQ file R2 (optional for paired-end).")
    parser.add_argument("--ref", required=True, help="Reference genome.")
    parser.add_argument("--output", default="output.bam", help="Output BAM filename.")

    args = parser.parse_args()
    start_time = time.time()
    peak_memory, cpu_percent = run_bwa_mem(args.fq1, args.fq2, args.ref, args.output)
    run_time = time.time() - start_time #getting run time to check REQ 06

#collect stats for json output file
    flagstat_metrics = alignment_stats(args.output)
    quality_metrics = quality_stats(args.output)
    alignment_results = flagstat_metrics  # Start with flagstat metrics
    alignment_results.update(quality_metrics)  # Add quality metrics

# Prepare final results
    results = {
        "alignment_stats": alignment_results,
        "computational_stats": {
            "total_runtime_seconds": round(run_time, 2),
            "peak_memory_MB": peak_memory,
            "cpu_usage_percent": cpu_percent
        }
    }
# Write to JSON
    json_filename = args.output.replace(".bam", "_results.json")
    with open(json_filename, "w") as f:
        json.dump(results, f, indent=4)

    print(f"\nResults written to {json_filename}")

if __name__ == "__main__":
    main()
