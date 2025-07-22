#!/usr/bin/env python3

import datetime

def parse_results(results_file):
    """Parse test results from a simple text file."""
    with open(results_file, "r") as f:
        lines = f.readlines()

    results = [line.strip() for line in lines if line.strip()]
    return results

def generate_markdown_report(results, output_file):
    """Generate a Markdown report from test results."""
    date_str = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    with open(output_file, "w") as f:
        f.write(f"# Alignment System Test Report\n\n")
        f.write(f"**Date of Test:** {date_str}\n\n")
        f.write(f"## Test Results Summary\n\n")

        f.write("| Test Case | Result |\n")
        f.write("|-----------|--------|\n")

        for idx, result in enumerate(results, start=1):
            f.write(f"| Test {idx} | {result} |\n")

        f.write("\n\n---\n\n")
        f.write("**Generated automatically by reporter.py.**\n")

    print(f"Markdown report generated: {output_file}")

def main():
    results_file = "test_results.txt"  # Simple text file with PASS/FAIL lines
    output_file = "report/example_report.md"

    results = parse_results(results_file)
    generate_markdown_report(results, output_file)

if __name__ == "__main__":
    main()
