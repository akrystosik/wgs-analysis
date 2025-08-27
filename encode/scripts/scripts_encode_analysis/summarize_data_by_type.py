#!/usr/bin/env python3
"""
summarize_data_by_type.py

Gathers file-size statistics for:
  /mnt/czi-sci-ai/intrinsic-variation-gene-ex/encode_analysis/data
Grouped by "fastq", "bam", "vcf", and "others".
For each type, prints:
  - Sum of sizes (GB)
  - Max file size (GB)
  - Average file size (GB)

Usage:
  chmod +x summarize_data_by_type.py
  ./summarize_data_by_type.py
"""

import os
import sys

DATA_DIR = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/encode_analysis/data"

# 1 GB in bytes (using 1024^3)
GB_DIVISOR = 1_073_741_824

def get_file_type(filename):
    """
    Returns a 'general type' string based on file extension or known pattern.
    Examples:
      - .fastq or .fastq.gz => "fastq"
      - .bam => "bam"
      - .vcf => "vcf"
      otherwise => "others"
    """
    fname_lower = filename.lower()
    if fname_lower.endswith(".fastq") or fname_lower.endswith(".fastq.gz"):
        return "fastq"
    elif fname_lower.endswith(".bam"):
        return "bam"
    elif fname_lower.endswith(".vcf") or fname_lower.endswith(".vcf.gz"):
        return "vcf"
    else:
        return "others"

def main():
    # Ensure the directory exists
    if not os.path.isdir(DATA_DIR):
        print(f"Error: Directory '{DATA_DIR}' does not exist.")
        sys.exit(1)

    # Data structure to hold stats per file type
    # We'll store: total size, max size, and count
    stats_by_type = {
        "fastq": {"total_size": 0, "max_size": 0, "count": 0},
        "bam":   {"total_size": 0, "max_size": 0, "count": 0},
        "vcf":   {"total_size": 0, "max_size": 0, "count": 0},
        "others":{"total_size": 0, "max_size": 0, "count": 0},
    }

    # Walk the directory tree
    for root, dirs, files in os.walk(DATA_DIR):
        for filename in files:
            filepath = os.path.join(root, filename)
            try:
                size_bytes = os.path.getsize(filepath)
            except OSError:
                # If there's a permission or other error, skip
                continue

            # Determine the general file type
            file_type = get_file_type(filename)

            # Update statistics
            stats = stats_by_type[file_type]
            stats["total_size"] += size_bytes
            stats["count"] += 1
            if size_bytes > stats["max_size"]:
                stats["max_size"] = size_bytes

    # Print summary
    print(f"File sizes in '{DATA_DIR}' by type (GB):\n")
    for ftype in ["fastq", "bam", "vcf", "others"]:
        stats = stats_by_type[ftype]
        count = stats["count"]
        if count == 0:
            # No files of this type
            print(f"{ftype}:\n  - No files found.\n")
            continue

        total_size_gb = stats["total_size"] / GB_DIVISOR
        max_size_gb   = stats["max_size"]   / GB_DIVISOR
        avg_size_gb   = (stats["total_size"] / count) / GB_DIVISOR

        print(f"{ftype}:")
        print(f"  Total number of files: {count}")
        print(f"  Sum of sizes:         {total_size_gb:.3f} GB")
        print(f"  Max file size:        {max_size_gb:.3f} GB")
        print(f"  Avg file size:        {avg_size_gb:.3f} GB\n")

if __name__ == "__main__":
    main()
