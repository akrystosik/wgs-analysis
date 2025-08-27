#!/bin/bash

# DeepMap Processing Pipeline
# Usage: ./deepmap.sh

set -e  # Exit on error

# Configuration
PYTHON_SCRIPT="scripts/process_deepmap.py"
CCLE_FILE="data/CCLE_mutations.csv"
SAMPLE_INFO="data/sample_info.csv"
OUTPUT_DIR="data/variants/k562_validation"
VCF_FILE="data/variants/deepsomatic.vcf"

# Check if Python script exists
if [ ! -f "$PYTHON_SCRIPT" ]; then
    echo "Error: Python script not found at $PYTHON_SCRIPT"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Run Python processing script
echo "Starting DeepMap validation pipeline..."
python3 "$PYTHON_SCRIPT" \
    --ccle "$CCLE_FILE" \
    --sample-info "$SAMPLE_INFO" \
    --vcf "$VCF_FILE" \
    --output "$OUTPUT_DIR"

# Check exit status
if [ $? -eq 0 ]; then
    echo "Validation pipeline completed successfully"
    echo "Results available in $OUTPUT_DIR"
else
    echo "Error: Pipeline failed"
    exit 1
fi