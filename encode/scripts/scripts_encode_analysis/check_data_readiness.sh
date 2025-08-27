#!/bin/bash

# Configuration
CELL_LINES=("K562" "GM23248" "HepG2" "ENCSR121TMQ")
BASE_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/encode_analysis"

echo "=== Data Readiness Check ==="
echo "Checking for required files..."

for cell_line in "${CELL_LINES[@]}"; do
    echo -e "\nChecking ${cell_line}:"
    
    # Check BAM
    if ls ${BASE_DIR}/data/bam/${cell_line}/*.marked_duplicates.bam 1> /dev/null 2>&1; then
        echo "✓ BAM file found"
    else
        echo "✗ BAM file missing"
    fi
    
    # Check DeepVariant
    if ls ${BASE_DIR}/data/variants/deepvariant/${cell_line}*.vcf.gz 1> /dev/null 2>&1; then
        echo "✓ DeepVariant calls found"
    else
        echo "✗ DeepVariant calls missing"
    fi
    
    # Check DeepSomatic
    if ls ${BASE_DIR}/data/variants/deepsomatic/${cell_line}*.vcf.gz 1> /dev/null 2>&1; then
        echo "✓ DeepSomatic calls found"
    else
        echo "✗ DeepSomatic calls missing"
    fi
done

# Check reference
if [ -f "${BASE_DIR}/data/reference/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta" ]; then
    echo -e "\n✓ Reference genome found"
else
    echo -e "\n✗ Reference genome missing"
fi
