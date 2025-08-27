#!/bin/bash

# Exit on error
set -e

# Set up directories
CHAIN_FILE="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/encode_analysis/data/reference/hg19ToHg38.over.chain.gz"
OUTPUT_DIR="data/variants/encode_lifted"
mkdir -p $OUTPUT_DIR

# Check for chain file
if [ ! -f "$CHAIN_FILE" ]; then
    echo "Chain file not found at $CHAIN_FILE"
    exit 1
fi

# Array of ENCODE VCFs
ENCODE_VCFS=(
    "ENCFF752OAX"
    "ENCFF785JVR"
    "ENCFF574MDJ"
    "ENCFF863MPP"
)

# Process each VCF
for vcf in "${ENCODE_VCFS[@]}"; do
    echo "Processing $vcf..."
    
    # Lift over to hg38
    CrossMap vcf \
        "$CHAIN_FILE" \
        data/variants/encode/${vcf}.vcf.gz \
        data/reference/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
        $OUTPUT_DIR/${vcf}.hg38.unsorted.vcf

    if [ ! -f "$OUTPUT_DIR/${vcf}.hg38.unsorted.vcf" ]; then
        echo "ERROR: Liftover failed for $vcf"
        continue
    fi

    echo "Sorting ${vcf}..."
    # Sort VCF and handle chromosomes
    (grep "^#" $OUTPUT_DIR/${vcf}.hg38.unsorted.vcf; 
     grep -v "^#" $OUTPUT_DIR/${vcf}.hg38.unsorted.vcf | sort -k1,1 -k2,2n) \
    > $OUTPUT_DIR/${vcf}.hg38.sorted.vcf

    echo "Compressing ${vcf}..."
    bgzip -f $OUTPUT_DIR/${vcf}.hg38.sorted.vcf
    
    echo "Indexing ${vcf}..."
    tabix -f $OUTPUT_DIR/${vcf}.hg38.sorted.vcf.gz
    
    # Generate stats
    echo "Generating stats for ${vcf}..."
    bcftools stats $OUTPUT_DIR/${vcf}.hg38.sorted.vcf.gz > $OUTPUT_DIR/${vcf}.stats
    
    # Clean up
    rm -f $OUTPUT_DIR/${vcf}.hg38.unsorted.vcf
done

# Create summary report
{
    echo "ENCODE VCF Liftover Summary"
    echo "=========================="
    echo
    echo "Original files (hg19):"
    for vcf in "${ENCODE_VCFS[@]}"; do
        echo -n "${vcf}: "
        bcftools view -H data/variants/encode/${vcf}.vcf.gz | wc -l
    done
    
    echo -e "\nLifted files (hg38):"
    for vcf in "${ENCODE_VCFS[@]}"; do
        if [ -f "$OUTPUT_DIR/${vcf}.hg38.sorted.vcf.gz" ]; then
            echo -n "${vcf}: "
            bcftools view -H $OUTPUT_DIR/${vcf}.hg38.sorted.vcf.gz | wc -l
        else
            echo "${vcf}: Failed to process"
        fi
    done
    
    echo -e "\nMapping Statistics:"
    for vcf in "${ENCODE_VCFS[@]}"; do
        if [ -f "$OUTPUT_DIR/${vcf}.stats" ]; then
            echo -e "\n${vcf}:"
            grep "^SN" $OUTPUT_DIR/${vcf}.stats | cut -f 3-
        fi
    done
} > $OUTPUT_DIR/liftover_summary.txt

echo "Processing complete! Check $OUTPUT_DIR/liftover_summary.txt for results"