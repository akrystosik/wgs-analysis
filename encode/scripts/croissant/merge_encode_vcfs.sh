#!/bin/bash
#
# Merge ENCODE cell line VCF files using bcftools
#
# This script merges the 9 ENCODE cell line combined variant files into a single 
# multi-sample VCF using standard computational biology tools (bcftools).
#

set -euo pipefail

# Define the base directory and output paths
BASE_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/wgs/encode/data/variants"
CELL_LINE_DIR="${BASE_DIR}/cell_line_analysis_v5"
OUTPUT_DIR="${BASE_DIR}/merged_encode"
OUTPUT_VCF="${OUTPUT_DIR}/encode_combined_variants.vcf.gz"

# Cell lines to merge
CELL_LINES=("A549" "Caki2" "GM23248" "HepG2" "K562" "NCI-H460" "Panc1" "T47D" "sknmc")

echo "Starting ENCODE VCF merge process..."
echo "Output directory: ${OUTPUT_DIR}"

# Create output directory
mkdir -p "${OUTPUT_DIR}"

# Check that all input files exist
echo "Checking input files..."
for cell_line in "${CELL_LINES[@]}"; do
    input_vcf="${CELL_LINE_DIR}/${cell_line}/final/combined_variants.vcf.gz"
    if [[ ! -f "${input_vcf}" ]]; then
        echo "Error: Missing VCF file for ${cell_line}: ${input_vcf}"
        exit 1
    fi
    if [[ ! -f "${input_vcf}.tbi" ]]; then
        echo "Error: Missing index file for ${cell_line}: ${input_vcf}.tbi"
        exit 1
    fi
    echo "  ✓ ${cell_line}: ${input_vcf}"
done

# Create a list of input files for bcftools merge
input_list="${OUTPUT_DIR}/input_files.txt"
echo "Creating input file list..."
for cell_line in "${CELL_LINES[@]}"; do
    echo "${CELL_LINE_DIR}/${cell_line}/final/combined_variants.vcf.gz" >> "${input_list}"
done

# Merge VCF files using bcftools
echo "Merging VCF files using bcftools..."
echo "Command: bcftools merge --file-list ${input_list} -Oz -o ${OUTPUT_VCF}"

bcftools merge \
    --file-list "${input_list}" \
    --output-type z \
    --output "${OUTPUT_VCF}"

# Create index for merged file
echo "Creating index for merged VCF..."
tabix -p vcf "${OUTPUT_VCF}"

# Generate summary statistics
echo "Generating summary statistics..."
echo "=== ENCODE VCF Merge Summary ===" > "${OUTPUT_DIR}/merge_summary.txt"
echo "Date: $(date)" >> "${OUTPUT_DIR}/merge_summary.txt"
echo "" >> "${OUTPUT_DIR}/merge_summary.txt"

echo "Input files:" >> "${OUTPUT_DIR}/merge_summary.txt"
for cell_line in "${CELL_LINES[@]}"; do
    input_vcf="${CELL_LINE_DIR}/${cell_line}/final/combined_variants.vcf.gz"
    variant_count=$(bcftools view -H "${input_vcf}" | wc -l)
    echo "  ${cell_line}: ${variant_count} variants" >> "${OUTPUT_DIR}/merge_summary.txt"
done

echo "" >> "${OUTPUT_DIR}/merge_summary.txt"
echo "Output file: ${OUTPUT_VCF}" >> "${OUTPUT_DIR}/merge_summary.txt"

# Get merged file statistics
merged_variants=$(bcftools view -H "${OUTPUT_VCF}" | wc -l)
merged_samples=$(bcftools query -l "${OUTPUT_VCF}" | wc -l)

echo "Merged variants: ${merged_variants}" >> "${OUTPUT_DIR}/merge_summary.txt"
echo "Merged samples: ${merged_samples}" >> "${OUTPUT_DIR}/merge_summary.txt"

echo "" >> "${OUTPUT_DIR}/merge_summary.txt"
echo "Sample names in merged file:" >> "${OUTPUT_DIR}/merge_summary.txt"
bcftools query -l "${OUTPUT_VCF}" >> "${OUTPUT_DIR}/merge_summary.txt"

# Display summary
echo ""
echo "=== Merge Complete ==="
echo "Merged VCF: ${OUTPUT_VCF}"
echo "Total variants: ${merged_variants}"
echo "Total samples: ${merged_samples}"
echo "Summary saved to: ${OUTPUT_DIR}/merge_summary.txt"

echo ""
echo "ENCODE VCF merge completed successfully!"