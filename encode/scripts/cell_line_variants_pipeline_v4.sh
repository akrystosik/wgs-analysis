#!/bin/bash
#encode_analysis/scripts/cell_line_variants_pipeline_v4.sh

# Purpose: Process variant calls from DeepSomatic and DeepVariant for cancer cell lines.
# This script combines variants from both callers, prioritizing DeepSomatic's
# somatic calls while preserving unique germline variants from DeepVariant.

set -e
set -o pipefail

# Base configuration
BASE_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/encode_analysis"
OUTPUT_ROOT="${BASE_DIR}/data/variants/cell_line_analysis_v4"

# Define our supported cell lines and their file paths
declare -A DEEPSOMATIC_FILES DEEPVARIANT_FILES
# Each cell line needs both DeepSomatic and DeepVariant results
DEEPSOMATIC_FILES["K562"]="data/variants/deepsomatic/K562_WGS_pair1.deepsomatic.with_germline.vcf.gz"
DEEPSOMATIC_FILES["Caki2"]="data/variants/deepsomatic/Caki2_somatic_full.vcf"
DEEPSOMATIC_FILES["HepG2"]="data/variants/deepsomatic/HepG2_WGS_pair1.vcf.gz"

DEEPVARIANT_FILES["K562"]="data/variants/deepvariant/K562_WGS_pair1.deepvariant.vcf.gz"
DEEPVARIANT_FILES["Caki2"]="data/variants/deepvariant/Caki2_WGS_pair1.vcf.gz"
DEEPVARIANT_FILES["HepG2"]="data/variants/deepvariant/HepG2_WGS_pair1.deepvariant.vcf.gz"

# Function to print timestamped log messages
log_step() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
    echo "----------------------------------------"
}

# Function to display usage information
usage() {
    echo "Usage: $0 -c <cell_line>"
    echo "Available cell lines: ${!DEEPSOMATIC_FILES[@]}"
    echo "Example: $0 -c HepG2"
    exit 1
}

# Process command line arguments
if [[ $# -eq 0 ]]; then
    usage
fi

CELL_LINE=""
while getopts "c:" opt; do
    case $opt in
        c) CELL_LINE="$OPTARG" ;;
        *) usage ;;
    esac
done

# Validate cell line selection
if [[ ! "${DEEPSOMATIC_FILES[$CELL_LINE]+isset}" ]]; then
    echo "Error: Invalid cell line. Choose from: ${!DEEPSOMATIC_FILES[@]}"
    exit 1
fi

# Setup for specific cell line
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
OUTDIR="${OUTPUT_ROOT}/${CELL_LINE}"
LOG_DIR="${BASE_DIR}/logs/${TIMESTAMP}_${CELL_LINE}"

# Input files with full paths
DEEPSOMATIC_VCF="${BASE_DIR}/${DEEPSOMATIC_FILES[$CELL_LINE]}"
DEEPVARIANT_VCF="${BASE_DIR}/${DEEPVARIANT_FILES[$CELL_LINE]}"

# Print file paths for verification
echo "Analysis Configuration:"
echo "======================"
echo "Cell Line: ${CELL_LINE}"
echo "DeepSomatic file: ${DEEPSOMATIC_VCF}"
echo "DeepVariant file: ${DEEPVARIANT_VCF}"
echo "Output directory: ${OUTDIR}"
echo "Log directory: ${LOG_DIR}"
echo "======================"

# Setup logging
mkdir -p "${LOG_DIR}"
exec 1> >(tee "${LOG_DIR}/pipeline_execution.log")
exec 2> >(tee "${LOG_DIR}/pipeline_errors.log")

log_step() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
    echo "-------------------------------------------"
}

# Create output structure
mkdir -p "${OUTDIR}"/{intermediate,final}

# Main pipeline steps
log_step "Processing ${CELL_LINE} (${sample_name})"

# Step 1: Extract all DeepVariant PASS calls
log_step "Extracting all DeepVariant PASS calls"
bcftools view -f PASS "${DEEPVARIANT_VCF}" \
    -Oz -o "${OUTDIR}/intermediate/deepvariant_all.vcf.gz"
bcftools index -t "${OUTDIR}/intermediate/deepvariant_all.vcf.gz"

# Step 2: Extract DeepSomatic variants by type
log_step "Extracting DeepSomatic variants by type"
bcftools view -f PASS "${DEEPSOMATIC_VCF}" \
    -Oz -o "${OUTDIR}/intermediate/deepsomatic_somatic.vcf.gz"
bcftools index -t "${OUTDIR}/intermediate/deepsomatic_somatic.vcf.gz"

bcftools view -i 'FILTER=="GERMLINE"' "${DEEPSOMATIC_VCF}" \
    -Oz -o "${OUTDIR}/intermediate/deepsomatic_germline.vcf.gz"
bcftools index -t "${OUTDIR}/intermediate/deepsomatic_germline.vcf.gz"

# Step 3: Create positions file for overlapping variants
log_step "Identifying overlapping positions"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' \
    "${OUTDIR}/intermediate/deepsomatic_somatic.vcf.gz" > \
    "${OUTDIR}/intermediate/somatic_positions.tsv"

# Step 4: Extract non-overlapping DeepVariant calls
log_step "Creating final variant set"
bcftools view \
    -T ^"${OUTDIR}/intermediate/somatic_positions.tsv" \
    "${OUTDIR}/intermediate/deepvariant_all.vcf.gz" \
    -Oz -o "${OUTDIR}/intermediate/deepvariant_unique.vcf.gz"
bcftools index -t "${OUTDIR}/intermediate/deepvariant_unique.vcf.gz"

# Step 5: Normalize sample names and combine
bcftools reheader \
    -s <(echo "${CELL_LINE}") \
    "${OUTDIR}/intermediate/deepvariant_unique.vcf.gz" \
    -o "${OUTDIR}/intermediate/deepvariant_unique_renamed.vcf.gz"

bcftools reheader \
    -s <(echo "${CELL_LINE}") \
    "${OUTDIR}/intermediate/deepsomatic_somatic.vcf.gz" \
    -o "${OUTDIR}/intermediate/deepsomatic_somatic_renamed.vcf.gz"

bcftools index -t "${OUTDIR}/intermediate/deepvariant_unique_renamed.vcf.gz"
bcftools index -t "${OUTDIR}/intermediate/deepsomatic_somatic_renamed.vcf.gz"

bcftools concat \
    -a \
    "${OUTDIR}/intermediate/deepvariant_unique_renamed.vcf.gz" \
    "${OUTDIR}/intermediate/deepsomatic_somatic_renamed.vcf.gz" \
    -Oz -o "${OUTDIR}/final/combined_variants.vcf.gz"
bcftools index -t "${OUTDIR}/final/combined_variants.vcf.gz"

# Generate summary statistics
log_step "Generating summary statistics"
{
    echo "=== Variant Analysis Summary v4 ==="
    echo "Timestamp: ${TIMESTAMP}"
    echo "Cell Line: ${CELL_LINE}"
    echo "Sample ID: ${sample_name}"
    echo ""
    echo "1. Total DeepVariant PASS calls: $(bcftools view -H ${OUTDIR}/intermediate/deepvariant_all.vcf.gz | wc -l)"
    echo ""
    echo "2. DeepSomatic variant counts:"
    echo "   - Somatic (PASS) variants: $(bcftools view -H ${OUTDIR}/intermediate/deepsomatic_somatic.vcf.gz | wc -l)"
    echo "   - Germline-labeled variants: $(bcftools view -H ${OUTDIR}/intermediate/deepsomatic_germline.vcf.gz | wc -l)"
    echo ""
    echo "3. Final combined set:"
    echo "   - Non-overlapping DeepVariant calls: $(bcftools view -H ${OUTDIR}/intermediate/deepvariant_unique.vcf.gz | wc -l)"
    echo "   - Total variants in final set: $(bcftools view -H ${OUTDIR}/final/combined_variants.vcf.gz | wc -l)"
} > "${OUTDIR}/final/analysis_summary.txt"

log_step "Pipeline completed successfully for ${CELL_LINE}"