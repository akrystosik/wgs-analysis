#!/bin/bash

# Exit on error and print commands for debugging
set -e
set -x

# Base directories - adjust these paths as needed
BASE_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/encode_analysis"
LOG_DIR="${BASE_DIR}/logs"
DATA_DIR="${BASE_DIR}/data"
REFERENCE="${DATA_DIR}/reference/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"

# Generate timestamp for unique run identification
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# Define sample-specific parameters
SAMPLE="K562"
BAM_FILE="${DATA_DIR}/bam/${SAMPLE}/${SAMPLE}_WGS_pair1.marked_duplicates.bam"

# Create directory structure for this run
RUN_DIR="${DATA_DIR}/variants/deepsomatic/runs/${TIMESTAMP}"
INTERMEDIATE_DIR="${RUN_DIR}/intermediate_results"
LOG_PREFIX="${LOG_DIR}/${SAMPLE}_deepsomatic_${TIMESTAMP}"

# Create necessary directories
mkdir -p ${LOG_DIR} ${RUN_DIR} ${INTERMEDIATE_DIR}

# Archive previous runs if they exist
LATEST_LINK="${DATA_DIR}/variants/deepsomatic/latest"
if [ -L "${LATEST_LINK}" ]; then
    mv "${LATEST_LINK}" "${DATA_DIR}/variants/deepsomatic/previous_${TIMESTAMP}"
fi

# Function to run DeepSomatic with enhanced parameters
run_deepsomatic_enhanced() {
    local output_prefix="${RUN_DIR}/${SAMPLE}_WGS_pair1.deepsomatic"
    
    # Run DeepSomatic with optimized parameters for cancer cell lines
    run_deepsomatic \
        --model_type=WGS_TUMOR_ONLY \
        --ref="${REFERENCE}" \
        --reads_tumor="${BAM_FILE}" \
        --output_vcf="${DATA_DIR}/variants/deepsomatic/${SAMPLE}_WGS_pair1.deepsomatic.with_germline.vcf.gz" \
        --output_gvcf="${DATA_DIR}/variants/deepsomatic/${SAMPLE}_WGS_pair1.deepsomatic.with_germline.g.vcf.gz" \
        --sample_name_tumor="${SAMPLE}_WGS_pair1" \
        --num_shards=64 \
        --logging_dir="${LOG_DIR}" \
        --intermediate_results_dir="${INTERMEDIATE_DIR}" \
        --use_default_pon_filtering=true \
        --dry_run=false \
        > "${LOG_PREFIX}.stdout" \
        2> "${LOG_PREFIX}.stderr"
}

# Function to separate and index variants by type
process_variants() {
    local input_vcf="${DATA_DIR}/variants/deepsomatic/${SAMPLE}_WGS_pair1.deepsomatic.with_germline.vcf.gz"
    
    # Extract germline variants
    bcftools view -f GERMLINE "${input_vcf}" -Oz -o "${RUN_DIR}/${SAMPLE}_germline.vcf.gz"
    
    # Extract somatic variants (PASS filter)
    bcftools view -f PASS "${input_vcf}" -Oz -o "${RUN_DIR}/${SAMPLE}_somatic.vcf.gz"
    
    # Index the VCF files
    bcftools index "${RUN_DIR}/${SAMPLE}_germline.vcf.gz"
    bcftools index "${RUN_DIR}/${SAMPLE}_somatic.vcf.gz"
}

# Main execution
echo "Starting DeepSomatic analysis for ${SAMPLE} at ${TIMESTAMP}"

# Run DeepSomatic
run_deepsomatic_enhanced

# Process and separate variants
process_variants

# Create symbolic link to latest run
ln -sf "${RUN_DIR}" "${LATEST_LINK}"

# Generate run summary
echo "=== Run Summary ===" > "${RUN_DIR}/run_summary.txt"
echo "Timestamp: ${TIMESTAMP}" >> "${RUN_DIR}/run_summary.txt"
echo "Sample: ${SAMPLE}" >> "${RUN_DIR}/run_summary.txt"
echo "" >> "${RUN_DIR}/run_summary.txt"
echo "=== Variant Statistics ===" >> "${RUN_DIR}/run_summary.txt"
bcftools stats "${RUN_DIR}/${SAMPLE}_germline.vcf.gz" >> "${RUN_DIR}/run_summary.txt"
bcftools stats "${RUN_DIR}/${SAMPLE}_somatic.vcf.gz" >> "${RUN_DIR}/run_summary.txt"

echo "Analysis complete. Results are in: ${RUN_DIR}"
echo "Check logs at: ${LOG_PREFIX}.{stdout,stderr}"