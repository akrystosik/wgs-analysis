#!/bin/bash
#encode_analysis/scripts/run_deepsomatic_analysis.sh
# Exit on error and print commands for debugging
set -e
set -x

# Base directories
BASE_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/encode_analysis"
REFERENCE="${BASE_DIR}/data/reference/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"

# Generate timestamp for unique run identification
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# Check if we received a sample name
if [ $# -ne 1 ]; then
    echo "Usage: $0 <sample_id>"
    echo "Example: $0 A549_WGS_pair1"
    exit 1
fi

# Handle sample naming consistently with deepvariant script
SAMPLE_ID="$1"  # This is the full ID (e.g., A549_WGS_pair1)
SAMPLE_BASE=$(echo "$SAMPLE_ID" | sed 's/_WGS_pair1//')  # Extract base name (e.g., A549)

# Define paths using consistent naming
BAM_FILE="${BASE_DIR}/data/bam/${SAMPLE_BASE}/${SAMPLE_ID}.marked_duplicates.bam"
OUTPUT_DIR="${BASE_DIR}/data/variants/deepsomatic"
LOG_DIR="${BASE_DIR}/logs/${TIMESTAMP}_${SAMPLE_ID}"
LOG_FILE="${LOG_DIR}/deepsomatic_run.log"

# Create necessary directories
mkdir -p "${LOG_DIR}" "${OUTPUT_DIR}"

# Function to log messages
log_message() {
    local message="[$(date '+%Y-%m-%d %H:%M:%S')] $1"
    echo "$message" | tee -a "${LOG_FILE}"
}

# Record the start of our analysis
{
    echo "=== DeepSomatic Analysis Run ==="
    echo "Sample ID: ${SAMPLE_ID}"
    echo "Base Sample Name: ${SAMPLE_BASE}"
    echo "Started: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "Command directory: $(pwd)"
    echo "Reference: ${REFERENCE}"
    echo "BAM file: ${BAM_FILE}"
    echo "==============================="
    echo ""
} > "${LOG_FILE}"

# Verify input files exist
if [ ! -f "${BAM_FILE}" ]; then
    log_message "ERROR: BAM file not found: ${BAM_FILE}"
    exit 1
fi

if [ ! -f "${REFERENCE}" ]; then
    log_message "ERROR: Reference file not found: ${REFERENCE}"
    exit 1
fi

# Run DeepSomatic with consistent output naming
log_message "Starting DeepSomatic analysis for ${SAMPLE_ID}"

/opt/deepvariant/bin/deepsomatic/run_deepsomatic \
    --model_type=WGS_TUMOR_ONLY \
    --ref="${REFERENCE}" \
    --reads_tumor="${BAM_FILE}" \
    --output_vcf="${OUTPUT_DIR}/${SAMPLE_ID}.deepsomatic.vcf.gz" \
    --output_gvcf="${OUTPUT_DIR}/${SAMPLE_ID}.deepsomatic.g.vcf.gz" \
    --sample_name_tumor="${SAMPLE_ID}" \
    --num_shards=128 \
    --intermediate_results_dir="${OUTPUT_DIR}/intermediate_${TIMESTAMP}" \
    --verbosity=2 2>&1 | tee -a "${LOG_FILE}"

# Record completion
log_message "DeepSomatic analysis completed"
log_message "Log file location: ${LOG_FILE}"

# Create a summary of the run
{
    echo ""
    echo "=== Run Summary ==="
    echo "Completed: $(date '+%Y-%m-%d %H:%M:%S')"
    echo "Output VCF: ${OUTPUT_DIR}/${SAMPLE_ID}.deepsomatic.vcf.gz"
    echo "Output gVCF: ${OUTPUT_DIR}/${SAMPLE_ID}.deepsomatic.g.vcf.gz"
    echo "Intermediate files: ${OUTPUT_DIR}/intermediate_${TIMESTAMP}"
} >> "${LOG_FILE}"