#!/bin/bash
#encode_analysis/scripts/run_deepvariant_analysis.sh
set -e

# Define our key file paths
BASE_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/encode_analysis"
REFERENCE="${BASE_DIR}/data/reference/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"

# Get current date and time for our log files
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# Simple function to print timestamped messages to both console and log file
log_message() {
    local message="[$(date '+%Y-%m-%d %H:%M:%S')] $1"
    echo "$message" | tee -a "${LOG_FILE}"
}

# Check if we received a sample name
if [ $# -ne 1 ]; then
    echo "Usage: $0 <sample_id>"
    echo "Example: $0 A549_WGS_pair1"
    exit 1
fi

SAMPLE_ID="$1"  # This is the full ID (e.g., A549_WGS_pair1)
SAMPLE_BASE=$(echo "$SAMPLE_ID" | sed 's/_WGS_pair1//')  # Extract base name (e.g., A549)

# Use proper directory and file naming
BAM_FILE="${BASE_DIR}/data/bam/${SAMPLE_BASE}/${SAMPLE_ID}.marked_duplicates.bam"
OUTPUT_DIR="${BASE_DIR}/data/variants/deepvariant"

# Create dated log directory for this run
LOG_DIR="${BASE_DIR}/logs/${TIMESTAMP}_${SAMPLE_ID}"
LOG_FILE="${LOG_DIR}/deepvariant_run.log"
mkdir -p "${LOG_DIR}"

# Record the start of our analysis with key information
{
    echo "=== Running DeepVariant with zombie process protection ==="
    echo "Sample: ${SAMPLE_ID}"
    echo "BAM: ${BAM_FILE}"
    echo "Started: $(date)"
} | tee -a "${LOG_FILE}"

# Create output directory for variants
mkdir -p "${OUTPUT_DIR}"

# Verify input files exist
if [ ! -f "${BAM_FILE}" ]; then
    log_message "ERROR: BAM file not found: ${BAM_FILE}"
    exit 1
fi

if [ ! -f "${REFERENCE}" ]; then
    log_message "ERROR: Reference file not found: ${REFERENCE}"
    exit 1
fi

# Verify BAM index exists, create if needed
if [ ! -f "${BAM_FILE}.bai" ]; then
    log_message "Creating BAM index"
    samtools index "${BAM_FILE}"
fi

# Set environment variables to force CPU usage and better manage memory
export TF_FORCE_GPU_ALLOW_GROWTH=true  # Limit GPU memory growth
export TF_CPP_MIN_LOG_LEVEL=2          # Reduce TensorFlow logging
export OMP_NUM_THREADS=16              # Limit OpenMP threads
export MALLOC_ARENA_MAX=2              # Limit memory fragmentation

# Reduce the number of shards for more stability
NUM_SHARDS=64

# Run DeepVariant with proper output naming - removing the problematic flag
/opt/deepvariant/bin/run_deepvariant \
    --model_type=WGS \
    --ref="${REFERENCE}" \
    --reads="${BAM_FILE}" \
    --output_vcf="${OUTPUT_DIR}/${SAMPLE_ID}.deepvariant.vcf.gz" \
    --output_gvcf="${OUTPUT_DIR}/${SAMPLE_ID}.deepvariant.g.vcf.gz" \
    --intermediate_results_dir="${OUTPUT_DIR}/intermediate_${TIMESTAMP}" \
    --num_shards=${NUM_SHARDS} \
    --verbosity=2 2>&1 | tee -a "${LOG_FILE}"

DV_EXIT=$?
echo "DeepVariant completed with status: $DV_EXIT" | tee -a "${LOG_FILE}"
echo "Finished: $(date)" | tee -a "${LOG_FILE}"

# Check if output was created properly
if [ ! -f "${OUTPUT_DIR}/${SAMPLE_ID}.deepvariant.vcf.gz" ]; then
    echo "ERROR: Output file not created" | tee -a "${LOG_FILE}"
    
    # Clean up zombie processes
    pkill -9 -f "make_examples" || true
    pkill -9 -f "call_variants" || true
    pkill -9 -f "postprocess_variants" || true
    
    exit 1
fi

# Create a summary of the run
{
    echo ""
    echo "=== Run Summary ==="
    echo "Exit status: $DV_EXIT"
    echo "Completed: $(date)"
    echo "Output VCF: ${OUTPUT_DIR}/${SAMPLE_ID}.deepvariant.vcf.gz"
    echo "Output gVCF: ${OUTPUT_DIR}/${SAMPLE_ID}.deepvariant.g.vcf.gz"
    echo "Intermediate files: ${OUTPUT_DIR}/intermediate_${TIMESTAMP}"
} | tee -a "${LOG_FILE}"

# Success!
exit 0