#!/bin/bash
# debug_deepsomatic.sh
set -e

# Define constants
SAMPLE_ID="$1"
BASE_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/encode_analysis"
BAM_FILE="${BASE_DIR}/data/bam/${SAMPLE_ID%_WGS_pair1}/${SAMPLE_ID}.marked_duplicates.bam"
REFERENCE="${BASE_DIR}/data/reference/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
OUTPUT_DIR="${BASE_DIR}/data/variants/deepsomatic"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_DIR="${BASE_DIR}/logs/debug_${TIMESTAMP}_${SAMPLE_ID}"
LOG_FILE="${LOG_DIR}/deepsomatic_debug.log"

# Create directories
mkdir -p "${LOG_DIR}" "${OUTPUT_DIR}"

# Log function
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "${LOG_FILE}"
}

# System resource monitoring function
monitor_resources() {
    while true; do
        echo "==== $(date) ====" >> "${LOG_DIR}/system_stats.log"
        echo "Memory usage:" >> "${LOG_DIR}/system_stats.log"
        free -h >> "${LOG_DIR}/system_stats.log"
        echo "CPU usage:" >> "${LOG_DIR}/system_stats.log"
        top -b -n 1 | head -n 20 >> "${LOG_DIR}/system_stats.log"
        echo "Disk usage:" >> "${LOG_DIR}/system_stats.log"
        df -h >> "${LOG_DIR}/system_stats.log"
        echo "Process count:" >> "${LOG_DIR}/system_stats.log"
        ps aux | grep -E "make_examples|call_variants|post" | wc -l >> "${LOG_DIR}/system_stats.log"
        echo "Zombie process count:" >> "${LOG_DIR}/system_stats.log"
        ps aux | grep defunct | wc -l >> "${LOG_DIR}/system_stats.log"
        sleep 30
    done
}

# Start monitor in background
log_message "Starting system resource monitoring"
monitor_resources &
MONITOR_PID=$!

# Ensure we clean up the monitor on exit
trap "kill $MONITOR_PID; log_message 'Monitoring stopped'" EXIT

# Set environment variables for better performance
export TF_FORCE_GPU_ALLOW_GROWTH=true
export TF_CPP_MIN_LOG_LEVEL=1
export OMP_NUM_THREADS=16
export MALLOC_ARENA_MAX=2

# Record start information
log_message "Starting DeepSomatic debug run for ${SAMPLE_ID}"
log_message "Using reference: ${REFERENCE}"
log_message "Using BAM file: ${BAM_FILE}"
log_message "Using reduced parallelism (32 shards)"
log_message "System information:"
uname -a >> "${LOG_FILE}"
cat /proc/meminfo | grep MemTotal >> "${LOG_FILE}"
cat /proc/cpuinfo | grep "model name" | head -1 >> "${LOG_FILE}"

# Run with reduced parallelism and timeout protection
log_message "Running DeepSomatic with output to: ${OUTPUT_DIR}/${SAMPLE_ID}.deepsomatic.vcf.gz"

timeout --kill-after=10m 12h /opt/deepvariant/bin/deepsomatic/run_deepsomatic \
    --model_type=WGS_TUMOR_ONLY \
    --ref="${REFERENCE}" \
    --reads_tumor="${BAM_FILE}" \
    --output_vcf="${OUTPUT_DIR}/${SAMPLE_ID}.deepsomatic.vcf.gz" \
    --output_gvcf="${OUTPUT_DIR}/${SAMPLE_ID}.deepsomatic.g.vcf.gz" \
    --sample_name_tumor="${SAMPLE_ID}" \
    --num_shards=32 \
    --intermediate_results_dir="${OUTPUT_DIR}/intermediate_${TIMESTAMP}" \
    --verbosity=2 2>&1 | tee -a "${LOG_FILE}"

EXIT_CODE=$?

if [ $EXIT_CODE -eq 124 ] || [ $EXIT_CODE -eq 137 ]; then
    log_message "ERROR: DeepSomatic timed out after 12 hours"
elif [ $EXIT_CODE -ne 0 ]; then
    log_message "ERROR: DeepSomatic failed with exit code ${EXIT_CODE}"
else
    log_message "DeepSomatic completed successfully"
fi

# Check for output
if [ -f "${OUTPUT_DIR}/${SAMPLE_ID}.deepsomatic.vcf.gz" ]; then
    log_message "Output VCF file created successfully"
    ls -lh "${OUTPUT_DIR}/${SAMPLE_ID}.deepsomatic.vcf.gz" >> "${LOG_FILE}"
else
    log_message "ERROR: Output VCF file was not created"
fi

# Clean up any zombie processes
pkill -9 -f "make_examples_somatic" || true

log_message "Debug run completed. See logs at ${LOG_DIR}"