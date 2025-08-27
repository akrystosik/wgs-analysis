#!/bin/bash
set -e

SAMPLE_ID="$1"
BASE_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/encode_analysis"
BAM_FILE="${BASE_DIR}/data/bam/${SAMPLE_ID%_WGS_pair1}/${SAMPLE_ID}.marked_duplicates.bam"
REFERENCE="${BASE_DIR}/data/reference/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
OUTPUT_DIR="${BASE_DIR}/data/variants/deepvariant"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_DIR="${BASE_DIR}/logs/${TIMESTAMP}_${SAMPLE_ID}"

mkdir -p "${LOG_DIR}" "${OUTPUT_DIR}/intermediate_${TIMESTAMP}"

echo "=== Running DeepVariant with zombie process protection ===" | tee -a "${LOG_DIR}/deepvariant_run.log"
echo "Sample: ${SAMPLE_ID}" | tee -a "${LOG_DIR}/deepvariant_run.log"
echo "BAM: ${BAM_FILE}" | tee -a "${LOG_DIR}/deepvariant_run.log"
echo "Started: $(date)" | tee -a "${LOG_DIR}/deepvariant_run.log"

# Run with proper resource monitoring
/opt/deepvariant/bin/run_deepvariant \
    --model_type=WGS \
    --ref="${REFERENCE}" \
    --reads="${BAM_FILE}" \
    --output_vcf="${OUTPUT_DIR}/${SAMPLE_ID}.deepvariant.vcf.gz" \
    --output_gvcf="${OUTPUT_DIR}/${SAMPLE_ID}.deepvariant.g.vcf.gz" \
    --intermediate_results_dir="${OUTPUT_DIR}/intermediate_${TIMESTAMP}" \
    --num_shards=64 \
    --runtime_report=true \
    --verbosity=2 2>&1 | tee -a "${LOG_DIR}/deepvariant_run.log"

RESULT=$?
echo "DeepVariant completed with status: ${RESULT}" | tee -a "${LOG_DIR}/deepvariant_run.log"
echo "Finished: $(date)" | tee -a "${LOG_DIR}/deepvariant_run.log"

# Check for output file
if [ -f "${OUTPUT_DIR}/${SAMPLE_ID}.deepvariant.vcf.gz" ]; then
    echo "Output file exists: ${OUTPUT_DIR}/${SAMPLE_ID}.deepvariant.vcf.gz" | tee -a "${LOG_DIR}/deepvariant_run.log"
    exit 0
else
    echo "ERROR: Output file not created" | tee -a "${LOG_DIR}/deepvariant_run.log"
    exit 1
fi
