#!/bin/bash
# run_complete_analysis.sh
#
# Purpose: Orchestrate the complete analysis pipeline for cancer cell lines,
# coordinating DeepSomatic and DeepVariant analysis with proper logging
# and error handling.

set -e

# Base configuration
BASE_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/encode_analysis"
SCRIPTS_DIR="${BASE_DIR}/scripts"
SAMPLE="K562"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_ROOT="${BASE_DIR}/logs"
RUN_LOG_DIR="${LOG_ROOT}/${TIMESTAMP}"

# Create comprehensive logging structure
mkdir -p "${RUN_LOG_DIR}"/{deepsomatic,variant_processing,summary}

# Initialize run manifest
cat > "${RUN_LOG_DIR}/manifest.txt" << EOF
Analysis Run Manifest
====================
Timestamp: ${TIMESTAMP}
Sample: ${SAMPLE}
Started: $(date)

Pipeline Configuration:
- DeepSomatic Mode: WGS_TUMOR_ONLY
- PoN Filtering: Enabled
- Version: v3 (Following recommended methodology)

Input Validation:
EOF

# Verify input files exist
verify_inputs() {
    local bam_file="${BASE_DIR}/data/bam/${SAMPLE}/${SAMPLE}_WGS_pair1.marked_duplicates.bam"
    local reference="${BASE_DIR}/data/reference/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
    
    [[ -f "${bam_file}" ]] && echo "BAM file: Found" >> "${RUN_LOG_DIR}/manifest.txt" || \
        { echo "BAM file: Missing" >> "${RUN_LOG_DIR}/manifest.txt"; return 1; }
    [[ -f "${reference}" ]] && echo "Reference: Found" >> "${RUN_LOG_DIR}/manifest.txt" || \
        { echo "Reference: Missing" >> "${RUN_LOG_DIR}/manifest.txt"; return 1; }
}

# Execute analysis with proper error handling
execute_pipeline() {
    local stage=$1
    local cmd=$2
    local log_prefix="${RUN_LOG_DIR}/${stage}"
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting ${stage}"
    if ${cmd} > "${log_prefix}.stdout" 2> "${log_prefix}.stderr"; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Completed ${stage}"
        return 0
    else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Failed ${stage}"
        return 1
    fi
}

# Main execution
echo "=== Starting Complete Analysis Pipeline for ${SAMPLE} ==="
echo "Timestamp: ${TIMESTAMP}"

# Verify inputs
if ! verify_inputs; then
    echo "Error: Required input files missing. Check manifest for details."
    exit 1
fi

# Step 1: Run DeepSomatic analysis
if ! execute_pipeline "deepsomatic" "${SCRIPTS_DIR}/run_deepsomatic_analysis.sh"; then
    echo "Error: DeepSomatic analysis failed. Check logs for details."
    exit 1
fi

# Allow time for file system synchronization
sleep 10

# Step 2: Run variant integration pipeline
if ! execute_pipeline "variant_processing" "${SCRIPTS_DIR}/cell_line_variants_pipeline_with_germlinev3.sh"; then
    echo "Error: Variant processing failed. Check logs for details."
    exit 1
fi

# Generate combined summary
SUMMARY_FILE="${RUN_LOG_DIR}/summary/complete_analysis_summary.txt"
echo "=== Combined Analysis Summary for ${SAMPLE} ===" > "${SUMMARY_FILE}"
echo "Timestamp: ${TIMESTAMP}" >> "${SUMMARY_FILE}"
echo "" >> "${SUMMARY_FILE}"

# Add summaries from both processes
cat "${BASE_DIR}/data/variants/deepsomatic/runs/latest/run_summary.txt" >> "${SUMMARY_FILE}"
echo "" >> "${SUMMARY_FILE}"
echo "=== Final Variant Statistics ===" >> "${SUMMARY_FILE}"
cat "${BASE_DIR}/data/variants/cell_line_analysis/${SAMPLE}/final/analysis_summary.txt" >> "${SUMMARY_FILE}"

# Create success marker
touch "${RUN_LOG_DIR}/pipeline_completed_successfully"

echo "Complete analysis finished successfully."
echo "Review complete summary at: ${SUMMARY_FILE}"
echo "All logs available in: ${RUN_LOG_DIR}"