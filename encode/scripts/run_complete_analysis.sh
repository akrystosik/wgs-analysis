#!/bin/bash
set -e

# Base configuration
BASE_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/encode_analysis"
SCRIPTS_DIR="${BASE_DIR}/scripts"
SAMPLE="K562"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# Ensure we're in the correct directory
cd ${BASE_DIR}

echo "=== Starting Complete Analysis Pipeline for ${SAMPLE} ==="
echo "Timestamp: ${TIMESTAMP}"

# Step 1: Run DeepSomatic analysis
echo "Running DeepSomatic analysis..."
${SCRIPTS_DIR}/run_deepsomatic_analysis.sh

# Wait for files to be properly written and indexed
sleep 10

# Step 2: Run variant integration pipeline
echo "Running variant integration pipeline..."
${SCRIPTS_DIR}/cell_line_variants_pipeline_with_germlinev2.sh

# Generate combined summary
SUMMARY_FILE="${BASE_DIR}/analysis_summary_${TIMESTAMP}.txt"
echo "=== Combined Analysis Summary for ${SAMPLE} ===" > ${SUMMARY_FILE}
echo "Timestamp: ${TIMESTAMP}" >> ${SUMMARY_FILE}
echo "" >> ${SUMMARY_FILE}

# Add summaries from both processes
cat "${BASE_DIR}/data/variants/deepsomatic/runs/latest/run_summary.txt" >> ${SUMMARY_FILE}
echo "" >> ${SUMMARY_FILE}
echo "=== Final Variant Statistics ===" >> ${SUMMARY_FILE}
cat "${BASE_DIR}/data/variants/cell_line_analysis/${SAMPLE}/final/analysis_summary.txt" >> ${SUMMARY_FILE}

echo "Complete analysis finished. Check ${SUMMARY_FILE} for results"

# Create a final report with key metrics
echo "Creating final report..."
echo "=== ${SAMPLE} Analysis Report - ${TIMESTAMP} ===" > ${BASE_DIR}/final_report_${TIMESTAMP}.txt
echo "Pipeline completed successfully" >> ${BASE_DIR}/final_report_${TIMESTAMP}.txt
echo "" >> ${BASE_DIR}/final_report_${TIMESTAMP}.txt
echo "Key locations:" >> ${BASE_DIR}/final_report_${TIMESTAMP}.txt
echo "- DeepSomatic results: ${BASE_DIR}/data/variants/deepsomatic/latest" >> ${BASE_DIR}/final_report_${TIMESTAMP}.txt
echo "- Final variants: ${BASE_DIR}/data/variants/cell_line_analysis/${SAMPLE}/final" >> ${BASE_DIR}/final_report_${TIMESTAMP}.txt
echo "- Analysis summary: ${SUMMARY_FILE}" >> ${BASE_DIR}/final_report_${TIMESTAMP}.txt

echo "Pipeline complete. Check final_report_${TIMESTAMP}.txt for details"