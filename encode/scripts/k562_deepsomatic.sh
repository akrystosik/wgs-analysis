#!/bin/bash
set -e  # Exit on error
set -x  # Print commands for debugging

LOG_DIR=/mnt/czi-sci-ai/intrinsic-variation-gene-ex/encode_analysis/logs
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_PREFIX=${LOG_DIR}/K562_deepsomatic_${TIMESTAMP}

# Ensure log directory exists
mkdir -p ${LOG_DIR}

# Clean up intermediate results
rm -rf /mnt/czi-sci-ai/intrinsic-variation-gene-ex/encode_analysis/data/variants/deepsomatic/intermediate_results/*

# Modified DeepSomatic command
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
    --dry_run=false

echo "Run completed. Check logs at:"
echo "STDOUT: ${LOG_PREFIX}.stdout"
echo "STDERR: ${LOG_PREFIX}.stderr"

# After completion, verify GERMLINE filter presence
echo "=== Checking VCF Header ==="
bcftools view -h /mnt/czi-sci-ai/intrinsic-variation-gene-ex/encode_analysis/data/variants/deepsomatic/K562_WGS_pair1.deepsomatic.with_germline.vcf.gz | grep "FILTER"

echo "=== Checking Filter Distribution ==="
bcftools view /mnt/czi-sci-ai/intrinsic-variation-gene-ex/encode_analysis/data/variants/deepsomatic/K562_WGS_pair1.deepsomatic.with_germline.vcf.gz | \
    grep -v "^#" | cut -f7 | sort | uniq -c
