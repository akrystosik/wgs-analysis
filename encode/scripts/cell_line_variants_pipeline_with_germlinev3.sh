#!/bin/bash
# cell_line_variants_pipeline_with_germlinev3.sh
#
# Purpose: Implement the recommended methodology for combining DeepVariant and 
# DeepSomatic calls for cancer cell lines without matched normal samples.
# The pipeline follows these key steps:
# 1. Extract GERMLINE-labeled variants from DeepSomatic
# 2. Use these positions to extract and keep DeepVariant calls
# 3. Extract PASS variants from DeepSomatic
# 4. Combine both sets while maintaining call information

set -e  # Exit on error
set -o pipefail  # Catch errors in pipelines

# Configuration
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
BASE_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/encode_analysis"
SAMPLE="K562"  # Can be parameterized for other cell lines
OUTDIR="${BASE_DIR}/data/variants/cell_line_analysis_v3/${SAMPLE}"
LOG_DIR="${BASE_DIR}/logs/${TIMESTAMP}"

# Input files
DEEPSOMATIC_VCF="${BASE_DIR}/data/variants/deepsomatic/${SAMPLE}_WGS_pair1.deepsomatic.with_germline.vcf.gz"
DEEPVARIANT_VCF="${BASE_DIR}/data/variants/deepvariant/${SAMPLE}_WGS_pair1.deepvariant.vcf.gz"


# Add after input files declaration:
# Validate input files
if [[ ! -f ${DEEPSOMATIC_VCF} || ! -f ${DEEPSOMATIC_VCF}.tbi ]]; then
    echo "ERROR: DeepSomatic VCF or index file not found"
    exit 1
fi

if [[ ! -f ${DEEPVARIANT_VCF} || ! -f ${DEEPVARIANT_VCF}.tbi ]]; then
    echo "ERROR: DeepVariant VCF or index file not found"
    exit 1
fi

# Set up logging
mkdir -p ${LOG_DIR}
exec 1> >(tee "${LOG_DIR}/pipeline_execution.log")
exec 2> >(tee "${LOG_DIR}/pipeline_errors.log")

# Function to log steps and timing
log_step() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
    echo "-------------------------------------------"
}
# Create output directory structure
mkdir -p ${OUTDIR}/{germline,somatic,final}
if [ -f ${OUTDIR}/germline/deepsomatic_germline.vcf.gz ]; then
    log_step "Archiving previous results"
    mkdir -p ${OUTDIR}/archive_${TIMESTAMP}
    mv ${OUTDIR}/germline/* ${OUTDIR}/archive_${TIMESTAMP}/
fi

# Step 1: Extract GERMLINE-labeled variants from DeepSomatic
log_step "Step 1: Extracting GERMLINE variants from DeepSomatic"
bcftools view \
    -i 'FILTER=="GERMLINE"' \
    ${DEEPSOMATIC_VCF} \
    -Oz -o ${OUTDIR}/germline/deepsomatic_germline_sites.vcf.gz
bcftools index -t ${OUTDIR}/germline/deepsomatic_germline_sites.vcf.gz

# Step 2: Extract DeepVariant calls at GERMLINE positions
log_step "Step 2: Extracting DeepVariant calls at GERMLINE positions"
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' \
    ${OUTDIR}/germline/deepsomatic_germline_sites.vcf.gz > \
    ${OUTDIR}/germline/germline_positions.tsv

bcftools view \
    -f PASS \
    -R ${OUTDIR}/germline/germline_positions.tsv \
    ${DEEPVARIANT_VCF} \
    -Oz -o ${OUTDIR}/germline/high_confidence_germline.vcf.gz
bcftools index -t ${OUTDIR}/germline/high_confidence_germline.vcf.gz

# Step 3: Extract PASS variants from DeepSomatic
log_step "Step 3: Extracting PASS variants from DeepSomatic"
bcftools view \
    -i 'FILTER=="PASS"' \
    ${DEEPSOMATIC_VCF} \
    -Oz -o ${OUTDIR}/somatic/deepsomatic_pass.vcf.gz
bcftools index -t ${OUTDIR}/somatic/deepsomatic_pass.vcf.gz

# Step 4: Combine the variants
log_step "Step 4: Combining variant sets"
bcftools concat \
    -a \
    -Oz \
    ${OUTDIR}/germline/high_confidence_germline.vcf.gz \
    ${OUTDIR}/somatic/deepsomatic_pass.vcf.gz \
    -o ${OUTDIR}/final/combined_variants.vcf.gz
bcftools index -t ${OUTDIR}/final/combined_variants.vcf.gz

# Generate summary statistics
log_step "Generating summary statistics"
{
    echo "=== Variant Analysis Summary ==="
    echo "Timestamp: ${TIMESTAMP}"
    echo "Sample: ${SAMPLE}"
    echo ""
    echo "DeepSomatic GERMLINE variants: $(bcftools view -H ${OUTDIR}/germline/deepsomatic_germline_sites.vcf.gz | wc -l)"
    echo "High-confidence germline variants (DeepVariant): $(bcftools view -H ${OUTDIR}/germline/high_confidence_germline.vcf.gz | wc -l)"
    echo "DeepSomatic PASS (somatic) variants: $(bcftools view -H ${OUTDIR}/somatic/deepsomatic_pass.vcf.gz | wc -l)"
    echo "Total variants in final set: $(bcftools view -H ${OUTDIR}/final/combined_variants.vcf.gz | wc -l)"
} > "${OUTDIR}/final/analysis_summary.txt"

log_step "Pipeline completed successfully"