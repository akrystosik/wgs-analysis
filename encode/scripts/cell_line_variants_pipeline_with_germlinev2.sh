#!/bin/bash
set -e

# Configuration
TIMESTAMP=$(date +%Y%m%d)
OUTDIR="data/variants/cell_line_analysis/K562"
DEEPSOMATIC_VCF="data/variants/deepsomatic/K562_WGS_pair1.deepsomatic.with_germline.vcf.gz"
DEEPVARIANT_VCF="data/variants/deepvariant/K562_WGS_pair1.deepvariant.vcf.gz"

# Create directories and archive old results
mkdir -p ${OUTDIR}/{germline,somatic,final}
if [ -f ${OUTDIR}/germline/deepsomatic_germline.vcf.gz ]; then
    mkdir -p ${OUTDIR}/archive_${TIMESTAMP}
    mv ${OUTDIR}/germline/* ${OUTDIR}/archive_${TIMESTAMP}/
fi

# Extract GERMLINE-only variants from DeepSomatic
echo "Extracting GERMLINE variants from DeepSomatic..."
bcftools view \
    -i 'FILTER=="GERMLINE"' \
    ${DEEPSOMATIC_VCF} \
    -Oz -o ${OUTDIR}/germline/deepsomatic_germline_only.vcf.gz
bcftools index -t ${OUTDIR}/germline/deepsomatic_germline_only.vcf.gz

# Extract PASS variants from DeepVariant
echo "Extracting PASS variants from DeepVariant..."
bcftools view \
    -f PASS \
    ${DEEPVARIANT_VCF} \
    -Oz -o ${OUTDIR}/germline/deepvariant_pass.vcf.gz
bcftools index -t ${OUTDIR}/germline/deepvariant_pass.vcf.gz

# Find overlapping germline sites
echo "Finding overlapping germline sites..."
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' ${OUTDIR}/germline/deepsomatic_germline_only.vcf.gz \
    | sort | uniq > ${OUTDIR}/germline/deepsomatic_germline.sites
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' ${OUTDIR}/germline/deepvariant_pass.vcf.gz \
    | sort | uniq > ${OUTDIR}/germline/deepvariant.sites
comm -12 \
    ${OUTDIR}/germline/deepsomatic_germline.sites \
    ${OUTDIR}/germline/deepvariant.sites \
    > ${OUTDIR}/germline/germline_overlap.sites

# Extract PASS-only variants from DeepSomatic
echo "Extracting PASS variants from DeepSomatic..."
bcftools view \
    -i 'FILTER=="PASS"' \
    ${DEEPSOMATIC_VCF} \
    -Oz -o ${OUTDIR}/somatic/deepsomatic_pass_only.vcf.gz
bcftools index -t ${OUTDIR}/somatic/deepsomatic_pass_only.vcf.gz

# Create high-confidence germline set
echo "Creating high-confidence germline set..."
bcftools view \
    -R ${OUTDIR}/germline/germline_overlap.sites \
    ${OUTDIR}/germline/deepvariant_pass.vcf.gz \
    -Oz -o ${OUTDIR}/final/high_confidence_germline.vcf.gz
bcftools index -t ${OUTDIR}/final/high_confidence_germline.vcf.gz

# Normalize sample names
echo "Normalizing sample names..."
bcftools reheader \
    -s <(echo "K562_WGS_pair1") \
    ${OUTDIR}/final/high_confidence_germline.vcf.gz \
    -o ${OUTDIR}/final/high_confidence_germline.normalized.vcf.gz
bcftools index -t ${OUTDIR}/final/high_confidence_germline.normalized.vcf.gz

# Combine the files
echo "Combining germline and somatic variants..."
bcftools concat \
    -a \
    -Oz \
    ${OUTDIR}/final/high_confidence_germline.normalized.vcf.gz \
    ${OUTDIR}/somatic/deepsomatic_pass_only.vcf.gz \
    -o ${OUTDIR}/final/combined_variants.vcf.gz
bcftools index -t ${OUTDIR}/final/combined_variants.vcf.gz

# Print analysis results
echo "=== K562 Variant Analysis ===" > ${OUTDIR}/final/analysis_summary.txt
echo "DeepSomatic GERMLINE variants: $(bcftools view -H ${OUTDIR}/germline/deepsomatic_germline_only.vcf.gz | wc -l)" >> ${OUTDIR}/final/analysis_summary.txt
echo "DeepSomatic PASS variants: $(bcftools view -H ${OUTDIR}/somatic/deepsomatic_pass_only.vcf.gz | wc -l)" >> ${OUTDIR}/final/analysis_summary.txt
echo "High-confidence germline variants: $(bcftools view -H ${OUTDIR}/final/high_confidence_germline.vcf.gz | wc -l)" >> ${OUTDIR}/final/analysis_summary.txt
echo "Total variants in final set: $(bcftools view -H ${OUTDIR}/final/combined_variants.vcf.gz | wc -l)" >> ${OUTDIR}/final/analysis_summary.txt

echo "Analysis complete. Check ${OUTDIR}/final/ for output files"