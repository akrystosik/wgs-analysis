#!/bin/bash
# This script performs the core variant analysis accounting for tumor-only sample limitations:

# - Extracts high-confidence germline variants from DeepVariant
# - Identifies potential somatic variants unique to DeepSomatic
# - Validates both sets against ENCODE

# This script creates our foundational BED files by processing the raw VCFs and identifying different variant categories.


OUTPUT_DIR="data/variants/tumor_only_analysis"
mkdir -p "$OUTPUT_DIR"

# Input files
DEEPSOMATIC="data/variants/deepsomatic/K562_WGS_pair1.deepsomatic.vcf.gz"
DEEPVARIANT="data/variants/deepvariant/K562_WGS_pair1.deepvariant.vcf.gz"
ENCODE_VCF="data/variants/encode_lifted/ENCFF752OAX.hg38.sorted.vcf.gz"

# First identify high-confidence germline variants from DeepVariant
bcftools view -f PASS "$DEEPVARIANT" | \
    bcftools query -f '%CHROM\t%POS\t%POS\n' > "$OUTPUT_DIR/germline_confident.bed"

# Get all PASS variants from DeepSomatic
bcftools view -f PASS "$DEEPSOMATIC" | \
    bcftools query -f '%CHROM\t%POS\t%POS\n' > "$OUTPUT_DIR/deepsomatic_all.bed"

# Find variants unique to DeepSomatic (potential somatic)
bedtools subtract -a "$OUTPUT_DIR/deepsomatic_all.bed" \
    -b "$OUTPUT_DIR/germline_confident.bed" > "$OUTPUT_DIR/potential_somatic.bed"

# Compare both sets with ENCODE validation
bedtools intersect -a "$OUTPUT_DIR/germline_confident.bed" \
    -b "$ENCODE_VCF" > "$OUTPUT_DIR/germline_validated.bed"
bedtools intersect -a "$OUTPUT_DIR/potential_somatic.bed" \
    -b "$ENCODE_VCF" > "$OUTPUT_DIR/somatic_validated.bed"