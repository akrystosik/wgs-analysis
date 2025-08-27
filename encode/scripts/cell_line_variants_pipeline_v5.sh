#!/bin/bash
set -e
set -o pipefail

CELL_LINE="$1"
BASE_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/encode_analysis"
VARIANTS_DIR="${BASE_DIR}/data/variants"
OUTDIR="${VARIANTS_DIR}/cell_line_analysis_v5/${CELL_LINE}"

mkdir -p "${OUTDIR}"/{intermediate,final}

DEEPVARIANT_VCF="${VARIANTS_DIR}/deepvariant/${CELL_LINE}_WGS_pair1.deepvariant.vcf.gz"
DEEPSOMATIC_VCF="${VARIANTS_DIR}/deepsomatic/${CELL_LINE}_WGS_pair1.deepsomatic.vcf.gz"

echo "Starting variant combination process for ${CELL_LINE}..."

# Create our header with the VARIANT_TYPE field definition
cat > "${OUTDIR}/intermediate/header.txt" << EOL
##fileformat=VCFv4.2
##INFO=<ID=VARIANT_TYPE,Number=1,Type=String,Description="Origin of variant call: somatic from DeepSomatic or germline from DeepVariant">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FILTER=<ID=PASS,Description="All filters passed">
EOL

# Process DeepSomatic variants and add somatic label
echo "Processing DeepSomatic variants..."
# First create a tab-delimited file with our annotations
bcftools view -f PASS "${DEEPSOMATIC_VCF}" | \
    bcftools query -f '%CHROM\t%POS\tsomatic\n' > "${OUTDIR}/intermediate/somatic_annotations.txt"

# Sort and compress for indexing
sort -k1,1 -k2,2n "${OUTDIR}/intermediate/somatic_annotations.txt" | \
    bgzip > "${OUTDIR}/intermediate/somatic_annotations.txt.gz"

# Index the annotations file
tabix -s 1 -b 2 -e 2 "${OUTDIR}/intermediate/somatic_annotations.txt.gz"

# Now annotate the VCF with our somatic labels
bcftools view -f PASS "${DEEPSOMATIC_VCF}" | \
    bcftools annotate \
    --header-lines "${OUTDIR}/intermediate/header.txt" \
    -a "${OUTDIR}/intermediate/somatic_annotations.txt.gz" \
    -c CHROM,POS,INFO/VARIANT_TYPE \
    --set-id +'%CHROM\_%POS\_%REF\_%ALT' \
    -Oz -o "${OUTDIR}/intermediate/deepsomatic_pass.vcf.gz"

bcftools index -t "${OUTDIR}/intermediate/deepsomatic_pass.vcf.gz"

# Create position file for filtering
echo "Identifying somatic variants for filtering..."
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' \
    "${OUTDIR}/intermediate/deepsomatic_pass.vcf.gz" > \
    "${OUTDIR}/intermediate/somatic_positions.tsv"

# Process DeepVariant calls with germline label
echo "Processing DeepVariant calls..."
# Create germline annotations
bcftools view -f PASS "${DEEPVARIANT_VCF}" | \
    bcftools view -T ^"${OUTDIR}/intermediate/somatic_positions.tsv" | \
    bcftools query -f '%CHROM\t%POS\tgermline\n' > "${OUTDIR}/intermediate/germline_annotations.txt"

# Sort and compress
sort -k1,1 -k2,2n "${OUTDIR}/intermediate/germline_annotations.txt" | \
    bgzip > "${OUTDIR}/intermediate/germline_annotations.txt.gz"

# Index
tabix -s 1 -b 2 -e 2 "${OUTDIR}/intermediate/germline_annotations.txt.gz"

# Annotate DeepVariant calls
bcftools view -f PASS "${DEEPVARIANT_VCF}" | \
    bcftools view -T ^"${OUTDIR}/intermediate/somatic_positions.tsv" | \
    bcftools annotate \
    --header-lines "${OUTDIR}/intermediate/header.txt" \
    -a "${OUTDIR}/intermediate/germline_annotations.txt.gz" \
    -c CHROM,POS,INFO/VARIANT_TYPE \
    --set-id +'%CHROM\_%POS\_%REF\_%ALT' \
    -Oz -o "${OUTDIR}/intermediate/deepvariant_unique.vcf.gz"

bcftools index -t "${OUTDIR}/intermediate/deepvariant_unique.vcf.gz"

# Standardize sample names
echo "Standardizing sample names..."
bcftools reheader \
    -s <(echo "${CELL_LINE}") \
    "${OUTDIR}/intermediate/deepvariant_unique.vcf.gz" \
    -o "${OUTDIR}/intermediate/deepvariant_unique_renamed.vcf.gz"

bcftools reheader \
    -s <(echo "${CELL_LINE}") \
    "${OUTDIR}/intermediate/deepsomatic_pass.vcf.gz" \
    -o "${OUTDIR}/intermediate/deepsomatic_pass_renamed.vcf.gz"

bcftools index -t "${OUTDIR}/intermediate/deepvariant_unique_renamed.vcf.gz"
bcftools index -t "${OUTDIR}/intermediate/deepsomatic_pass_renamed.vcf.gz"

# Combine variants
echo "Creating final combined variant set..."
bcftools concat \
    -a \
    "${OUTDIR}/intermediate/deepvariant_unique_renamed.vcf.gz" \
    "${OUTDIR}/intermediate/deepsomatic_pass_renamed.vcf.gz" \
    -Oz -o "${OUTDIR}/final/combined_variants.vcf.gz"

bcftools index -t "${OUTDIR}/final/combined_variants.vcf.gz"

# Generate summary statistics
echo "Generating analysis summary..."
{
    echo "=== Variant Analysis Summary v5 ==="
    echo "Generated: $(date '+%Y-%m-%d %H:%M:%S')"
    echo
    echo "Cell Line: ${CELL_LINE}"
    echo
    echo "1. Input Processing:"
    echo "   - DeepVariant PASS variants: $(bcftools view -H "${DEEPVARIANT_VCF}" -f PASS | wc -l)"
    echo "   - DeepSomatic somatic (PASS) variants: $(bcftools view -H "${OUTDIR}/intermediate/deepsomatic_pass.vcf.gz" | wc -l)"
    echo "   - DeepSomatic germline-labeled variants: $(bcftools view -H "${OUTDIR}/intermediate/deepsomatic_germline.vcf.gz" | wc -l)"
    echo
    echo "2. Final Set Composition:"
    echo "   - Non-overlapping DeepVariant calls: $(bcftools view -H "${OUTDIR}/intermediate/deepvariant_unique.vcf.gz" | wc -l)"
    echo "   - DeepSomatic PASS calls: $(bcftools view -H "${OUTDIR}/intermediate/deepsomatic_pass.vcf.gz" | wc -l)"
    echo "   - Total variants in combined set: $(bcftools view -H "${OUTDIR}/final/combined_variants.vcf.gz" | wc -l)"

    echo "Variant Type Verification:"
    echo "- Somatic variants in final VCF: $(bcftools view -H -i 'INFO/VARIANT_TYPE=="somatic"' "${OUTDIR}/final/combined_variants.vcf.gz" | wc -l)"
    echo "- Germline variants in final VCF: $(bcftools view -H -i 'INFO/VARIANT_TYPE=="germline"' "${OUTDIR}/final/combined_variants.vcf.gz" | wc -l)"

} > "${OUTDIR}/final/analysis_summary.txt"

echo "Variant combination complete. Results in ${OUTDIR}/final/"
echo "See ${OUTDIR}/final/analysis_summary.txt for detailed statistics"