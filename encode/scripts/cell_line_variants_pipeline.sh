#!/bin/bash
set -e  # Exit on error

# Clear previous outputs
echo "=== Clearing previous outputs ==="
rm -rf data/variants/cell_line_analysis/K562/germline/*
rm -rf data/variants/cell_line_analysis/K562/somatic/*
rm -rf data/variants/cell_line_analysis/K562/final/*
mkdir -p data/variants/cell_line_analysis/K562/{germline,somatic,final}

# 1. Get high-confidence variants from DeepSomatic (PASS only)
echo "=== Processing DeepSomatic variants ==="
bcftools view \
    -i 'FILTER=="PASS" && FORMAT/VAF[0:0] >= 0.4' \
    data/variants/deepsomatic/K562_WGS_pair1.deepsomatic.vcf.gz \
    -Oz -o data/variants/cell_line_analysis/K562/germline/deepsomatic_germline.vcf.gz

# Index the file
bcftools index -t data/variants/cell_line_analysis/K562/germline/deepsomatic_germline.vcf.gz

# 2. Get PASS variants from DeepVariant
echo "=== Processing DeepVariant variants ==="
bcftools view \
    -f PASS \
    data/variants/deepvariant/K562_WGS_pair1.deepvariant.vcf.gz \
    -Oz -o data/variants/cell_line_analysis/K562/germline/deepvariant_pass.vcf.gz

# Index the file
bcftools index -t data/variants/cell_line_analysis/K562/germline/deepvariant_pass.vcf.gz

# 3. Prepare sites-only VCFs for intersection
echo "=== Preparing sites for intersection ==="
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' \
    data/variants/cell_line_analysis/K562/germline/deepvariant_pass.vcf.gz \
    > data/variants/cell_line_analysis/K562/germline/deepvariant.sites

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' \
    data/variants/cell_line_analysis/K562/germline/deepsomatic_germline.vcf.gz \
    > data/variants/cell_line_analysis/K562/germline/deepsomatic.sites

# 4. Find intersection using simple file comparison
echo "=== Finding intersection ==="
sort data/variants/cell_line_analysis/K562/germline/deepvariant.sites | uniq > data/variants/cell_line_analysis/K562/germline/deepvariant.sites.sorted
sort data/variants/cell_line_analysis/K562/germline/deepsomatic.sites | uniq > data/variants/cell_line_analysis/K562/germline/deepsomatic.sites.sorted

# Find common sites
comm -12 \
    data/variants/cell_line_analysis/K562/germline/deepvariant.sites.sorted \
    data/variants/cell_line_analysis/K562/germline/deepsomatic.sites.sorted \
    > data/variants/cell_line_analysis/K562/germline/common.sites

# 5. Extract variants at common sites from DeepVariant
echo "=== Creating final VCF ==="
bcftools view \
    -R data/variants/cell_line_analysis/K562/germline/common.sites \
    data/variants/cell_line_analysis/K562/germline/deepvariant_pass.vcf.gz \
    -Oz -o data/variants/cell_line_analysis/K562/germline/high_confidence_germline.vcf.gz

# Index the final file
bcftools index -t data/variants/cell_line_analysis/K562/germline/high_confidence_germline.vcf.gz

# 6. Print analysis results
echo "=== K562 Variant Analysis ==="
echo "DeepVariant PASS variants: $(bcftools view -H data/variants/cell_line_analysis/K562/germline/deepvariant_pass.vcf.gz | wc -l)"
echo "DeepSomatic high-VAF variants: $(bcftools view -H data/variants/cell_line_analysis/K562/germline/deepsomatic_germline.vcf.gz | wc -l)"
echo "Intersection variants: $(wc -l < data/variants/cell_line_analysis/K562/germline/common.sites)"

# 7. Show first few overlapping variants
echo -e "\n=== Sample of Overlapping Variants ==="
head -n 5 data/variants/cell_line_analysis/K562/germline/common.sites