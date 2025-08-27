#!/bin/bash
set -e

# Define key variables and paths
TIMESTAMP=$(date +%Y%m%d)
OUTDIR="data/variants/cell_line_analysis/K562"
DEEPSOMATIC_VCF="data/variants/deepsomatic/K562_WGS_pair1.deepsomatic.with_germline.vcf.gz"
DEEPVARIANT_VCF="data/variants/deepvariant/K562_WGS_pair1.deepvariant.vcf.gz"

# Set up directory structure and handle existing files
mkdir -p ${OUTDIR}/{germline,somatic,final}
if [ -f ${OUTDIR}/germline/deepsomatic_germline.vcf.gz ]; then
    mkdir -p ${OUTDIR}/archive_${TIMESTAMP}
    mv ${OUTDIR}/germline/* ${OUTDIR}/archive_${TIMESTAMP}/
fi

# Step 1: Extract variants from DeepSomatic (both GERMLINE and PASS filters)
echo "Extracting DeepSomatic variants..."
bcftools view \
    -i 'FILTER=="GERMLINE" || FILTER=="PASS"' \
    ${DEEPSOMATIC_VCF} \
    -Oz -o ${OUTDIR}/germline/deepsomatic_germline.vcf.gz
bcftools index -t ${OUTDIR}/germline/deepsomatic_germline.vcf.gz

# Step 2: Extract PASS variants from DeepVariant
echo "Extracting DeepVariant variants..."
bcftools view \
    -f PASS \
    ${DEEPVARIANT_VCF} \
    -Oz -o ${OUTDIR}/germline/deepvariant_pass.vcf.gz
bcftools index -t ${OUTDIR}/germline/deepvariant_pass.vcf.gz

# Step 3: Extract variant sites for comparison
echo "Extracting and sorting variant sites..."
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' ${OUTDIR}/germline/deepvariant_pass.vcf.gz \
    > ${OUTDIR}/germline/deepvariant.sites
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' ${OUTDIR}/germline/deepsomatic_germline.vcf.gz \
    > ${OUTDIR}/germline/deepsomatic.sites

# Step 4: Sort sites for comparison
sort ${OUTDIR}/germline/deepvariant.sites | uniq > ${OUTDIR}/germline/deepvariant.sites.sorted
sort ${OUTDIR}/germline/deepsomatic.sites | uniq > ${OUTDIR}/germline/deepsomatic.sites.sorted

# Step 5: Find overlapping sites
echo "Finding overlapping sites..."
comm -12 \
    ${OUTDIR}/germline/deepvariant.sites.sorted \
    ${OUTDIR}/germline/deepsomatic.sites.sorted \
    > ${OUTDIR}/germline/common.sites

# Step 6: Create final high-confidence VCF
echo "Creating final VCF..."
bcftools view \
    -R ${OUTDIR}/germline/common.sites \
    ${OUTDIR}/germline/deepvariant_pass.vcf.gz \
    -Oz -o ${OUTDIR}/germline/high_confidence_germline.vcf.gz
bcftools index -t ${OUTDIR}/germline/high_confidence_germline.vcf.gz

# Step 7: Find variants unique to DeepSomatic
echo "Analyzing unique DeepSomatic variants..."
bcftools isec \
    -C \
    ${OUTDIR}/germline/deepsomatic_germline.vcf.gz \
    ${OUTDIR}/germline/deepvariant_pass.vcf.gz \
    -p ${OUTDIR}/germline/unique_analysis

# Step 8: Print analysis results
echo "=== K562 Variant Analysis ==="
echo "DeepVariant PASS variants: $(bcftools view -H ${OUTDIR}/germline/deepvariant_pass.vcf.gz | wc -l)"
echo "DeepSomatic GERMLINE+PASS variants: $(bcftools view -H ${OUTDIR}/germline/deepsomatic_germline.vcf.gz | wc -l)"
echo "Intersection variants at matching sites: $(wc -l < ${OUTDIR}/germline/common.sites)"

# Step 9: Analyze DeepSomatic filter distribution
echo "=== DeepSomatic Filter Distribution ==="
bcftools view ${OUTDIR}/germline/unique_analysis/0000.vcf | grep -v "^#" | cut -f7 | sort | uniq -c

# Step 10: Create detailed analysis log
{
    echo "=== Analysis Run ${TIMESTAMP} ==="
    echo "DeepSomatic VCF: $(bcftools --version | head -n1)"
    
    echo -e "\nVersion Information:"
    bcftools view -h ${DEEPSOMATIC_VCF} | grep "##" | grep -i "version"
    
    echo -e "\nOverall DeepSomatic Filter Distribution:"
    bcftools view ${OUTDIR}/germline/deepsomatic_germline.vcf.gz | grep -v "^#" | cut -f7 | sort | uniq -c
    
    echo -e "\nFilter Distribution for Variants Unique to DeepSomatic:"
    bcftools view ${OUTDIR}/germline/unique_analysis/0000.vcf | grep -v "^#" | cut -f7 | sort | uniq -c
    
    echo -e "\nSite Overlap Summary:"
    echo "Total sites in overlap: $(wc -l < ${OUTDIR}/germline/common.sites)"
    echo -e "\nSample of overlapping sites:"
    head -n 5 ${OUTDIR}/germline/common.sites
} > ${OUTDIR}/germline/analysis_log_${TIMESTAMP}.txt

echo "Analysis complete. Results written to ${OUTDIR}/germline/analysis_log_${TIMESTAMP}.txt"