#!/bin/bash

# Set up directories
OUTPUT_DIR="data/variants/comparisons_new"
mkdir -p $OUTPUT_DIR/debug

# Check VCF headers and formats
echo "Checking VCF formats..."
{
    echo "DeepSomatic VCF Format:"
    echo "======================"
    bcftools view -h data/variants/deepsomatic/K562_WGS_pair1.deepsomatic.vcf.gz | grep -E "^##FORMAT|^##INFO|^##FILTER"
    
    echo -e "\nDeepVariant VCF Format:"
    echo "====================="
    bcftools view -h data/variants/deepvariant/K562_WGS_pair1.deepvariant.vcf.gz | grep -E "^##FORMAT|^##INFO|^##FILTER"
    
    echo -e "\nChromosome distribution in DeepSomatic:"
    bcftools view data/variants/deepsomatic/K562_WGS_pair1.deepsomatic.vcf.gz | grep -v "^#" | cut -f1 | sort | uniq -c
    
    echo -e "\nChromosome distribution in DeepVariant:"
    bcftools view data/variants/deepvariant/K562_WGS_pair1.deepvariant.vcf.gz | grep -v "^#" | cut -f1 | sort | uniq -c
} > $OUTPUT_DIR/debug/vcf_format.txt

# Test VAF extraction
echo "Testing VAF extraction..."
{
    echo "Sample DeepSomatic variants with VAF:"
    bcftools view data/variants/deepsomatic/K562_WGS_pair1.deepsomatic.vcf.gz | \
        bcftools query -f '%CHROM\t%POS\t%FILTER\t[%VAF]\n' | head -n 5
    
    echo -e "\nSample GERMLINE variants:"
    bcftools view data/variants/deepsomatic/K562_WGS_pair1.deepsomatic.vcf.gz | \
        bcftools query -f '%CHROM\t%POS\t%FILTER\t[%VAF]\n' -i 'FILTER="GERMLINE"' | head -n 5
    
    echo -e "\nSample heterozygous variants from DeepVariant:"
    bcftools view data/variants/deepvariant/K562_WGS_pair1.deepvariant.vcf.gz | \
        bcftools query -f '%CHROM\t%POS\t[%GT]\n' -i 'GT="0/1"' | head -n 5
} > $OUTPUT_DIR/debug/variant_samples.txt

echo "Debug output saved to $OUTPUT_DIR/debug/"