#!/bin/bash

# Define sample ID and input VCF
SAMPLE="GTEX-1117F"
VCF="/mnt/czi-sci-ai/intrinsic-variation-gene-ex-3/project_gene_regulation/data/GTEx/WGS/phg001796.v1.GTEx_v9_WGS_phased.genotype-calls-vcf.c1/GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv.SHAPEIT2_phased.vcf.gz"
OUTPUT_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/wgs/data/qc/variants/gtex_counts"

# Create output directory
mkdir -p $OUTPUT_DIR

# Initialize counters
total_snps=0
total_indels=0

# Process all standard chromosomes with correct chromosome format
for chr in {1..22} X; do
    echo "Processing chromosome chr${chr}..."
    
    # Count SNPs (use grep to filter non-reference genotypes)
    bcftools view -s $SAMPLE -r chr${chr} -v snps $VCF > $OUTPUT_DIR/temp_snps.vcf
    snp_count=$(grep -v "^#" $OUTPUT_DIR/temp_snps.vcf | grep -v "0|0" | grep -v "\\.|\\." | wc -l)
    echo "  SNPs: $snp_count"
    total_snps=$((total_snps + snp_count))
    
    # Count indels (use grep to filter non-reference genotypes)
    bcftools view -s $SAMPLE -r chr${chr} -v indels $VCF > $OUTPUT_DIR/temp_indels.vcf
    indel_count=$(grep -v "^#" $OUTPUT_DIR/temp_indels.vcf | grep -v "0|0" | grep -v "\\.|\\." | wc -l)
    echo "  Indels: $indel_count"
    total_indels=$((total_indels + indel_count))
    
    # Remove temporary files
    rm -f $OUTPUT_DIR/temp_snps.vcf $OUTPUT_DIR/temp_indels.vcf
    
    # Save chromosome-specific counts
    echo -e "Type\tCount" > $OUTPUT_DIR/${SAMPLE}_chr${chr}.txt
    echo -e "SNPs\t$snp_count" >> $OUTPUT_DIR/${SAMPLE}_chr${chr}.txt
    echo -e "Indels\t$indel_count" >> $OUTPUT_DIR/${SAMPLE}_chr${chr}.txt
done

# Save total counts
echo -e "Type\tCount" > $OUTPUT_DIR/${SAMPLE}_total.txt
echo -e "SNPs\t$total_snps" >> $OUTPUT_DIR/${SAMPLE}_total.txt
echo -e "Indels\t$total_indels" >> $OUTPUT_DIR/${SAMPLE}_total.txt

echo "Processing complete!"
echo "Total SNPs: $total_snps"
echo "Total Indels: $total_indels"