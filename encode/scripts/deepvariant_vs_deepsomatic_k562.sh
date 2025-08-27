#!/bin/bash

# Paths and file definitions
PAIR_ID="K562_WGS_pair1"
DEEPSOMATIC_VCF="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/encode_analysis/data/variants/deepsomatic/${PAIR_ID}.deepsomatic.vcf.gz"
DEEPVARIANT_VCF="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/encode_analysis/data/variants/deepvariant/${PAIR_ID}.deepvariant.vcf.gz"
OUTPUT_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/encode_analysis/data/variants/analysis"
POTENTIAL_HET_SITES="${OUTPUT_DIR}/${PAIR_ID}_potential_het_sites.txt"
DEEPSOMATIC_FILTERED="${OUTPUT_DIR}/${PAIR_ID}_deepsomatic_filtered.vcf.gz"
DEEPVARIANT_FILTERED="${OUTPUT_DIR}/${PAIR_ID}_deepvariant_filtered.vcf.gz"
DISCORDANT_REGIONS="${OUTPUT_DIR}/${PAIR_ID}_discordant_regions.txt"
DISCORDANT_DETAILS="${OUTPUT_DIR}/${PAIR_ID}_discordant_details.txt"

mkdir -p ${OUTPUT_DIR}

# Step 1: Analyze VAF Distribution in DeepSomatic
echo "Analyzing VAF distribution for DeepSomatic..."
bcftools view -f PASS ${DEEPSOMATIC_VCF} | \
bcftools query -f '[%VAF]\n' -i 'VAF > 0' | \
awk '{if($1<0.1) a["0-0.1"]++; else if($1<0.2) a["0.1-0.2"]++; else if($1<0.3) a["0.2-0.3"]++; else if($1<0.4) a["0.3-0.4"]++; else if($1<0.5) a["0.4-0.5"]++; else if($1<0.6) a["0.5-0.6"]++; else if($1<0.7) a["0.6-0.7"]++; else if($1<0.8) a["0.7-0.8"]++; else if($1<0.9) a["0.8-0.9"]++; else a["0.9-1.0"]++} END {for(i in a) print i, a[i]}' > ${OUTPUT_DIR}/${PAIR_ID}_deepsomatic_vaf_distribution.txt

# Step 2: Filter High-Confidence Variants
echo "Filtering high-confidence variants..."
bcftools filter -i 'FILTER="PASS" && DP >= 10 && QUAL >= 30' ${DEEPSOMATIC_VCF} -Oz -o ${DEEPSOMATIC_FILTERED}
bcftools index ${DEEPSOMATIC_FILTERED}

bcftools filter -i 'FILTER="PASS" && DP >= 10 && QUAL >= 30' ${DEEPVARIANT_VCF} -Oz -o ${DEEPVARIANT_FILTERED}
bcftools index ${DEEPVARIANT_FILTERED}

# Step 3: Identify Discordant Regions
echo "Identifying discordant regions between DeepVariant and DeepSomatic..."
vcftools --gzvcf ${DEEPVARIANT_FILTERED} \
         --gzdiff ${DEEPSOMATIC_FILTERED} \
         --diff-site --out ${OUTPUT_DIR}/${PAIR_ID}_discordant

cut -f1,2 ${OUTPUT_DIR}/${PAIR_ID}_discordant.diff.sites | tail -n +2 > ${DISCORDANT_REGIONS}

# Step 4: Detailed Analysis of Discordant Regions
echo "Analyzing discordant regions in detail..."
bcftools view -R ${DISCORDANT_REGIONS} ${DEEPSOMATIC_VCF} | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t[%AD]\t[%VAF]\n' > ${OUTPUT_DIR}/${PAIR_ID}_deepsomatic_discordant.txt

bcftools view -R ${DISCORDANT_REGIONS} ${DEEPVARIANT_VCF} | \
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t[%AD]\t%QUAL\n' > ${OUTPUT_DIR}/${PAIR_ID}_deepvariant_discordant.txt

# Step 5: IGV Region File Generation
echo "Preparing IGV region file..."
awk '{print $1":"$2"-"$2}' ${DISCORDANT_REGIONS} > ${OUTPUT_DIR}/${PAIR_ID}_igv_regions.txt

# Step 6: Extract Example Variants for Visualization
echo "Extracting example variants for visualization..."
head -n 5 ${DISCORDANT_REGIONS} | while read -r line; do
    echo "DeepSomatic:"
    grep -w "$line" ${OUTPUT_DIR}/${PAIR_ID}_deepsomatic_discordant.txt
    echo "DeepVariant:"
    grep -w "$line" ${OUTPUT_DIR}/${PAIR_ID}_deepvariant_discordant.txt
done > ${DISCORDANT_DETAILS}

# Step 7: Summarize Results
echo "Analysis complete. Summary:"
echo "1. VAF distribution saved to: ${OUTPUT_DIR}/${PAIR_ID}_deepsomatic_vaf_distribution.txt"
echo "2. Filtered DeepSomatic VCF: ${DEEPSOMATIC_FILTERED}"
echo "3. Filtered DeepVariant VCF: ${DEEPVARIANT_FILTERED}"
echo "4. Discordant regions: ${DISCORDANT_REGIONS}"
echo "5. Discordant details: ${DISCORDANT_DETAILS}"
echo "6. IGV regions file: ${OUTPUT_DIR}/${PAIR_ID}_igv_regions.txt"



# Step 8: Generate IGV Report
echo -e "chr1\t1000000\t2000000" > regions.bed
echo "Generating IGV Report for filtered DeepSomatic and deepvariant VCF..."

create_report /mnt/czi-sci-ai/intrinsic-variation-gene-ex/encode_analysis/data/variants/analysis/K562_WGS_pair1_deepsomatic_filtered.vcf.gz \
  --genome hg38 \
  --output /mnt/czi-sci-ai/intrinsic-variation-gene-ex/encode_analysis/data/variants/analysis/K562_WGS_pair1_igv_report.html \
  --tracks /mnt/czi-sci-ai/intrinsic-variation-gene-ex/encode_analysis/data/bam/K562/K562_WGS_pair1.marked_duplicates.bam \
  /mnt/czi-sci-ai/intrinsic-variation-gene-ex/encode_analysis/data/variants/analysis/K562_WGS_pair1_deepvariant_filtered.vcf.gz \
  /mnt/czi-sci-ai/intrinsic-variation-gene-ex/encode_analysis/data/variants/analysis/K562_WGS_pair1_deepsomatic_filtered.vcf.gz \
  --bed /mnt/czi-sci-ai/intrinsic-variation-gene-ex/encode_analysis/data/variants/analysis/regions.bed
