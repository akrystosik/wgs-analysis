#!/bin/bash

echo "Starting K562 DeepVariant vs DeepSomatic Analysis"
echo "==============================================="

# Define input files
DEEPVARIANT="data/variants/deepvariant/K562_WGS_pair1.deepvariant.vcf.gz"
DEEPSOMATIC="data/variants/deepsomatic/K562_WGS_pair1.deepsomatic.vcf.gz"

# 1. Initial Comparison
echo "1. Performing Initial Site Comparison"
vcftools --gzvcf $DEEPVARIANT \
    --gzdiff $DEEPSOMATIC \
    --diff-site \
    --out k562_comparison

# 2. Generate Discordant Sites Analysis
echo "2. Analyzing Discordant Sites"
vcftools --gzvcf $DEEPVARIANT \
    --gzdiff $DEEPSOMATIC \
    --diff-site-discordance \
    --out k562_discordant_sites

# Extract regions and analyze discordant variants
cut -f1,2 k562_discordant_sites.diff.sites | tail -n +2 > k562_regions.txt

echo "3. Analyzing DeepVariant Discordant Calls"
bcftools view -R k562_regions.txt $DEEPVARIANT | \
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t%QUAL\n' > k562_discordant_deepvariant.txt

# 3. Analyze DeepSomatic PASS variants
echo "4. Analyzing DeepSomatic PASS Variants"
echo "4.1 Distribution of GT calls for PASS variants:"
bcftools view -f PASS $DEEPSOMATIC | \
    bcftools query -f '[%GT]\n' | sort | uniq -c

# 4. VAF Distribution for 1/1 calls
echo "4.2 VAF Distribution for 1/1 PASS variants:"
bcftools view -f PASS $DEEPSOMATIC | \
    bcftools query -f '[%GT]\t[%VAF]\n' | \
    awk '$1=="1/1" {if($2<0.1) a["0-0.1"]++; 
         else if($2<0.2) a["0.1-0.2"]++; 
         else if($2<0.3) a["0.2-0.3"]++; 
         else if($2<0.4) a["0.3-0.4"]++; 
         else if($2<0.5) a["0.4-0.5"]++; 
         else if($2<0.6) a["0.5-0.6"]++; 
         else if($2<0.7) a["0.6-0.7"]++; 
         else if($2<0.8) a["0.7-0.8"]++; 
         else if($2<0.9) a["0.8-0.9"]++; 
         else a["0.9-1.0"]++}
         END{for(i in a)print i,a[i]}' | sort

# 5. Analyze potential heterozygous variants (VAF 0.4-0.6)
echo "5. Analyzing Potential Heterozygous Variants"
echo "5.1 Extracting potential heterozygous sites (VAF 0.4-0.6):"
bcftools view -f PASS $DEEPSOMATIC | \
    bcftools query -f '%CHROM\t%POS\n' -i 'GT="1/1" && VAF >= 0.4 && VAF <= 0.6' > k562_potential_het_sites.txt

echo "5.2 DeepVariant genotypes at these sites:"
bcftools view -R k562_potential_het_sites.txt $DEEPVARIANT | \
    bcftools query -f '[%GT]\n' | sort | uniq -c

# 6. Coverage Analysis
echo "6. Coverage Distribution Analysis"
bcftools view -f PASS $DEEPSOMATIC | \
    bcftools query -f '[%AD]\n' | \
    awk -F',' '{sum=$1+$2; print sum}' | sort -n | \
    awk '{if($1<10) a["0-10"]++; 
         else if($1<20) a["10-20"]++; 
         else if($1<30) a["20-30"]++; 
         else if($1<40) a["30-40"]++; 
         else a["40+"]++}
         END{for(i in a)print i,a[i]}' | sort

echo "Analysis Complete!"