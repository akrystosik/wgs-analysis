#!/bin/bash

# Input files
DEEPSOMATIC="data/variants/deepsomatic/K562_WGS_pair1.deepsomatic.vcf.gz"
DEEPVARIANT="data/variants/deepvariant/K562_WGS_pair1.deepvariant.vcf.gz"
ENCODE_VCF="ENCFF752OAX.vcf.gz"  # Primary ENCODE VCF
OUTPUT_DIR="k562_analysis_next"

mkdir -p $OUTPUT_DIR

echo "Starting K562 detailed analysis..."

# 1. ENCODE VCF Comparison
echo "1. Comparing with ENCODE validated calls..."

# Extract high-confidence calls from our data
bcftools view -f PASS $DEEPSOMATIC | \
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%VAF]\n' -i 'VAF >= 0.8' \
    > $OUTPUT_DIR/our_high_conf.txt

# Compare with ENCODE calls
bcftools isec -p $OUTPUT_DIR/encode_comparison \
    $DEEPSOMATIC $ENCODE_VCF

# Generate summary of overlap
echo "ENCODE comparison summary:" > $OUTPUT_DIR/encode_summary.txt
wc -l $OUTPUT_DIR/encode_comparison/000* >> $OUTPUT_DIR/encode_summary.txt

# 2. Chromosome-specific Analysis
echo "2. Analyzing chromosome-specific patterns..."

# Get distribution by chromosome for each caller
echo "DeepSomatic chromosome distribution:" > $OUTPUT_DIR/chr_patterns.txt
bcftools view -f PASS $DEEPSOMATIC | \
    bcftools query -f '%CHROM\t[%VAF]\n' | \
    awk -F'\t' '{
        sum[$1] += 1
        vaf_sum[$1] += $2
    }
    END {
        print "Chrom\tCount\tAvg_VAF"
        for (chr in sum) 
            printf "%s\t%d\t%.3f\n", chr, sum[chr], vaf_sum[chr]/sum[chr]
    }' | sort -k1,1 >> $OUTPUT_DIR/chr_patterns.txt

echo -e "\nDeepVariant chromosome distribution:" >> $OUTPUT_DIR/chr_patterns.txt
bcftools view $DEEPVARIANT | \
    bcftools query -f '%CHROM\t[%GT]\n' | \
    awk -F'\t' '{count[$1]++; if($2=="0/1") het[$1]++; if($2=="1/1") hom[$1]++}
    END {
        print "Chrom\tTotal\tHet\tHom"
        for (chr in count) 
            printf "%s\t%d\t%d\t%d\n", chr, count[chr], het[chr], hom[chr]
    }' | sort -k1,1 >> $OUTPUT_DIR/chr_patterns.txt

# 3. Discordant Calls Analysis
echo "3. Investigating discordant calls..."

# Create regions where callers disagree
bcftools view -f PASS $DEEPSOMATIC | \
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%VAF]\n' -i 'VAF >= 0.4' \
    > $OUTPUT_DIR/deepsomatic_calls.txt

bcftools view $DEEPVARIANT | \
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' \
    > $OUTPUT_DIR/deepvariant_calls.txt

# Find discordant sites
awk 'NR==FNR{a[$1,$2]=$5;next} 
    ($1,$2) in a{
        if(($5=="0/1" && a[$1,$2]>=0.8) || 
           ($5=="1/1" && a[$1,$2]<0.8) ||
           ($5=="0/0" && a[$1,$2]>=0.4))
            print $1"\t"$2"\t"a[$1,$2]"\t"$5
    }' \
    $OUTPUT_DIR/deepsomatic_calls.txt \
    $OUTPUT_DIR/deepvariant_calls.txt > $OUTPUT_DIR/discordant_calls.txt

# Create BED file for discordant regions
awk -v OFS='\t' '{
    start=$2 - 100
    end=$2 + 100
    if(start < 0) start=0
    print $1, start, end, "DS_VAF=" $3 ";DV_GT=" $4
}' $OUTPUT_DIR/discordant_calls.txt | \
    sort -k1,1 -k2,2n > $OUTPUT_DIR/discordant_regions.bed

# Generate summary statistics
echo "Generating final summaries..."

# Discordant calls summary
echo "Discordant calls analysis:" > $OUTPUT_DIR/discordant_summary.txt
echo "Total discordant sites: $(wc -l < $OUTPUT_DIR/discordant_calls.txt)" >> $OUTPUT_DIR/discordant_summary.txt
echo -e "\nBreakdown by type:" >> $OUTPUT_DIR/discordant_summary.txt
awk '{print $3, $4}' $OUTPUT_DIR/discordant_calls.txt | \
    sort | uniq -c >> $OUTPUT_DIR/discordant_summary.txt

echo "Analysis complete! Check $OUTPUT_DIR for results"