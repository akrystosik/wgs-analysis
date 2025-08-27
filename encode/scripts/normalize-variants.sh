#!/bin/bash

# Input files
ENCODE_VCF="data/variants/encode/ENCFF752OAX.vcf.gz"
DEEPSOMATIC="data/variants/deepsomatic/K562_WGS_pair1.deepsomatic.vcf.gz"
DEEPVARIANT="data/variants/deepvariant/K562_WGS_pair1.deepvariant.vcf.gz"
OUTPUT_DIR="k562_normalized_comparison"

mkdir -p $OUTPUT_DIR

echo "Normalizing variants..."

# 1. Normalize each VCF
for vcf in $ENCODE_VCF $DEEPSOMATIC $DEEPVARIANT; do
    base=$(basename $vcf .vcf.gz)
    echo "Processing $base..."
    
    # Split multi-allelic sites + left-align indels
    bcftools norm -m-any -f data/reference/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta $vcf | \
    bcftools view -i 'FILTER=="PASS"' | \
    bcftools norm -f data/reference/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta > \
    $OUTPUT_DIR/${base}_normalized.vcf

    # Extract variants to BED format for visualization
    bcftools query -f '%CHROM\t%POS\t%POS\t%REF>%ALT\n' $OUTPUT_DIR/${base}_normalized.vcf | \
    sort -k1,1 -k2,2n > $OUTPUT_DIR/${base}_variants.bed
done

# 2. Find overlapping regions of interest
bedtools intersect -wa -wb \
    -a $OUTPUT_DIR/K562_WGS_pair1.deepsomatic_variants.bed \
    -b $OUTPUT_DIR/K562_WGS_pair1.deepvariant_variants.bed | \
    cut -f1-4 | sort -u > $OUTPUT_DIR/caller_concordant.bed

# 3. Create regions file for IGV report focusing on:
# - Regions where both callers agree
# - Add 100bp padding
awk -v OFS='\t' '{
    start=$2-100;
    if(start<0) start=0;
    print $1,start,$3+100,$4
}' $OUTPUT_DIR/caller_concordant.bed > $OUTPUT_DIR/igv_regions.bed

# 4. Generate IGV report config
cat > $OUTPUT_DIR/igv_config.json <<EOL
{
    "reference": {
        "id": "hg38",
        "name": "Human (GRCh38/hg38)",
        "fastaURL": "data/reference/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta",
        "indexed": true
    },
    "tracks": [
        {
            "name": "DeepSomatic",
            "type": "variant",
            "format": "vcf",
            "url": "${DEEPSOMATIC}",
            "indexed": true
        },
        {
            "name": "DeepVariant",
            "type": "variant",
            "format": "vcf",
            "url": "${DEEPVARIANT}",
            "indexed": true
        },
        {
            "name": "ENCODE",
            "type": "variant",
            "format": "vcf",
            "url": "${ENCODE_VCF}",
            "indexed": true
        }
    ]
}
EOL

# 5. Generate summary statistics
{
    echo "Normalized Variant Analysis Summary"
    echo "=================================="
    echo
    echo "Variant counts after normalization:"
    for vcf in $OUTPUT_DIR/*_normalized.vcf; do
        echo "$(basename $vcf): $(grep -v '^#' $vcf | wc -l) variants"
    done
    echo
    echo "Concordant regions: $(wc -l < $OUTPUT_DIR/caller_concordant.bed)"
} > $OUTPUT_DIR/analysis_summary.txt

echo "Analysis complete! Results in $OUTPUT_DIR/"