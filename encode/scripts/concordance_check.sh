#!/bin/bash

OUTPUT_DIR="data/variants/three_way_analysis"
mkdir -p "$OUTPUT_DIR"

DEEPSOMATIC="data/variants/deepsomatic/K562_WGS_pair1.deepsomatic.vcf.gz"
DEEPVARIANT="data/variants/deepvariant/K562_WGS_pair1.deepvariant.vcf.gz"
ENCODE_VCF="data/variants/encode_lifted/ENCFF752OAX.hg38.sorted.vcf.gz"
STATS_FILE="$OUTPUT_DIR/three_way_comparison.txt"

echo "K562 Three-Way Variant Comparison" > "$STATS_FILE"
echo "===============================" >> "$STATS_FILE"

for chr in {1..22} X; do
    echo "Analyzing chr$chr"
    comp_dir="$OUTPUT_DIR/chr$chr"
    mkdir -p "$comp_dir"
    
    bcftools view -r chr$chr "$DEEPSOMATIC" | bgzip > "$comp_dir/ds.vcf.gz"
    bcftools view -r chr$chr "$DEEPVARIANT" | bgzip > "$comp_dir/dv.vcf.gz"
    bcftools view -r chr$chr "$ENCODE_VCF" | bgzip > "$comp_dir/encode.vcf.gz"
    
    bcftools index "$comp_dir/ds.vcf.gz"
    bcftools index "$comp_dir/dv.vcf.gz"
    bcftools index "$comp_dir/encode.vcf.gz"

    bcftools isec -p "$comp_dir/isec" \
        "$comp_dir/ds.vcf.gz" \
        "$comp_dir/dv.vcf.gz" \
        "$comp_dir/encode.vcf.gz"

    echo "Checking isec output files:"
    ls -l "$comp_dir/isec/"       

    total_ds=$(bcftools view "$comp_dir/ds.vcf.gz" | grep -vc "^#")
    total_dv=$(bcftools view "$comp_dir/dv.vcf.gz" | grep -vc "^#")
    total_encode=$(bcftools view "$comp_dir/encode.vcf.gz" | grep -vc "^#")
    
    # Each file represents variants unique to that combination
    unique_ds=$(bcftools view "$comp_dir/isec/0000.vcf" | grep -vc "^#")
    unique_dv=$(bcftools view "$comp_dir/isec/0001.vcf" | grep -vc "^#")
    unique_encode=$(bcftools view "$comp_dir/isec/0002.vcf" | grep -vc "^#")
    shared_all=$(bcftools view "$comp_dir/isec/0003.vcf" | grep -vc "^#")

    echo "Chr$chr Analysis:" >> "$STATS_FILE"
    echo "Total variants:" >> "$STATS_FILE"
    echo "  DeepSomatic: $total_ds" >> "$STATS_FILE"
    echo "  DeepVariant: $total_dv" >> "$STATS_FILE"
    echo "  ENCODE: $total_encode" >> "$STATS_FILE"
    echo "Variant distribution:" >> "$STATS_FILE"
    echo "  Unique to DeepSomatic: $unique_ds" >> "$STATS_FILE"
    echo "  Unique to DeepVariant: $unique_dv" >> "$STATS_FILE"
    echo "  Unique to ENCODE: $unique_encode" >> "$STATS_FILE"
    echo "  Shared by all: $shared_all" >> "$STATS_FILE"
    echo "" >> "$STATS_FILE"

    rm -rf "$comp_dir"
done