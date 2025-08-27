#!/bin/bash

OUTPUT_DIR="data/variants/three_way_analysis"
mkdir -p "$OUTPUT_DIR"
STATS_FILE="$OUTPUT_DIR/three_way_comparison.txt"

# Input files
DEEPSOMATIC="data/variants/deepsomatic/K562_WGS_pair1.deepsomatic.vcf.gz"
DEEPVARIANT="data/variants/deepvariant/K562_WGS_pair1.deepvariant.vcf.gz"
ENCODE_VCF="data/variants/encode_lifted/ENCFF752OAX.hg38.sorted.vcf.gz"

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

    # Two-way comparisons
    bcftools isec -p "$comp_dir/ds_dv" "$comp_dir/ds.vcf.gz" "$comp_dir/dv.vcf.gz"
    bcftools isec -p "$comp_dir/ds_encode" "$comp_dir/ds.vcf.gz" "$comp_dir/encode.vcf.gz"
    bcftools isec -p "$comp_dir/dv_encode" "$comp_dir/dv.vcf.gz" "$comp_dir/encode.vcf.gz"

    # Get stats
    total_ds=$(bcftools view "$comp_dir/ds.vcf.gz" | grep -vc "^#")
    total_dv=$(bcftools view "$comp_dir/dv.vcf.gz" | grep -vc "^#")
    total_encode=$(bcftools view "$comp_dir/encode.vcf.gz" | grep -vc "^#")
    
    ds_dv_shared=$(bcftools view "$comp_dir/ds_dv/0002.vcf" | grep -vc "^#")
    ds_encode_shared=$(bcftools view "$comp_dir/ds_encode/0002.vcf" | grep -vc "^#")
    dv_encode_shared=$(bcftools view "$comp_dir/dv_encode/0002.vcf" | grep -vc "^#")

    echo "Chr$chr Analysis:" >> "$STATS_FILE"
    echo "Total variants:" >> "$STATS_FILE"
    echo "  DeepSomatic: $total_ds" >> "$STATS_FILE"
    echo "  DeepVariant: $total_dv" >> "$STATS_FILE"
    echo "  ENCODE: $total_encode" >> "$STATS_FILE"
    echo "Shared variants:" >> "$STATS_FILE"
    echo "  DS-DV: $ds_dv_shared" >> "$STATS_FILE"
    echo "  DS-ENCODE: $ds_encode_shared" >> "$STATS_FILE"
    echo "  DV-ENCODE: $dv_encode_shared" >> "$STATS_FILE"
    echo "" >> "$STATS_FILE"

    rm -rf "$comp_dir"
done
