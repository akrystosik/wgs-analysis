#!/bin/bash

OUTPUT_DIR="data/variants/detailed_concordance"
mkdir -p "$OUTPUT_DIR"
STATS_FILE="$OUTPUT_DIR/concordance_summary.txt"

DEEPSOMATIC="data/variants/deepsomatic/K562_WGS_pair1.deepsomatic.vcf.gz"
DEEPVARIANT="data/variants/deepvariant/K562_WGS_pair1.deepvariant.vcf.gz"
ENCODE_VCF="data/variants/encode_lifted/ENCFF752OAX.hg38.sorted.vcf.gz"

echo "Three-Way Concordance Summary" > "$STATS_FILE"
echo "============================" >> "$STATS_FILE"
echo "" >> "$STATS_FILE"

total_ds=0
total_dv=0
total_encode=0
total_shared=0

for chr in {1..22} X; do
    echo "Processing chr$chr"
    
    bcftools view -r chr$chr "$DEEPSOMATIC" | \
        bcftools filter -i 'VAF >= 0.45 && VAF <= 0.55 || VAF >= 0.90' | \
        bcftools query -f '%CHROM\t%POS\t%POS\n' > "$OUTPUT_DIR/ds.bed"

    bcftools view -r chr$chr "$DEEPVARIANT" | \
        bcftools filter -i 'GT="0/1" || GT="1/1"' | \
        bcftools query -f '%CHROM\t%POS\t%POS\n' > "$OUTPUT_DIR/dv.bed"

    bcftools view -r chr$chr "$ENCODE_VCF" | \
        bcftools query -f '%CHROM\t%POS\t%POS\n' > "$OUTPUT_DIR/encode.bed"

    bedtools intersect -a "$OUTPUT_DIR/ds.bed" -b "$OUTPUT_DIR/dv.bed" > "$OUTPUT_DIR/ds_dv.bed"
    bedtools intersect -a "$OUTPUT_DIR/ds_dv.bed" -b "$OUTPUT_DIR/encode.bed" > "$OUTPUT_DIR/shared.bed"

    ds_count=$(wc -l < "$OUTPUT_DIR/ds.bed")
    dv_count=$(wc -l < "$OUTPUT_DIR/dv.bed")
    encode_count=$(wc -l < "$OUTPUT_DIR/encode.bed")
    shared_count=$(wc -l < "$OUTPUT_DIR/shared.bed")

    total_ds=$((total_ds + ds_count))
    total_dv=$((total_dv + dv_count))
    total_encode=$((total_encode + encode_count))
    total_shared=$((total_shared + shared_count))

    overlap_pct=$(awk "BEGIN {printf \"%.2f\", ($shared_count * 100 / $encode_count)}")

    {
        echo "Chr$chr:"
        echo "  DeepSomatic: $ds_count"
        echo "  DeepVariant: $dv_count"
        echo "  ENCODE: $encode_count"
        echo "  Shared: $shared_count"
        echo "  Overlap %: ${overlap_pct}%"
        echo ""
    } >> "$STATS_FILE"
done

overall_pct=$(awk "BEGIN {printf \"%.2f\", ($total_shared * 100 / $total_encode)}")

{
    echo "Totals:"
    echo "  DeepSomatic: $total_ds"
    echo "  DeepVariant: $total_dv"
    echo "  ENCODE: $total_encode"
    echo "  Shared: $total_shared"
    echo "  Overall Overlap: ${overall_pct}%"
} >> "$STATS_FILE"