#!/bin/bash

OUTPUT_DIR="data/variants/tumor_only_analysis"
# Check the depmap bed file format and first few lines
head -n 5 data/variants/tumor_only_analysis/depmap_variants.bed
# Check one of the chromosome-specific DepMap files
head -n 5 data/variants/tumor_only_analysis/igv_report/chr1_depmap.bed
# Check if DepMap entries made it into the combined regions
grep "DEPMAP" data/variants/tumor_only_analysis/igv_report/regions_of_interest.bed | head -n 5
cat data/variants/tumor_only_analysis/igv_report/visualization_summary.txt

IGV_DIR="$OUTPUT_DIR/igv_report"
mkdir -p "$IGV_DIR"

# Define input VCFs and BAM
DEEPSOMATIC="data/variants/deepsomatic/K562_WGS_pair1.deepsomatic.vcf.gz"
DEEPVARIANT="data/variants/deepvariant/K562_WGS_pair1.deepvariant.vcf.gz"
ENCODE_VCF="data/variants/encode_lifted/ENCFF752OAX.hg38.sorted.vcf.gz"
BAM_FILE="data/bam/K562/K562_WGS_pair1.marked_duplicates.bam"

# DepMap BED created via depmap_to_bed.sh
DEP_MAP_BED="$OUTPUT_DIR/depmap_variants.bed"

# Per-chromosome region files
for chr in {1..22} X; do
    echo "Processing chr${chr}..."

    # 1) Validated Germline
    grep "^chr${chr}[[:space:]]" "$OUTPUT_DIR/germline_validated.bed" | \
      awk -v chr=$chr '{print $1"\t"$2-100"\t"$2+100"\tVALIDATED_GERMLINE_"NR"\t1\t+"}' | \
      head -n 50 > "$IGV_DIR/chr${chr}_germline.bed"

    # 2) Validated Somatic
    grep "^chr${chr}[[:space:]]" "$OUTPUT_DIR/somatic_validated.bed" | \
      awk -v chr=$chr '{print $1"\t"$2-100"\t"$2+100"\tVALIDATED_SOMATIC_"NR"\t2\t+"}' | \
      head -n 50 > "$IGV_DIR/chr${chr}_somatic.bed"

    # 3) Discordant Calls
    grep "^chr${chr}[[:space:]]" "$OUTPUT_DIR/potential_somatic.bed" | \
      awk -v chr=$chr '{print $1"\t"$2-100"\t"$2+100"\tDISCORDANT_"NR"\t3\t+"}' | \
      head -n 50 > "$IGV_DIR/chr${chr}_discordant.bed"

  # 4) DepMap Variants - Show ALL variants, not just first 50
      grep "^chr${chr}[[:space:]]" "$DEP_MAP_BED" | \
        awk -v chr=$chr '{start=$2; end=$3; print $1"\t"(start-100)"\t"(end+100)"\tDEPMAP_"NR"\t4\t+"}' \
        > "$IGV_DIR/chr${chr}_depmap.bed"      
done

# Create DepMap-specific regions file
cat "$IGV_DIR"/chr*_depmap.bed | sort -k1,1 -k2,2n > "$IGV_DIR/depmap_regions.bed"

# Combine all regions
cat "$IGV_DIR"/chr*_*.bed | sort -k1,1 -k2,2n > "$IGV_DIR/regions_of_interest.bed"

# Create IGV report
create_report \
    "$IGV_DIR/depmap_regions.bed" \
    data/reference/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
    --flanking 1000 \
    --tracks \
        "$BAM_FILE" \
        "$DEEPVARIANT" \
        "$DEEPSOMATIC" \
        "$ENCODE_VCF" \
        "$DEP_MAP_BED" \
    --output "$IGV_DIR/depmap_regions_variants_report.html"

# Generate summary statistics
echo "IGV Visualization Region Summary" > "$IGV_DIR/visualization_summary.txt"
echo "=============================" >> "$IGV_DIR/visualization_summary.txt"
echo "" >> "$IGV_DIR/visualization_summary.txt"

echo "Regions selected per category:" >> "$IGV_DIR/visualization_summary.txt"
echo "  Validated germline: $(grep -c VALIDATED_GERMLINE "$IGV_DIR/regions_of_interest.bed")" >> "$IGV_DIR/visualization_summary.txt"
echo "  Validated somatic: $(grep -c VALIDATED_SOMATIC "$IGV_DIR/regions_of_interest.bed")" >> "$IGV_DIR/visualization_summary.txt"
echo "  Discordant calls:  $(grep -c DISCORDANT       "$IGV_DIR/regions_of_interest.bed")" >> "$IGV_DIR/visualization_summary.txt"
echo "  DepMap calls:      $(grep -c DEPMAP           "$IGV_DIR/regions_of_interest.bed")" >> "$IGV_DIR/visualization_summary.txt"
echo "" >> "$IGV_DIR/visualization_summary.txt"

echo "Distribution by chromosome:" >> "$IGV_DIR/visualization_summary.txt"
cut -f1 "$IGV_DIR/regions_of_interest.bed" | sort | uniq -c >> "$IGV_DIR/visualization_summary.txt"

echo "IGV report generated at: $IGV_DIR/variants_report.html"
