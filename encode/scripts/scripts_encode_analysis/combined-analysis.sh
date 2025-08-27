#!/bin/bash

# Directories
OUTPUT_DIR="data/variants/comparisons_new"
IGV_DIR="$OUTPUT_DIR/igv_report"
REGION_DIR="$OUTPUT_DIR/regions"
STATS_DIR="$OUTPUT_DIR/stats"

mkdir -p "$IGV_DIR"
mkdir -p "$REGION_DIR"
mkdir -p "$STATS_DIR"

# Input files
DEEPSOMATIC="data/variants/deepsomatic/K562_WGS_pair1.deepsomatic.vcf.gz"
DEEPVARIANT="data/variants/deepvariant/K562_WGS_pair1.deepvariant.vcf.gz"
ENCODE_VCFS=(
    "data/variants/encode_lifted/ENCFF752OAX.hg38.sorted.vcf.gz"
    "data/variants/encode_lifted/ENCFF785JVR.hg38.sorted.vcf.gz"
    "data/variants/encode_lifted/ENCFF574MDJ.hg38.sorted.vcf.gz"
)

# Chromosomes to iterate through
CHROMS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 \
chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 \
chr21 chr22 chrX chrY chrM"

################################################################################
# 1. Compare DeepSomatic and DeepVariant
################################################################################
echo "Comparing DeepSomatic and DeepVariant..."
mkdir -p "$STATS_DIR/callers"
bcftools isec \
    -p "$STATS_DIR/callers" \
    "$DEEPSOMATIC" \
    "$DEEPVARIANT"

# Generate initial stats
{
    echo "Caller Comparison Statistics"
    echo "==========================="
    echo "Unique to DeepSomatic: $(grep -vc "^#" $STATS_DIR/callers/0000.vcf)"
    echo "Unique to DeepVariant: $(grep -vc "^#" $STATS_DIR/callers/0001.vcf)"
    echo "Shared variants: $(grep -vc "^#" $STATS_DIR/callers/0002.vcf)"
} > "$STATS_DIR/caller_comparison.txt"

################################################################################
# 2. Compare with ENCODE VCFs
################################################################################
echo "Comparing with ENCODE VCFs..."
mkdir -p "$STATS_DIR/encode"
mkdir -p "$STATS_DIR/encode_temp"

for encode_vcf in "${ENCODE_VCFS[@]}"; do
    base=$(basename "$encode_vcf" .hg38.sorted.vcf.gz)
    echo "Processing $base..."
    
    mkdir -p "$STATS_DIR/encode/$base"
    
    # Compare with DeepSomatic
    bcftools isec \
        -p "$STATS_DIR/encode/$base/deepsomatic" \
          -w1 \
          -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX \
        "$DEEPSOMATIC" \
        "$encode_vcf"
    
    # Compare with DeepVariant
    bcftools isec \
        -p "$STATS_DIR/encode/$base/deepvariant" \
        "$DEEPVARIANT" \
        "$encode_vcf"
        
    # Generate stats
    {
        echo "$base Comparison Results"
        echo "======================="
        echo "DeepSomatic comparison:"
        echo "  Unique to DeepSomatic: $(grep -vc "^#" $STATS_DIR/encode/$base/deepsomatic/0000.vcf)"
        echo "  Unique to ENCODE: $(grep -vc "^#" $STATS_DIR/encode/$base/deepsomatic/0001.vcf)"
        echo "  Shared variants: $(grep -vc "^#" $STATS_DIR/encode/$base/deepsomatic/0002.vcf)"
        echo
        echo "DeepVariant comparison:"
        echo "  Unique to DeepVariant: $(grep -vc "^#" $STATS_DIR/encode/$base/deepvariant/0000.vcf)"
        echo "  Unique to ENCODE: $(grep -vc "^#" $STATS_DIR/encode/$base/deepvariant/0001.vcf)"
        echo "  Shared variants: $(grep -vc "^#" $STATS_DIR/encode/$base/deepvariant/0002.vcf)"
    } > "$STATS_DIR/encode/$base.txt"
done

################################################################################
# 3. Extract regions by chromosome
################################################################################
# Temporary files
TMP_GERMLINE="$REGION_DIR/tmp_germline.bed"
TMP_HET_LIKE="$REGION_DIR/tmp_het_like.bed"
TMP_DV_HET="$REGION_DIR/tmp_dv_het.bed"
TMP_ENCODE="$REGION_DIR/tmp_encode.bed"

> "$TMP_GERMLINE"
> "$TMP_HET_LIKE"
> "$TMP_DV_HET"
> "$TMP_ENCODE"

# Extract regions for each chromosome
for CHR in $CHROMS; do
    # High-confidence GERMLINE (VAF >= 0.90)
    bcftools view -f PASS -r "$CHR" "$DEEPSOMATIC" | \
        bcftools query -f '%CHROM\t%POS\t%FILTER\t[%VAF]\n' | \
        awk '$4 >= 0.90' | shuf | head -n 100 | \
        awk -v OFS='\t' '{
            start=$2-100;
            if(start<0) start=0;
            print $1,start,$2+100,"GERMLINE_VAF_"$4,1,"+"
        }' >> "$TMP_GERMLINE"
    
    # DeepSomatic heterozygous-like
    bcftools view -f PASS -r "$CHR" "$DEEPSOMATIC" | \
        bcftools query -f '%CHROM\t%POS\t%FILTER\t[%VAF]\n' | \
        awk '$4 >= 0.45 && $4 <= 0.55' | shuf | head -n 100 | \
        awk -v OFS='\t' '{
            start=$2-100;
            if(start<0) start=0;
            print $1,start,$2+100,"HET_VAF_"$4,2,"+"
        }' >> "$TMP_HET_LIKE"
    
    # DeepVariant heterozygous
    bcftools view -r "$CHR" "$DEEPVARIANT" | \
        bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t[%AD]\n' -i 'GT="0/1"' | \
        shuf | head -n 100 | \
        awk -v OFS='\t' '{
            start=$2-100;
            if(start<0) start=0;
            print $1,start,$2+100,"DV_HET_"$3"_"$4"_AD"$6,3,"+"
        }' >> "$TMP_DV_HET"

    # ENCODE concordant regions
    for encode_vcf in "${ENCODE_VCFS[@]}"; do
        base=$(basename "$encode_vcf" .hg38.sorted.vcf.gz)
        comp_dir="$STATS_DIR/encode_temp/${base}"
        mkdir -p "$comp_dir"
        
        # Create temporary compressed files
        bcftools view -r "$CHR" "$DEEPSOMATIC" | bgzip > "$comp_dir/ds.vcf.gz"
        bcftools view -r "$CHR" "$encode_vcf" | bgzip > "$comp_dir/encode.vcf.gz"
        bcftools index "$comp_dir/ds.vcf.gz"
        bcftools index "$comp_dir/encode.vcf.gz"

        # Check content before intersection
        echo "Variants in $CHR for DeepSomatic:"
        bcftools view "$comp_dir/ds.vcf.gz" | wc -l
        echo "Variants in $CHR for $base:"
        bcftools view "$comp_dir/encode.vcf.gz" | wc -l

        bcftools isec \
            -p "$comp_dir" \
            "$comp_dir/ds.vcf.gz" \
            "$comp_dir/encode.vcf.gz"

        # Check intersection results
        echo "Concordant variants in $CHR between DeepSomatic and $base:"
        if [ -f "$comp_dir/0002.vcf" ]; then
            wc -l "$comp_dir/0002.vcf"
            bcftools view "$comp_dir/0002.vcf" | \
                bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' | \
                awk -v OFS='\t' -v base="$base" '{
                    start=$2-100;
                    if(start<0) start=0;
                    print $1,start,$2+100,"ENCODE_"base"_"NR,4,"+"
                }' | head -n 33 >> "$TMP_ENCODE"
        else
            echo "No concordant variants found"
        fi
                    
        rm -rf "$comp_dir" 2>/dev/null || true
    done


done

################################################################################
# 4. Create combined BED and IGV report
################################################################################
# Combine regions
{
    echo -e "#chrom\tstart\tend\tname\tscore\tstrand"
    cat "$TMP_GERMLINE"
    cat "$TMP_HET_LIKE"
    cat "$TMP_DV_HET"
    cat "$TMP_ENCODE"
} > "$REGION_DIR/combined_regions.bed"

# Create IGV report
echo "Creating IGV report..."
create_report \
    "$REGION_DIR/combined_regions.bed" \
    data/reference/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
    --flanking 1000 \
    --tracks \
        "$DEEPSOMATIC" \
        "$DEEPVARIANT" \
        "${ENCODE_VCFS[@]}" \
    --output "$IGV_DIR/variants_report.html"

# Generate final summary
{
    echo "Region Summary"
    echo "============="
    echo
    echo "High-confidence GERMLINE regions (VAF >= 0.90):"
    wc -l "$TMP_GERMLINE"
    echo
    echo "DeepSomatic heterozygous-like regions (0.45 <= VAF <= 0.55):"
    wc -l "$TMP_HET_LIKE"
    echo
    echo "DeepVariant heterozygous:"
    wc -l "$TMP_DV_HET"
    echo
    echo "ENCODE concordant regions:"
    wc -l "$TMP_ENCODE"
    echo
    echo "Total combined regions:"
    grep -v "^#" "$REGION_DIR/combined_regions.bed" | wc -l
    
    echo -e "\nChromosome distribution in combined BED:"
    cut -f1 "$REGION_DIR/combined_regions.bed" | grep -v "^#" | sort | uniq -c
} > "$REGION_DIR/region_summary.txt"

echo "Analysis complete!"
echo "IGV report created at $IGV_DIR/variants_report.html"
echo "Region summary available at $REGION_DIR/region_summary.txt"
echo "Detailed comparisons in $STATS_DIR/"