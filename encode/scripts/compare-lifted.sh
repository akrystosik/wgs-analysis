#!/bin/bash

# Exit on error
set -e

# Input files
DEEPSOMATIC="data/variants/deepsomatic/K562_WGS_pair1.deepsomatic.vcf.gz"
DEEPVARIANT="data/variants/deepvariant/K562_WGS_pair1.deepvariant.vcf.gz"
ENCODE_VCF="data/variants/encode_lifted/ENCFF752OAX.hg38.sorted.vcf.gz"

# Output directory - use new path to avoid NFS locks
OUTPUT_DIR="data/variants/comparisons_new"
mkdir -p $OUTPUT_DIR

echo "Starting step-by-step comparison..."

# Step 1: Compare DeepSomatic and DeepVariant
echo "Step 1: Comparing DeepSomatic and DeepVariant..."
mkdir -p $OUTPUT_DIR/callers
bcftools isec \
    -p $OUTPUT_DIR/callers \
    $DEEPSOMATIC \
    $DEEPVARIANT

# Verify Step 1
echo "Verifying Step 1 output..."
ls -l $OUTPUT_DIR/callers/

# Generate initial stats
echo "Generating initial comparison stats..."
{
    echo "Initial Comparison: DeepSomatic vs DeepVariant"
    echo "----------------------------------------"
    echo "Unique to DeepSomatic: $(grep -vc "^#" $OUTPUT_DIR/callers/0000.vcf)"
    echo "Unique to DeepVariant: $(grep -vc "^#" $OUTPUT_DIR/callers/0001.vcf)"
    echo "Shared variants: $(grep -vc "^#" $OUTPUT_DIR/callers/0002.vcf)"
} > $OUTPUT_DIR/initial_stats.txt

# Step 2: Compare with ENCODE
echo "Step 2: Comparing with ENCODE..."
mkdir -p $OUTPUT_DIR/encode_comparison/{deepsomatic,deepvariant}

# Compare DeepSomatic with ENCODE
echo "Comparing DeepSomatic with ENCODE..."
bcftools isec \
    -p $OUTPUT_DIR/encode_comparison/deepsomatic \
    $DEEPSOMATIC \
    $ENCODE_VCF

# Compare DeepVariant with ENCODE
echo "Comparing DeepVariant with ENCODE..."
bcftools isec \
    -p $OUTPUT_DIR/encode_comparison/deepvariant \
    $DEEPVARIANT \
    $ENCODE_VCF

# Generate ENCODE comparison stats
echo "Generating ENCODE comparison stats..."
{
    echo -e "\nENCODE Comparison Results"
    echo "-------------------------"
    echo "DeepSomatic vs ENCODE:"
    echo "  Unique to DeepSomatic: $(grep -vc "^#" $OUTPUT_DIR/encode_comparison/deepsomatic/0000.vcf)"
    echo "  Unique to ENCODE: $(grep -vc "^#" $OUTPUT_DIR/encode_comparison/deepsomatic/0001.vcf)"
    echo "  Shared variants: $(grep -vc "^#" $OUTPUT_DIR/encode_comparison/deepsomatic/0002.vcf)"
    echo
    echo "DeepVariant vs ENCODE:"
    echo "  Unique to DeepVariant: $(grep -vc "^#" $OUTPUT_DIR/encode_comparison/deepvariant/0000.vcf)"
    echo "  Unique to ENCODE: $(grep -vc "^#" $OUTPUT_DIR/encode_comparison/deepvariant/0001.vcf)"
    echo "  Shared variants: $(grep -vc "^#" $OUTPUT_DIR/encode_comparison/deepvariant/0002.vcf)"
} >> $OUTPUT_DIR/initial_stats.txt

# Step 3: Create visualization regions
echo "Creating visualization regions..."
{
    echo -e "#chrom\tstart\tend\tname\tscore\tstrand"
    
    # Add regions where DeepSomatic and ENCODE agree
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' \
        $OUTPUT_DIR/encode_comparison/deepsomatic/0002.vcf | \
        head -n 500 | \
        awk -v OFS='\t' '{
            start=$2-100;
            if(start<0) start=0;
            print $1,start,$2+100,"DS_ENCODE_" NR,1,"+"
        }'
    
    # Add regions where DeepVariant and ENCODE agree
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' \
        $OUTPUT_DIR/encode_comparison/deepvariant/0002.vcf | \
        head -n 500 | \
        awk -v OFS='\t' '{
            start=$2-100;
            if(start<0) start=0;
            print $1,start,$2+100,"DV_ENCODE_" NR,2,"+"
        }'
} > $OUTPUT_DIR/visualization_regions.bed

echo "Analysis complete! Check $OUTPUT_DIR/initial_stats.txt for results"
echo "Visualization regions in $OUTPUT_DIR/visualization_regions.bed"