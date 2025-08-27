#!/bin/bash

# Input files
DEEPSOMATIC="data/variants/deepsomatic/K562_WGS_pair1.deepsomatic.vcf.gz"
DEEPVARIANT="data/variants/deepvariant/K562_WGS_pair1.deepvariant.vcf.gz"
REGIONS_BED="k562_regions/focused/top_regions.bed"
OUTPUT_DIR="k562_comparison"

mkdir -p $OUTPUT_DIR

echo "Analyzing variants in selected regions..."

# 1. Extract DeepSomatic variants in regions
bcftools view -R $REGIONS_BED $DEEPSOMATIC | \
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%VAF]\t%FILTER\n' \
    > $OUTPUT_DIR/deepsomatic_variants.txt

# 2. Extract DeepVariant variants in same regions
bcftools view -R $REGIONS_BED $DEEPVARIANT | \
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\t%QUAL\n' \
    > $OUTPUT_DIR/deepvariant_variants.txt

# 3. Create comparison report
echo "Creating IGV report with both variant sets..."
create_report \
  $REGIONS_BED \
  /mnt/czi-sci-ai/intrinsic-variation-gene-ex/encode_analysis/data/reference/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta \
  --flanking 1000 \
  --output k562_igv_report/k562_comparison_report.html \
  --tracks \
    /mnt/czi-sci-ai/intrinsic-variation-gene-ex/encode_analysis/data/bam/K562/K562_WGS_pair1.marked_duplicates.bam \
    $DEEPSOMATIC \
    $DEEPVARIANT

# 4. Generate summary statistics
echo "Generating comparison statistics..."
echo "DeepSomatic variants in regions:" $(wc -l < $OUTPUT_DIR/deepsomatic_variants.txt)
echo "DeepVariant variants in regions:" $(wc -l < $OUTPUT_DIR/deepvariant_variants.txt)

# 5. Analyze concordance
echo "Analyzing variant concordance..."
paste $OUTPUT_DIR/deepsomatic_variants.txt $OUTPUT_DIR/deepvariant_variants.txt | \
    awk -F'\t' '{
        if($1==$7 && $2==$8 && $3==$9 && $4==$10) 
            print "CONCORDANT\t" $0
        else 
            print "DISCORDANT\t" $0
    }' > $OUTPUT_DIR/variant_comparison.txt

echo "Analysis complete! Check k562_igv_report/k562_comparison_report.html for visualization"