#!/bin/bash

# Input files
DEEPSOMATIC="data/variants/deepsomatic/K562_WGS_pair1.deepsomatic.vcf.gz"
DEEPVARIANT="data/variants/deepvariant/K562_WGS_pair1.deepvariant.vcf.gz"

# ENCODE VCFs array
ENCODE_VCFS=(
    "ENCFF752OAX"
    "ENCFF785JVR"
    "ENCFF574MDJ"
)

OUTPUT_DIR="k562_corrected_comparison"
mkdir -p $OUTPUT_DIR

echo "Starting corrected K562 ENCODE validation analysis..."

for encode_vcf in "${ENCODE_VCFS[@]}"; do
    echo "Processing ${encode_vcf}..."
    
    # Create output directory for this comparison
    comp_dir="$OUTPUT_DIR/${encode_vcf}"
    mkdir -p "$comp_dir"
    
    # Compare DeepSomatic with ENCODE
    echo "Comparing DeepSomatic with ${encode_vcf}..."
    bcftools isec \
        -p "$comp_dir/deepsomatic" \
        -n =2 \
        $DEEPSOMATIC \
        "data/variants/encode/${encode_vcf}.vcf.gz"
    
    # Compare DeepVariant with ENCODE
    echo "Comparing DeepVariant with ${encode_vcf}..."
    bcftools isec \
        -p "$comp_dir/deepvariant" \
        -n =2 \
        $DEEPVARIANT \
        "data/variants/encode/${encode_vcf}.vcf.gz"
    
    # Generate statistics
    echo "Generating statistics for ${encode_vcf}..."
    {
        echo "=== ${encode_vcf} Analysis ==="
        echo
        echo "DeepSomatic intersection:"
        wc -l "$comp_dir/deepsomatic/0000.vcf" "$comp_dir/deepsomatic/0001.vcf" 2>/dev/null | grep -v "total"
        echo
        echo "DeepVariant intersection:"
        wc -l "$comp_dir/deepvariant/0000.vcf" "$comp_dir/deepvariant/0001.vcf" 2>/dev/null | grep -v "total"
        echo
        # Compare actual variants between callers
        echo "Comparing variants between callers..."
        comm -12 \
            <(bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$comp_dir/deepsomatic/0002.vcf" 2>/dev/null | sort) \
            <(bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$comp_dir/deepvariant/0002.vcf" 2>/dev/null | sort) \
            > "$comp_dir/shared_variants.txt"
        echo "Variants found by both callers: $(wc -l < "$comp_dir/shared_variants.txt")"
    } > "$comp_dir/comparison_stats.txt"
    
    # Generate VAF distribution for DeepSomatic concordant variants
    echo "Analyzing VAF distribution for ${encode_vcf}..."
    bcftools query -f '[%VAF]\n' "$comp_dir/deepsomatic/0002.vcf" 2>/dev/null | \
    awk 'BEGIN {bins=10}
         {
             if($1!="." && $1>=0 && $1<=1) {
                 bin = int($1*bins);
                 count[bin]++;
                 total++;
             }
         }
         END {
             print "VAF Distribution for DeepSomatic concordant variants:";
             for (i=0; i<bins; i++)
                 printf("%.1f-%.1f: %d (%.1f%%)\n", 
                        i/bins, (i+1)/bins, 
                        count[i], 
                        count[i]/total*100);
         }' > "$comp_dir/vaf_distribution.txt"
done

# Generate summary report
{
    echo "K562 Variant Calling Analysis Summary"
    echo "===================================="
    echo
    echo "Analysis Date: $(date)"
    echo
    for encode_vcf in "${ENCODE_VCFS[@]}"; do
        echo "=== ${encode_vcf} Results ==="
        cat "$OUTPUT_DIR/${encode_vcf}/comparison_stats.txt"
        echo
        echo "VAF Distribution:"
        cat "$OUTPUT_DIR/${encode_vcf}/vaf_distribution.txt"
        echo
        echo "----------------------------------------"
    done
} > "$OUTPUT_DIR/analysis_summary.txt"

echo "Analysis complete! Check $OUTPUT_DIR/analysis_summary.txt for results"