#!/bin/bash

# Cell lines to compare
CELL_LINES=(
    "K562:data/variants/deepsomatic/K562_WGS_pair1.deepsomatic.vcf.gz:data/variants/deepvariant/K562_WGS_pair1.deepvariant.vcf.gz"
    "Caki2:data/variants/deepsomatic/Caki2_somatic_full.vcf:data/variants/deepvariant/Caki2_WGS_pair1.vcf.gz"
)

# ENCODE VCFs for K562
ENCODE_VCFS=(
    "ENCFF752OAX"
    "ENCFF785JVR"
    "ENCFF574MDJ"
    "ENCFF863MPP"
)

OUTPUT_DIR="cell_line_comparison"
mkdir -p $OUTPUT_DIR

echo "Starting cell line comparison analysis..."

for cell_line_info in "${CELL_LINES[@]}"; do
    IFS=':' read -r cell_line deepsomatic deepvariant <<< "$cell_line_info"
    echo "Processing $cell_line..."
    
    mkdir -p "$OUTPUT_DIR/$cell_line"
    
    # Basic caller comparison
    echo "Comparing DeepSomatic and DeepVariant for $cell_line..."
    bcftools isec \
        -p "$OUTPUT_DIR/$cell_line/caller_comparison" \
        $deepsomatic \
        $deepvariant
    
    # Get VAF distribution for DeepSomatic
    echo "Analyzing VAF distribution..."
    bcftools query -f '[%VAF]\n' $deepsomatic | \
    awk '{
        if($1>=0 && $1<=1) {
            bin=int($1*10);
            count[bin]++;
            total++;
        }
    }
    END {
        for(i=0;i<10;i++)
            printf("%.1f-%.1f: %d (%.1f%%)\n",
                   i/10,(i+1)/10,count[i],count[i]/total*100)
    }' > "$OUTPUT_DIR/$cell_line/vaf_distribution.txt"

    # If K562, compare with ENCODE
    if [ "$cell_line" == "K562" ]; then
        for encode_vcf in "${ENCODE_VCFS[@]}"; do
            echo "Comparing with $encode_vcf..."
            bcftools isec \
                -p "$OUTPUT_DIR/$cell_line/encode_${encode_vcf}" \
                $deepsomatic \
                data/variants/encode_lifted/${encode_vcf}.hg38.sorted.vcf.gz
        done
    fi
done

# Generate comparison report
{
    echo "Cell Line Comparison Analysis"
    echo "==========================="
    echo
    
    for cell_line_info in "${CELL_LINES[@]}"; do
        IFS=':' read -r cell_line deepsomatic deepvariant <<< "$cell_line_info"
        echo "$cell_line Analysis:"
        echo "-----------------"
        
        # Get caller comparison stats
        DS_ONLY=$(grep -vc "^#" "$OUTPUT_DIR/$cell_line/caller_comparison/0000.vcf")
        DV_ONLY=$(grep -vc "^#" "$OUTPUT_DIR/$cell_line/caller_comparison/0001.vcf")
        SHARED=$(grep -vc "^#" "$OUTPUT_DIR/$cell_line/caller_comparison/0002.vcf")
        
        echo "Caller Comparison:"
        echo "  Unique to DeepSomatic: $DS_ONLY"
        echo "  Unique to DeepVariant: $DV_ONLY"
        echo "  Shared variants: $SHARED"
        echo
        
        echo "VAF Distribution:"
        cat "$OUTPUT_DIR/$cell_line/vaf_distribution.txt"
        echo
        
        if [ "$cell_line" == "K562" ]; then
            echo "ENCODE Comparisons:"
            for encode_vcf in "${ENCODE_VCFS[@]}"; do
                echo "  $encode_vcf:"
                echo "    Shared variants: $(grep -vc "^#" "$OUTPUT_DIR/$cell_line/encode_${encode_vcf}/0002.vcf")"
            done
        fi
        echo
    done
} > "$OUTPUT_DIR/comparison_report.txt"

echo "Analysis complete! Check $OUTPUT_DIR/comparison_report.txt for results"