#!/bin/bash

OUTPUT_DIR="k562_debug_comparison"
mkdir -p $OUTPUT_DIR

# 1. First, let's check the actual content of our normalized files
echo "Checking first 5 lines of each normalized file..."
echo "ENCODE variants:" > $OUTPUT_DIR/debug_output.txt
head -n 5 k562_improved_comparison/encode_normalized.txt >> $OUTPUT_DIR/debug_output.txt
echo -e "\nDeepSomatic variants:" >> $OUTPUT_DIR/debug_output.txt
head -n 5 k562_improved_comparison/deepsomatic_normalized.txt >> $OUTPUT_DIR/debug_output.txt

# 2. Count variants per caller
echo -e "\nVariant counts:" >> $OUTPUT_DIR/debug_output.txt
echo "ENCODE variants: $(wc -l k562_improved_comparison/encode_normalized.txt)" >> $OUTPUT_DIR/debug_output.txt
echo "DeepSomatic variants: $(wc -l k562_improved_comparison/deepsomatic_normalized.txt)" >> $OUTPUT_DIR/debug_output.txt

# 3. Try a simpler intersection approach
echo -e "\nDirect intersection test:" >> $OUTPUT_DIR/debug_output.txt
# Create position-only files for each
cut -f1,2 k562_improved_comparison/encode_normalized.txt | sort > $OUTPUT_DIR/encode_pos.txt
cut -f1,2 k562_improved_comparison/deepsomatic_normalized.txt | sort > $OUTPUT_DIR/deepsomatic_pos.txt

# Find exact position matches
comm -12 $OUTPUT_DIR/encode_pos.txt $OUTPUT_DIR/deepsomatic_pos.txt > $OUTPUT_DIR/common_positions.txt
echo "Positions in both files: $(wc -l < $OUTPUT_DIR/common_positions.txt)" >> $OUTPUT_DIR/debug_output.txt

# 4. Check a specific region in detail
echo -e "\nDetailed check of chr1:10000-20000:" >> $OUTPUT_DIR/debug_output.txt
echo "ENCODE variants in region:" >> $OUTPUT_DIR/debug_output.txt
awk '$1=="chr1" && $2>=10000 && $2<=20000' k562_improved_comparison/encode_normalized.txt >> $OUTPUT_DIR/debug_output.txt
echo -e "\nDeepSomatic variants in region:" >> $OUTPUT_DIR/debug_output.txt
awk '$1=="chr1" && $2>=10000 && $2<=20000' k562_improved_comparison/deepsomatic_normalized.txt >> $OUTPUT_DIR/debug_output.txt

# 5. Try a stricter comparison for one chromosome
echo -e "\nStrict comparison for chr1:" >> $OUTPUT_DIR/debug_output.txt
CHROM="chr1"
# Extract chr1 variants
awk -v chr=$CHROM '$1==chr' k562_improved_comparison/encode_normalized.txt > $OUTPUT_DIR/encode_chr1.txt
awk -v chr=$CHROM '$1==chr' k562_improved_comparison/deepsomatic_normalized.txt > $OUTPUT_DIR/deepsomatic_chr1.txt

# Compare exact matches (position and alleles)
awk 'NR==FNR{a[$1,$2,$3,$4]++;next} ($1,$2,$3,$4) in a' \
    $OUTPUT_DIR/encode_chr1.txt $OUTPUT_DIR/deepsomatic_chr1.txt > $OUTPUT_DIR/exact_matches_chr1.txt

echo "Chr1 exact matches (pos + alleles): $(wc -l < $OUTPUT_DIR/exact_matches_chr1.txt)" >> $OUTPUT_DIR/debug_output.txt

echo "Debug analysis complete! Check $OUTPUT_DIR/debug_output.txt for results"