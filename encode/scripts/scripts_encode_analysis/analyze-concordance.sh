#!/bin/bash

OUTPUT_DIR="k562_normalized_comparison"

echo "Analyzing normalized variants..."

# 1. Check normalized VCF stats
for vcf in $OUTPUT_DIR/*_normalized.vcf; do
    echo "Stats for $(basename $vcf):"
    bcftools stats $vcf | grep "^SN"
done > $OUTPUT_DIR/normalized_stats.txt

# 2. Extract variants in BED format for overlapping analysis
for vcf in $OUTPUT_DIR/*_normalized.vcf; do
    base=$(basename $vcf .vcf)
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\n' $vcf > \
        $OUTPUT_DIR/${base}_variants.txt
done

# 3. Find overlapping variants
python3 - <<EOF
import pandas as pd

# Read variant files
encode = pd.read_csv('$OUTPUT_DIR/ENCFF752OAX_normalized_variants.txt', sep='\t', 
                    names=['CHROM', 'POS', 'REF', 'ALT', 'QUAL'])
ds = pd.read_csv('$OUTPUT_DIR/K562_WGS_pair1.deepsomatic_normalized_variants.txt', 
                 sep='\t', names=['CHROM', 'POS', 'REF', 'ALT', 'QUAL'])
dv = pd.read_csv('$OUTPUT_DIR/K562_WGS_pair1.deepvariant_normalized_variants.txt', 
                 sep='\t', names=['CHROM', 'POS', 'REF', 'ALT', 'QUAL'])

# Find overlapping positions
encode_pos = set(zip(encode.CHROM, encode.POS))
ds_pos = set(zip(ds.CHROM, ds.POS))
dv_pos = set(zip(dv.CHROM, dv.POS))

# Calculate overlaps
all_overlap = encode_pos & ds_pos & dv_pos
ds_dv_overlap = ds_pos & dv_pos

# Write results
with open('$OUTPUT_DIR/overlap_analysis.txt', 'w') as f:
    f.write(f'Total variants:\n')
    f.write(f'ENCODE: {len(encode):,}\n')
    f.write(f'DeepSomatic: {len(ds):,}\n')
    f.write(f'DeepVariant: {len(dv):,}\n\n')
    f.write(f'Overlaps:\n')
    f.write(f'All callers: {len(all_overlap):,}\n')
    f.write(f'DeepSomatic-DeepVariant: {len(ds_dv_overlap):,}\n')

# Create BED file for visualization
with open('$OUTPUT_DIR/overlap_regions.bed', 'w') as f:
    for chrom, pos in all_overlap:
        f.write(f'{chrom}\t{pos-100}\t{pos+100}\tAll_callers\n')
    for chrom, pos in ds_dv_overlap - all_overlap:
        f.write(f'{chrom}\t{pos-100}\t{pos+100}\tDS_DV_only\n')
EOF

echo "Analysis complete! Check $OUTPUT_DIR/overlap_analysis.txt for results"