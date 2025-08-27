#!/bin/bash

DEP_MAP_IN="data/deepmap/K562_mutations.csv"
DEP_MAP_BED="data/variants/tumor_only_analysis/depmap_variants.bed"

awk -F '\t' 'NR>1 {
    chrom = $2;
    pos   = $3;
    start = pos - 1;
    end   = pos;
    print chrom"\t"start"\t"end"\tDepMap_"NR"\t0\t+"
}' "$DEP_MAP_IN" > "$DEP_MAP_BED"

echo "Created DepMap BED with $(wc -l < "$DEP_MAP_BED") variants."
