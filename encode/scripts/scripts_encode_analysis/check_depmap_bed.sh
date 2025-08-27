#!/bin/bash

DEP_MAP_BED="data/variants/tumor_only_analysis/depmap_variants.bed"

echo "=== DepMap BED Info ==="
echo "DepMap BED: $DEP_MAP_BED"
echo

# Show first few lines
echo "[HEAD of DepMap BED]"
head -n 20 "$DEP_MAP_BED"
echo

# Look for duplicates by full coordinate (chr+start+end+name)
echo "[Check for duplicate lines in DepMap BED]"
sort -k1,1 -k2,2n -k3,3n -k4,4 "$DEP_MAP_BED" | uniq -c | awk '$1>1 {print "DUPLICATE:", $0}'
echo

# Look for duplicates by 'Name' alone, e.g. "DepMap_9"
echo "[Check for repeated 'Name' fields]"
cut -f4 "$DEP_MAP_BED" | sort | uniq -c | sort -nr | head -n 10
