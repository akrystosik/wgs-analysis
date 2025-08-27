#!/bin/bash

ROIBED="data/variants/tumor_only_analysis/igv_report/regions_of_interest.bed"

echo "=== IGV Regions of Interest Info ==="
echo "regions_of_interest.bed: $ROIBED"
echo

# Show lines specifically labeled DepMap
echo "[Lines labeled 'DEPMAP' in ROI BED]"
grep -n "DEPMAP_" "$ROIBED"
echo

# Check for duplicates by chromosome+start+end+name
echo "[Check for duplicates among DepMap lines]"
grep "DEPMAP_" "$ROIBED" | sort -k1,1 -k2,2n -k3,3n -k4,4 | uniq -c | awk '$1>1 {print "DUPLICATE:", $0}'
echo

# Quick summary: how many total distinct DepMap lines?
echo "[Count distinct DepMap lines]"
grep "DEPMAP_" "$ROIBED" | sort -u | wc -l
