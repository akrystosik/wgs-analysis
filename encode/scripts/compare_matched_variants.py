#!/usr/bin/env python3

import pandas as pd

# Adjust these paths to your actual files
matched_file = "data/variants/tumor_only_analysis/somatic_matches.bed"   # The 18 "matched" variants
depmap_bed   = "data/variants/tumor_only_analysis/depmap_variants.bed"   # ~289 lines from DepMap
roi_bed      = "data/variants/tumor_only_analysis/igv_report/regions_of_interest.bed"

bed_cols = ["Chrom","Start","End","Name","Score","Strand"]
matched_df = pd.read_csv(matched_file, sep='\t', names=bed_cols)
depmap_df  = pd.read_csv(depmap_bed,   sep='\t', names=bed_cols)
roi_df     = pd.read_csv(roi_bed,      sep='\t', names=bed_cols)

print(f"=== Basic Counts ===")
print(f"  matched.bed:  {len(matched_df)} variants")
print(f"  depmap.bed:   {len(depmap_df)} variants")
print(f"  ROI.bed:      {len(roi_df)} variants")
print()

def variant_key_df(df):
    return df[["Chrom","Start","End"]].drop_duplicates()

matched_keys = variant_key_df(matched_df)
depmap_keys  = variant_key_df(depmap_df)
roi_keys     = variant_key_df(roi_df)

match_depmap = pd.merge(matched_keys, depmap_keys, on=["Chrom","Start","End"], how="inner")
match_roi    = pd.merge(matched_keys, roi_keys,    on=["Chrom","Start","End"], how="inner")
depmap_roi   = pd.merge(depmap_keys,  roi_keys,    on=["Chrom","Start","End"], how="inner")

print("=== Overlap Summary ===")
print(f"Overlaps between matched & DepMap:  {len(match_depmap)}")
print(f"Overlaps between matched & ROI:     {len(match_roi)}")
print(f"Overlaps between DepMap & ROI:      {len(depmap_roi)}")
print()

# Identify which matched variants are missing from DepMap or ROI
matched_not_in_depmap = pd.merge(matched_keys, depmap_keys, on=["Chrom","Start","End"], how="left", indicator=True)
matched_not_in_roi    = pd.merge(matched_keys, roi_keys,    on=["Chrom","Start","End"], how="left", indicator=True)

num_missing_depmap = (matched_not_in_depmap["_merge"]=="left_only").sum()
num_missing_roi    = (matched_not_in_roi["_merge"]=="left_only").sum()

print("=== 'Matched' variants that are missing in DepMap/ROI ===")
print(f"Missing from DepMap: {num_missing_depmap}")
print(f"Missing from ROI:    {num_missing_roi}")
print()

# Truncate output so we only see the first 5 missing lines
MAX_SHOW = 5

if num_missing_depmap > 0:
    missing_rows = matched_not_in_depmap[matched_not_in_depmap["_merge"]=="left_only"].copy()
    # Show at most 5 rows
    to_show = missing_rows.head(MAX_SHOW)
    print("Matched variants missing in DepMap (showing up to 5):")
    print(to_show[["Chrom","Start","End"]].to_string(index=False))
    if len(missing_rows) > MAX_SHOW:
        print(f"... {len(missing_rows) - MAX_SHOW} more omitted ...")
    print()

if num_missing_roi > 0:
    missing_rows = matched_not_in_roi[matched_not_in_roi["_merge"]=="left_only"].copy()
    # Show at most 5 rows
    to_show = missing_rows.head(MAX_SHOW)
    print("Matched variants missing in ROI (showing up to 5):")
    print(to_show[["Chrom","Start","End"]].to_string(index=False))
    if len(missing_rows) > MAX_SHOW:
        print(f"... {len(missing_rows) - MAX_SHOW} more omitted ...")
    print()

print("Done. Summary only—no massive output!")
