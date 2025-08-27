#!/usr/bin/env python3

import pandas as pd
import numpy as np
from pathlib import Path

def analyze_highaf_missed():
    """Analyze high-AF variants we're missing from DepMap"""
    
    # Load high-AF missed variants
    missed = pd.read_csv("analysis_output/missed_highaf_variants.bed", sep='\t')
    print(f"Analyzing {len(missed)} high-AF missed variants")
    
    # Group by chromosomes and look for patterns
    chrom_counts = missed['Chromosome'].value_counts().sort_index()
    print("\nChromosomal distribution of missed variants:")
    for chrom, count in chrom_counts.items():
        print(f"{chrom}: {count} variants")
        
    # Analyze variant properties
    print("\nVariant Properties:")
    print("Type distribution:")
    print(missed['Variant Type'].value_counts())
    
    print("\nAllele Fraction distribution:")
    af_bins = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    af_dist = pd.cut(missed['Allele Fraction'], bins=af_bins)
    print(af_dist.value_counts().sort_index())
    
    # Create BED file for IGV visualization
    create_igv_batch(missed)
    
    # Look for potential systematic issues
    identify_systematic_issues(missed)

def create_igv_batch(missed):
    """Create IGV batch script to visualize missed variants"""
    output_dir = Path("analysis_output")
    
    # Create IGV batch script
    with open(output_dir / "inspect_missed_variants.igv.bat", 'w') as f:
        f.write("new\ngenome hg38\n")
        f.write("snapshotDirectory analysis_output/igv_snapshots\n")
        
        # Sort variants by allele fraction
        missed_sorted = missed.sort_values('Allele Fraction', ascending=False)
        
        # Create goto commands for top 50 highest AF variants
        for idx, row in missed_sorted.head(50).iterrows():
            pos = int(row['Position'])
            f.write(f"goto {row['Chromosome']}:{pos-50}-{pos+50}\n")
            f.write(f"snapshot {row['Chromosome']}_{pos}.png\n")
    
    print("\nCreated IGV batch script for visualization")

def identify_systematic_issues(missed):
    """Look for patterns that might indicate systematic issues"""
    print("\nInvestigating potential systematic issues:")
    
    # Check for clustering
    missed['binned_pos'] = missed['Position'] // 1000000  # 1Mb bins
    clusters = missed.groupby(['Chromosome', 'binned_pos']).size()
    clusters = clusters[clusters > 1].sort_values(ascending=False)
    
    if len(clusters) > 0:
        print("\nFound variant clusters (>1 variant per Mb):")
        for (chrom, bin_pos), count in clusters.items():
            start = bin_pos * 1000000
            end = (bin_pos + 1) * 1000000
            print(f"{chrom}:{start}-{end}: {count} variants")
    
    # Look for base change patterns
    if all(col in missed.columns for col in ['Ref Allele', 'Alt Allele']):
        print("\nBase change patterns:")
        missed['change'] = missed['Ref Allele'] + '>' + missed['Alt Allele']
        print(missed['change'].value_counts())

def suggest_next_steps(clusters_found):
    """Suggest next steps based on analysis"""
    print("\nRecommended next steps:")
    
    print("1. Variant Caller Parameters:")
    print("   - Review minimum coverage thresholds")
    print("   - Check allele fraction cutoffs")
    print("   - Verify match with DepMap's genome build")
    
    if clusters_found:
        print("2. Investigate clustered regions:")
        print("   - Check for systematic alignment issues")
        print("   - Look for problematic sequence contexts")
    
    print("3. Validation Steps:")
    print("   - IGV inspection of highest-AF missed variants")
    print("   - Compare raw alignment files if available")
    print("   - Consider running alternative variant caller")

def main():
    print("=== Analysis of High-AF Missed DepMap Variants ===")
    
    # Run analysis
    analyze_highaf_missed()
    
    # Output files created
    print("\nOutput files:")
    print("1. analysis_output/inspect_missed_variants.igv.bat")
    print("   - IGV batch script for visualizing top 50 missed variants")
    print("2. analysis_output/igv_snapshots/")
    print("   - Directory for IGV snapshots")

if __name__ == "__main__":
    main()