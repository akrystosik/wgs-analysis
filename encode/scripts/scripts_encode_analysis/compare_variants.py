#!/usr/bin/env python3

import pandas as pd
import numpy as np
from pathlib import Path

def load_depmap(depmap_file='data/deepmap/K562_mutations.csv'):
    """Load DepMap variants and standardize format"""
    print("Loading DepMap variants...")
    
    # Read with tab delimiter
    depmap = pd.read_csv(depmap_file, sep='\t')
    
    # Basic stats
    print(f"Loaded {len(depmap)} DepMap variants")
    print("\nVariant types in DepMap:")
    print(depmap['Variant Type'].value_counts())
    
    return depmap

def load_deepsomatic_bed(bed_file='data/variants/tumor_only_analysis/potential_somatic.bed'):
    """Load and process DeepSomatic BED variants"""
    print("\nLoading DeepSomatic BED variants...")
    
    # First peek at the file
    with open(bed_file, 'r') as f:
        first_lines = [next(f) for _ in range(5)]
    print("\nFirst few lines of BED file:")
    for line in first_lines:
        print(line.strip())
    
    # Read BED file
    deep = pd.read_csv(bed_file, sep='\t', header=None,
                      names=['Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand'])
    
    print(f"\nLoaded {len(deep)} DeepSomatic variants")
    print("\nColumn info:")
    for col in deep.columns:
        print(f"{col}: {deep[col].dtype}")
        if col in ['Name', 'Score']:
            print(f"First few unique values: {deep[col].unique()[:5]}")
    
    return deep

def compare_variants(depmap_df, deep_df, window_size=5):
    """Compare variants between datasets"""
    print(f"\nComparing variants (window size: {window_size}bp)...")
    
    # Track matches by chromosome
    matches_by_chr = {}
    window_matches_by_chr = {}
    
    # Process each chromosome separately
    for chrom in depmap_df['Chromosome'].unique():
        depmap_chr = depmap_df[depmap_df['Chromosome'] == chrom]
        deep_chr = deep_df[deep_df['Chromosome'] == chrom]
        
        # Create set of positions for this chromosome
        deep_positions = set(deep_chr['Start'])
        
        # Count matches
        exact_matches = 0
        window_matches = 0
        
        for pos in depmap_chr['Position']:
            # Check exact match
            if pos in deep_positions:
                exact_matches += 1
                continue
                
            # Check window match
            for offset in range(-window_size, window_size + 1):
                if pos + offset in deep_positions:
                    window_matches += 1
                    break
        
        matches_by_chr[chrom] = exact_matches
        window_matches_by_chr[chrom] = window_matches
    
    # Summarize results
    print("\nResults by chromosome:")
    print("Chrom\tDepMap\tDeepSomatic\tExact\tWindow")
    print("-" * 50)
    
    total_exact = 0
    total_window = 0
    
    for chrom in sorted(matches_by_chr.keys()):
        depmap_count = len(depmap_df[depmap_df['Chromosome'] == chrom])
        deep_count = len(deep_df[deep_df['Chromosome'] == chrom])
        exact = matches_by_chr[chrom]
        window = window_matches_by_chr[chrom]
        
        print(f"{chrom}\t{depmap_count}\t{deep_count}\t{exact}\t{window}")
        
        total_exact += exact
        total_window += window
    
    print("-" * 50)
    print(f"Total\t{len(depmap_df)}\t{len(deep_df)}\t{total_exact}\t{total_window}")
    
    return {
        'exact_matches': total_exact,
        'window_matches': total_window,
        'match_rate': total_exact / len(depmap_df),
        'window_match_rate': total_window / len(depmap_df)
    }

def main():
    # Load datasets
    depmap = load_depmap()
    deep = load_deepsomatic_bed()
    
    # Compare variants
    stats = compare_variants(depmap, deep)
    
    # Print summary
    print("\nOverall Summary:")
    print(f"DepMap variants: {len(depmap):,}")
    print(f"DeepSomatic variants: {len(deep):,}")
    print(f"Exact matches: {stats['exact_matches']} ({stats['match_rate']:.2%})")
    print(f"Window matches: {stats['window_matches']} ({stats['window_match_rate']:.2%})")
    print(f"Potential false positive rate: {1 - (stats['exact_matches']/len(deep)):.2%}")

if __name__ == "__main__":
    main()
    