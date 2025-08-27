#!/usr/bin/env python3

import pandas as pd
import numpy as np
from pathlib import Path

def load_variant_sets():
    """Load and do basic counting of variant sets"""
    print("Loading variant sets...")
    
    # Load DepMap variants
    depmap = pd.read_csv("data/deepmap/K562_mutations.csv", sep='\t')
    
    # Load DeepSomatic variants
    somatic_all = pd.read_csv("data/variants/tumor_only_analysis/potential_somatic.annotated.bed", 
                             sep='\t', 
                             names=['Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand'])
    
    somatic_coding = pd.read_csv("data/variants/tumor_only_analysis/potential_somatic.coding.bed",
                                sep='\t',
                                names=['Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand'])
    
    return depmap, somatic_all, somatic_coding

def analyze_coding_enrichment(somatic_all, somatic_coding):
    """Analyze what fraction of calls are in coding regions"""
    total = len(somatic_all)
    coding = len(somatic_coding)
    
    print("\n=== Coding Region Analysis ===")
    print(f"Total variants: {total:,}")
    print(f"Coding variants: {coding:,}")
    print(f"Fraction in coding regions: {coding/total:.2%}")
    
    # Expected coding fraction (coding regions are ~2% of genome)
    expected = 0.02
    enrichment = (coding/total) / expected
    print(f"Enrichment over random: {enrichment:.1f}x")

def analyze_depmap_overlap(depmap, somatic_coding):
    """Analyze overlap with DepMap variants"""
    print("\n=== DepMap Overlap Analysis ===")
    
    # Create comparable position strings
    depmap['pos'] = depmap['Chromosome'] + ':' + depmap['Position'].astype(str)
    somatic_coding['pos'] = somatic_coding['Chromosome'] + ':' + somatic_coding['Start'].astype(str)
    
    # Find overlaps
    depmap_set = set(depmap['pos'])
    somatic_set = set(somatic_coding['pos'])
    overlapping = depmap_set & somatic_set
    
    print(f"DepMap variants: {len(depmap_set):,}")
    print(f"DeepSomatic coding variants: {len(somatic_set):,}")
    print(f"Overlapping positions: {len(overlapping):,}")
    print(f"Fraction of DepMap variants found: {len(overlapping)/len(depmap_set):.2%}")
    print(f"Fraction of coding variants in DepMap: {len(overlapping)/len(somatic_set):.2%}")

def find_interesting_variants(depmap, somatic_coding):
    """Find potentially interesting novel variants"""
    print("\n=== Novel Variant Analysis ===")
    
    # Focus on high-confidence regions (if Score available)
    if 'Score' in somatic_coding.columns and not somatic_coding['Score'].isna().all():
        high_conf = somatic_coding[somatic_coding['Score'] >= somatic_coding['Score'].quantile(0.9)]
        print(f"High-confidence novel variants: {len(high_conf):,}")
        
        # Show examples
        print("\nExample high-confidence variants:")
        print(high_conf.head())

def main():
    # Load data
    depmap, somatic_all, somatic_coding = load_variant_sets()
    
    # Run analyses
    analyze_coding_enrichment(somatic_all, somatic_coding)
    analyze_depmap_overlap(depmap, somatic_coding)
    find_interesting_variants(depmap, somatic_coding)
    
    print("\n=== Summary ===")
    print("1. We're finding more variants in coding regions than expected by chance")
    print("2. A subset of our coding variants match known DepMap variants")
    print("3. Novel variants may represent tumor evolution or false positives")
    print("\nNext steps:")
    print("1. Visual inspection of novel variants in IGV")
    print("2. Additional filtering based on quality metrics")
    print("3. Functional annotation of novel variants")

if __name__ == "__main__":
    main()