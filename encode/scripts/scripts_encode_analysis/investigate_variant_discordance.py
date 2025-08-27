#!/usr/bin/env python3

import pandas as pd
import numpy as np
from pathlib import Path

def load_datasets():
    """Load and standardize variant datasets"""
    # Load DepMap variants
    depmap = pd.read_csv("data/deepmap/K562_mutations.csv", sep='\t')
    
    # Load our coding variants
    coding = pd.read_csv("data/variants/tumor_only_analysis/potential_somatic.coding.bed",
                        sep='\t',
                        names=['Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand'])
    
    print(f"Loaded {len(depmap)} DepMap variants and {len(coding)} coding variants")
    return depmap, coding

def analyze_missed_depmap_variants(depmap, coding):
    """Analyze DepMap variants we're missing"""
    # Create position keys
    depmap['pos_key'] = depmap['Chromosome'] + ':' + depmap['Position'].astype(str)
    coding['pos_key'] = coding['Chromosome'] + ':' + coding['Start'].astype(str)
    
    # Find missed variants
    found_variants = set(coding['pos_key'])
    missed_variants = depmap[~depmap['pos_key'].isin(found_variants)]
    
    print("\n=== Analysis of Missed DepMap Variants ===")
    print(f"Total missed: {len(missed_variants)} / {len(depmap)} ({len(missed_variants)/len(depmap):.1%})")
    
    # Analyze by variant type
    print("\nMissed variants by type:")
    print(missed_variants['Variant Type'].value_counts())
    
    # Analyze by allele fraction
    if 'Allele Fraction' in missed_variants.columns:
        print("\nAllele fraction of missed variants:")
        print(missed_variants['Allele Fraction'].describe())
    
    return missed_variants

def analyze_novel_coding_variants(depmap, coding):
    """Analyze our novel coding variants"""
    # Create position keys
    depmap['pos_key'] = depmap['Chromosome'] + ':' + depmap['Position'].astype(str)
    coding['pos_key'] = coding['Chromosome'] + ':' + coding['Start'].astype(str)
    
    # Find novel variants
    depmap_positions = set(depmap['pos_key'])
    novel_variants = coding[~coding['pos_key'].isin(depmap_positions)]
    
    print("\n=== Analysis of Novel Coding Variants ===")
    print(f"Total novel: {len(novel_variants)} / {len(coding)} ({len(novel_variants)/len(coding):.1%})")
    
    # Distribution by chromosome
    chrom_dist = novel_variants['Chromosome'].value_counts()
    print("\nChromosomal distribution of novel variants:")
    for chrom, count in chrom_dist.items():
        print(f"{chrom}: {count:,} variants")
    
    return novel_variants

def suggest_investigation_strategy(missed_variants, novel_variants):
    """Suggest strategy to investigate discordant variants"""
    print("\n=== Investigation Strategy ===")
    
    # For missed DepMap variants
    print("\n1. Priority DepMap variants to investigate:")
    high_af_missed = missed_variants[missed_variants['Allele Fraction'] > 0.3]
    print(f"- {len(high_af_missed)} high allele fraction variants")
    
    # For novel variants
    print("\n2. Novel variants to investigate:")
    if 'Score' in novel_variants.columns:
        high_conf = novel_variants[novel_variants['Score'] >= novel_variants['Score'].quantile(0.9)]
        print(f"- Top {len(high_conf)} high-confidence novel variants")
    
    # Output investigation BED files
    output_dir = Path("analysis_output")
    output_dir.mkdir(exist_ok=True)
    
    # Save high-AF missed variants
    if len(high_af_missed) > 0:
        out_file = output_dir / "missed_highaf_variants.bed"
        high_af_missed.to_csv(out_file, sep='\t', index=False)
        print(f"\nSaved high-AF missed variants to: {out_file}")
    
    # Output commands for IGV inspection
    print("\nNext steps:")
    print("1. Check genome build compatibility")
    print("2. Inspect high-AF missed variants in IGV")
    print("3. Review variant calling parameters")
    print("4. Consider adding matched normal or panel of normals")

def main():
    # Load data
    depmap, coding = load_datasets()
    
    # Analyze discordant variants
    missed_variants = analyze_missed_depmap_variants(depmap, coding)
    novel_variants = analyze_novel_coding_variants(depmap, coding)
    
    # Suggest investigation strategy
    suggest_investigation_strategy(missed_variants, novel_variants)

if __name__ == "__main__":
    main()