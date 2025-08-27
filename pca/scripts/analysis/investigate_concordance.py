#!/usr/bin/env python3
"""
Investigate Concordance Analysis
Detailed examination of the concordance between self-reported ethnicity and genomic ancestry
"""

import pandas as pd
import numpy as np
from pathlib import Path

# Set up file paths
base_dir = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/wgs")
pca_dir = base_dir / "pca"

def investigate_concordance():
    """Investigate concordance analysis in detail"""
    print("Investigating concordance analysis...")
    
    # Load complete results
    complete_file = pca_dir / "results" / "combined_analysis" / "complete_ancestry_results.csv"
    data = pd.read_csv(complete_file)
    
    print(f"Total samples: {len(data)}")
    
    # Check datasets with self-reported ethnicity
    datasets_with_ethnicity = data[data['self_reported_ethnicity'].notna()]
    print(f"Samples with self-reported ethnicity: {len(datasets_with_ethnicity)}")
    
    for dataset in data['dataset'].unique():
        if pd.isna(dataset):
            continue
            
        subset = data[data['dataset'] == dataset]
        with_ethnicity = subset[subset['self_reported_ethnicity'].notna()]
        
        print(f"\n{dataset} dataset:")
        print(f"  Total samples: {len(subset)}")
        print(f"  With self-reported ethnicity: {len(with_ethnicity)}")
        
        if len(with_ethnicity) > 0:
            print("  Self-reported ethnicity distribution:")
            ethnicity_counts = with_ethnicity['self_reported_ethnicity'].value_counts()
            for ethnicity, count in ethnicity_counts.head(10).items():
                print(f"    {ethnicity}: {count}")
    
    # Focus on concordance analysis
    print("\n" + "="*60)
    print("CONCORDANCE ANALYSIS INVESTIGATION")
    print("="*60)
    
    # Recreate the analysis logic from combined_ancestry_analysis.py
    analysis_data = data[
        (data['self_reported_ethnicity'].notna()) & 
        (data['final_consensus_ancestry'] != 'Unknown') &
        (data['final_consensus_ancestry'] != 'Uncertain')
    ].copy()
    
    print(f"Samples for concordance analysis: {len(analysis_data)}")
    
    # Standardize ethnicity labels (same mapping as original)
    ethnicity_mapping = {
        'white': 'EUR',
        'european': 'EUR',
        'caucasian': 'EUR',
        'black or african american': 'AFR', 
        'african american': 'AFR',
        'african': 'AFR',
        'asian': 'EAS',
        'east asian': 'EAS',
        'chinese': 'EAS',
        'japanese': 'EAS',
        'korean': 'EAS',
        'south asian': 'SAS',
        'indian': 'SAS',
        'hispanic': 'AMR',
        'latino': 'AMR',
        'mixed': 'MIXED',
        'other': 'OTHER',
        'multiethnic': 'MIXED',
        'american indian or alaska native': 'AMR',
        'native hawaiian or other pacific islander': 'OTHER',
        'unknown or not reported': 'OTHER'
    }
    
    analysis_data['self_reported_ancestry'] = analysis_data['self_reported_ethnicity'].str.lower().map(ethnicity_mapping)
    
    print(f"Samples with mapped ancestry: {analysis_data['self_reported_ancestry'].notna().sum()}")
    
    # Show unmapped ethnicity labels
    unmapped = analysis_data[analysis_data['self_reported_ancestry'].isna()]
    if len(unmapped) > 0:
        print("\nUnmapped ethnicity labels:")
        unmapped_counts = unmapped['self_reported_ethnicity'].value_counts()
        for label, count in unmapped_counts.items():
            print(f"  {label}: {count}")
    
    # Calculate concordance by dataset
    print("\n" + "-"*50)
    print("DATASET-SPECIFIC CONCORDANCE ANALYSIS")
    print("-"*50)
    
    for dataset in analysis_data['dataset'].unique():
        if pd.isna(dataset):
            continue
            
        subset = analysis_data[
            (analysis_data['dataset'] == dataset) & 
            (analysis_data['self_reported_ancestry'].notna())
        ]
        
        if len(subset) == 0:
            continue
            
        print(f"\n{dataset.upper()} DATASET:")
        print(f"Samples for analysis: {len(subset)}")
        
        # Show breakdown of comparisons
        print("\nSelf-reported vs Genomic Ancestry Comparison:")
        comparison_counts = subset.groupby(['self_reported_ancestry', 'final_consensus_ancestry']).size().reset_index(name='count')
        
        for _, row in comparison_counts.iterrows():
            concordant = "✓" if row['self_reported_ancestry'] == row['final_consensus_ancestry'] else "✗"
            print(f"  {row['self_reported_ancestry']} → {row['final_consensus_ancestry']}: {row['count']} samples {concordant}")
        
        # Calculate concordance
        concordant = subset['self_reported_ancestry'] == subset['final_consensus_ancestry']
        concordance_rate = concordant.sum() / len(subset) * 100
        
        print(f"\nConcordance: {concordant.sum()}/{len(subset)} = {concordance_rate:.1f}%")
        
        # Show discordant cases in detail
        discordant = subset[~concordant]
        if len(discordant) > 0:
            print(f"\nDiscordant cases ({len(discordant)} samples):")
            discord_summary = discordant.groupby(['self_reported_ancestry', 'final_consensus_ancestry']).size().reset_index(name='count')
            for _, row in discord_summary.iterrows():
                print(f"  {row['self_reported_ancestry']} → {row['final_consensus_ancestry']}: {row['count']}")

if __name__ == "__main__":
    investigate_concordance()