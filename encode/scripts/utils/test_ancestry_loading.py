#!/usr/bin/env python3
"""
Test script to verify MAGE ancestry data loading
"""

import pandas as pd
import numpy as np

def test_load_mage_ancestry_data():
    """Test loading MAGE ancestry assignment data"""
    print("Testing MAGE ancestry data loading...")
    
    # Load the 1000 Genomes population data
    pop_file = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/mage/20130606_g1k_3202_samples_ped_population.txt"
    pop_df = pd.read_csv(pop_file, sep=' ')
    print(f"Population file columns: {list(pop_df.columns)}")
    print(f"Population file shape: {pop_df.shape}")
    print("First few rows:")
    print(pop_df.head())
    
    # Load MAGE metadata
    metadata_file = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/mage/sample.metadata.MAGE.v1.0.txt"
    metadata_df = pd.read_csv(metadata_file, sep='\t')
    print(f"\nMetadata file columns: {list(metadata_df.columns)}")
    print(f"Metadata file shape: {metadata_df.shape}")
    print("First few rows:")
    print(metadata_df.head())
    
    # Create sample to ancestry mapping
    sample_to_ancestry = {}
    sample_to_population = {}
    
    # From population file (HG/NA prefixed samples)
    for _, row in pop_df.iterrows():
        sample_id = row['SampleID']
        sample_to_ancestry[sample_id] = row['Superpopulation']
        sample_to_population[sample_id] = row['Population']
    
    # From metadata file (GM/NA prefixed samples)
    for _, row in metadata_df.iterrows():
        sample_id = row['sample_kgpID']  # NA12843 format
        sample_to_ancestry[sample_id] = row['continentalGroup']
        sample_to_population[sample_id] = row['population']
    
    print(f"\nLoaded ancestry data for {len(sample_to_ancestry)} samples")
    
    # Print ancestry distribution
    ancestry_counts = pd.Series(list(sample_to_ancestry.values())).value_counts()
    print("Ancestry distribution:")
    for ancestry, count in ancestry_counts.items():
        print(f"  {ancestry}: {count}")
    
    # Show some examples
    print("\nSample examples:")
    sample_examples = list(sample_to_ancestry.items())[:10]
    for sample_id, ancestry in sample_examples:
        population = sample_to_population.get(sample_id, 'Unknown')
        print(f"  {sample_id}: {ancestry} ({population})")
    
    return sample_to_ancestry, sample_to_population

if __name__ == "__main__":
    sample_to_ancestry, sample_to_population = test_load_mage_ancestry_data()
    print(f"\n✅ Successfully loaded ancestry data for {len(sample_to_ancestry)} samples")