#!/usr/bin/env python3
"""
Enhanced ancestry loader that combines MAGE data with RNA-seq ethnicity metadata
"""

import pandas as pd
import numpy as np
import os

def load_mage_ancestry_data():
    """Load original MAGE ancestry data"""
    sample_to_ancestry = {}
    sample_to_population = {}
    
    # Load the 1000 Genomes population data
    pop_file = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/mage/20130606_g1k_3202_samples_ped_population.txt"
    if os.path.exists(pop_file):
        pop_df = pd.read_csv(pop_file, sep=' ')
        for _, row in pop_df.iterrows():
            sample_id = row['SampleID']
            sample_to_ancestry[sample_id] = row['Superpopulation']
            sample_to_population[sample_id] = row['Population']
    
    # Load MAGE metadata
    metadata_file = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/mage/sample.metadata.MAGE.v1.0.txt"
    if os.path.exists(metadata_file):
        metadata_df = pd.read_csv(metadata_file, sep='\t')
        for _, row in metadata_df.iterrows():
            sample_id = row['sample_kgpID']  # NA12843 format
            sample_to_ancestry[sample_id] = row['continentalGroup']
            sample_to_population[sample_id] = row['population']
    
    return sample_to_ancestry, sample_to_population

def load_rnaseq_ethnicity_data():
    """Load RNA-seq ethnicity metadata"""
    ethnicity_mapping = {}
    
    # Load the complete ethnicity mapping
    ethnicity_file = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/reports/sample_ethnicity_mapping_complete.csv"
    if os.path.exists(ethnicity_file):
        print(f"Loading RNA-seq ethnicity data from {ethnicity_file}...")
        df = pd.read_csv(ethnicity_file)
        
        for _, row in df.iterrows():
            sample_id = row['sample_id']
            ethnicity = row['ethnicity']
            dataset = row['dataset']
            
            # Map ethnicity terms to superpopulation codes
            ethnicity_to_ancestry = {
                'white': 'EUR',
                'black or african american': 'AFR', 
                'asian': 'EAS',
                'hispanic or latino': 'AMR',
                'american indian or alaska native': 'AMR',
                'multiethnic': 'MIXED',
                'more than one race': 'MIXED',
                'native hawaiian or other pacific islander': 'EAS',
                'unknown or not reported': 'UNKNOWN',
                'other': 'UNKNOWN'
            }
            
            ancestry = ethnicity_to_ancestry.get(ethnicity.lower(), 'UNKNOWN')
            ethnicity_mapping[sample_id] = {
                'ancestry': ancestry,
                'ethnicity': ethnicity,
                'dataset': dataset
            }
    
    return ethnicity_mapping

def create_comprehensive_ancestry_mapping(vcf_samples):
    """Create comprehensive ancestry mapping combining all sources"""
    print("Creating comprehensive ancestry mapping...")
    
    # Load MAGE ancestry data
    mage_ancestry, mage_population = load_mage_ancestry_data()
    print(f"Loaded {len(mage_ancestry)} samples from MAGE data")
    
    # Load RNA-seq ethnicity data
    rnaseq_ethnicity = load_rnaseq_ethnicity_data()
    print(f"Loaded {len(rnaseq_ethnicity)} samples from RNA-seq data")
    
    # Create comprehensive mapping
    comprehensive_mapping = {}
    source_counts = {'MAGE': 0, 'RNASEQ': 0, 'UNKNOWN': 0}
    
    for sample_id in vcf_samples:
        if sample_id in mage_ancestry:
            # Use MAGE data (most reliable)
            comprehensive_mapping[sample_id] = {
                'ancestry': mage_ancestry[sample_id],
                'population': mage_population.get(sample_id, 'Unknown'),
                'source': 'MAGE'
            }
            source_counts['MAGE'] += 1
            
        elif sample_id in rnaseq_ethnicity:
            # Use RNA-seq ethnicity data
            data = rnaseq_ethnicity[sample_id]
            comprehensive_mapping[sample_id] = {
                'ancestry': data['ancestry'],
                'population': data['ethnicity'],
                'source': f"RNASEQ_{data['dataset']}"
            }
            source_counts['RNASEQ'] += 1
            
        else:
            # Unknown ancestry
            comprehensive_mapping[sample_id] = {
                'ancestry': 'UNKNOWN',
                'population': 'Unknown',
                'source': 'UNKNOWN'
            }
            source_counts['UNKNOWN'] += 1
    
    # Print mapping statistics
    print(f"\nAncestry mapping results:")
    print(f"  Total VCF samples: {len(vcf_samples)}")
    print(f"  MAGE mapped: {source_counts['MAGE']} ({source_counts['MAGE']/len(vcf_samples)*100:.1f}%)")
    print(f"  RNA-seq mapped: {source_counts['RNASEQ']} ({source_counts['RNASEQ']/len(vcf_samples)*100:.1f}%)")
    print(f"  Unknown: {source_counts['UNKNOWN']} ({source_counts['UNKNOWN']/len(vcf_samples)*100:.1f}%)")
    print(f"  Total mapped: {source_counts['MAGE'] + source_counts['RNASEQ']} ({(source_counts['MAGE'] + source_counts['RNASEQ'])/len(vcf_samples)*100:.1f}%)")
    
    # Print ancestry distribution
    ancestry_counts = {}
    for mapping in comprehensive_mapping.values():
        ancestry = mapping['ancestry']
        ancestry_counts[ancestry] = ancestry_counts.get(ancestry, 0) + 1
    
    print(f"\nAncestry distribution:")
    for ancestry, count in sorted(ancestry_counts.items()):
        print(f"  {ancestry}: {count} samples ({count/len(vcf_samples)*100:.1f}%)")
    
    return comprehensive_mapping

def test_enhanced_ancestry_mapping():
    """Test the enhanced ancestry mapping"""
    print("="*60)
    print("TESTING ENHANCED ANCESTRY MAPPING")
    print("="*60)
    
    # Get sample VCF samples (first 20 for testing)
    test_samples = [
        "NA06985", "NA07000", "002_S_0413", "002_S_0729", "GTEX-1117F",
        "HG00096", "HG00100", "003_S_0907", "GTEX-111CU", "A549"
    ]
    
    mapping = create_comprehensive_ancestry_mapping(test_samples)
    
    print(f"\nSample mappings:")
    for sample_id in test_samples:
        info = mapping.get(sample_id, {})
        ancestry = info.get('ancestry', 'UNKNOWN')
        population = info.get('population', 'Unknown')
        source = info.get('source', 'UNKNOWN')
        print(f"  {sample_id}: {ancestry} ({population}) [source: {source}]")
    
    return mapping

if __name__ == "__main__":
    test_enhanced_ancestry_mapping()