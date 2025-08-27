#!/usr/bin/env python3
"""
Improved ancestry loader that addresses the 31.4% mapping bottleneck
"""

import pandas as pd
import numpy as np
import os

def load_1000g_ancestry_data():
    """Load 1000 Genomes ancestry data"""
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
    
    # Load MAGE metadata for additional mappings
    metadata_file = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/mage/sample.metadata.MAGE.v1.0.txt"
    if os.path.exists(metadata_file):
        metadata_df = pd.read_csv(metadata_file, sep='\t')
        for _, row in metadata_df.iterrows():
            sample_id = row['sample_kgpID']  # NA12843 format
            sample_to_ancestry[sample_id] = row['continentalGroup']
            sample_to_population[sample_id] = row['population']
    
    return sample_to_ancestry, sample_to_population

def load_adni_ethnicity_mapping():
    """Load ADNI ethnicity mapping from CZI-compliant CSV"""
    adni_mapping = {}
    
    ethnicity_file = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/reports/subject_ethnicity_mapping_czi_compliant.csv"
    if os.path.exists(ethnicity_file):
        df = pd.read_csv(ethnicity_file)
        
        # HANCESTRO to superpopulation mapping
        hancestro_to_superpop = {
            'HANCESTRO:0005': 'EUR',  # white
            'HANCESTRO:0010': 'AFR',  # black or african american
            'HANCESTRO:0014': 'EAS',  # asian (if present)
            'HANCESTRO:0016': 'SAS',  # south asian (if present)
            'HANCESTRO:0019': 'AMR',  # hispanic or latino (if present)
            'unknown': 'UNKNOWN'
        }
        
        for _, row in df.iterrows():
            sample_id = row['subject_id']
            hancestro_term = row['self_reported_ethnicity_ontology_term_id']
            
            # Map to superpopulation
            ancestry = hancestro_to_superpop.get(hancestro_term, 'UNKNOWN')
            adni_mapping[sample_id] = {
                'ancestry': ancestry,
                'population': row['ethnicity_label_for_reference'],
                'source': 'ADNI_ETHNICITY',
                'hancestro_term': hancestro_term
            }
    
    return adni_mapping

def create_comprehensive_ancestry_mapping(vcf_samples):
    """Create comprehensive ancestry mapping for all VCF samples"""
    
    print("Creating comprehensive ancestry mapping...")
    
    # Initialize mapping dictionary
    ancestry_mapping = {}
    
    # Load 1000 Genomes data
    kg_ancestry, kg_population = load_1000g_ancestry_data()
    print(f"Loaded 1000 Genomes data for {len(kg_ancestry)} samples")
    
    # Load ADNI ethnicity data
    adni_mapping = load_adni_ethnicity_mapping()
    print(f"Loaded ADNI ethnicity data for {len(adni_mapping)} samples")
    
    # Process each VCF sample
    mapped_count = 0
    
    for sample_id in vcf_samples:
        if sample_id in kg_ancestry:
            # 1000 Genomes sample (HG, GM, NA prefixes)
            ancestry_mapping[sample_id] = {
                'ancestry': kg_ancestry[sample_id],
                'population': kg_population[sample_id],
                'source': 'MAGE_1000G',
                'confidence': 'high'
            }
            mapped_count += 1
            
        elif sample_id in adni_mapping:
            # ADNI sample with ethnicity mapping
            ancestry_mapping[sample_id] = adni_mapping[sample_id]
            ancestry_mapping[sample_id]['confidence'] = 'medium'
            mapped_count += 1
            
        else:
            # Unknown sample (likely ENCODE cell lines)
            ancestry_mapping[sample_id] = {
                'ancestry': 'UNKNOWN',
                'population': 'UNKNOWN',
                'source': 'UNKNOWN',
                'confidence': 'low'
            }
    
    print(f"Successfully mapped {mapped_count}/{len(vcf_samples)} samples ({mapped_count/len(vcf_samples)*100:.1f}%)")
    
    # Print mapping statistics
    ancestry_counts = {}
    source_counts = {}
    
    for sample_id, data in ancestry_mapping.items():
        ancestry = data['ancestry']
        source = data['source']
        
        ancestry_counts[ancestry] = ancestry_counts.get(ancestry, 0) + 1
        source_counts[source] = source_counts.get(source, 0) + 1
    
    print("\nAncestry distribution:")
    for ancestry, count in sorted(ancestry_counts.items()):
        print(f"  {ancestry}: {count} samples")
    
    print("\nData source distribution:")
    for source, count in sorted(source_counts.items()):
        print(f"  {source}: {count} samples")
    
    return ancestry_mapping

def save_ancestry_mapping_report(ancestry_mapping, output_file):
    """Save detailed ancestry mapping report"""
    
    # Convert to DataFrame for analysis
    mapping_data = []
    for sample_id, data in ancestry_mapping.items():
        mapping_data.append({
            'sample_id': sample_id,
            'ancestry': data['ancestry'],
            'population': data['population'],
            'source': data['source'],
            'confidence': data['confidence']
        })
    
    df = pd.DataFrame(mapping_data)
    df.to_csv(output_file, index=False)
    print(f"Ancestry mapping report saved to: {output_file}")
    
    return df

if __name__ == "__main__":
    # Test the improved mapping
    import subprocess
    
    print("Testing improved ancestry mapping...")
    
    # Load VCF samples using bcftools
    vcf_file = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/wgs/data/variants/all_vcfs/merged.all.biallelic.maf0.01.vcf.gz"
    result = subprocess.run(['bcftools', 'query', '-l', vcf_file], capture_output=True, text=True)
    samples = result.stdout.strip().split('\n')
    
    print(f"Total VCF samples: {len(samples)}")
    
    # Create comprehensive mapping
    ancestry_mapping = create_comprehensive_ancestry_mapping(samples)
    
    # Save report
    output_file = "../../results/comprehensive_ancestry_mapping.csv"
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    mapping_df = save_ancestry_mapping_report(ancestry_mapping, output_file)
    
    print(f"\\nImproved ancestry mapping completed!")
    print(f"Mapping success rate: {len([s for s in ancestry_mapping.values() if s['ancestry'] != 'UNKNOWN'])/len(ancestry_mapping)*100:.1f}%")