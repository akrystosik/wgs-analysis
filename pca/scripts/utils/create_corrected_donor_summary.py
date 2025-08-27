#!/usr/bin/env python3
"""
Create Corrected Donor-Level Summary
Fixes the tissue duplication issue and creates a proper donor-level summary
"""

import pandas as pd
import numpy as np
from pathlib import Path

# Set up file paths
base_dir = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/wgs")
pca_dir = base_dir / "pca"
rnaseq_dir = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq")
results_dir = pca_dir / "results" / "updated_metadata"
results_dir.mkdir(parents=True, exist_ok=True)

def create_corrected_summary():
    """Create corrected donor-level summary from existing results"""
    print("Creating corrected donor-level summary...")
    
    # 1. Load the correct donor-level reference projection results (2,331 donors)
    ref_proj_file = pca_dir / "results" / "reference_projection" / "reference_projection_results.csv"
    donor_data = pd.read_csv(ref_proj_file)
    print(f"Loaded reference projection data: {len(donor_data)} unique donors")
    
    # Fix GM23248 classification (it's ENCODE, not MAGE)
    gm_mask = donor_data['IID'] == 'GM23248'
    if gm_mask.any():
        donor_data.loc[gm_mask, 'data_source'] = 'ENCODE'
        donor_data.loc[gm_mask, 'sample_type'] = 'cell_line'
        print("Fixed GM23248: reclassified from MAGE to ENCODE")
    
    # 2. Load ethnicity mapping (donor-level)
    ethnicity_file = rnaseq_dir / "reports" / "sample_ethnicity_mapping_with_ontology.csv"
    if ethnicity_file.exists():
        ethnicity_data = pd.read_csv(ethnicity_file)
        print(f"Loaded ethnicity mapping: {len(ethnicity_data)} entries")
        
        # Deduplicate ethnicity data to donor level (in case it has tissue duplicates)
        ethnicity_donor = ethnicity_data.drop_duplicates(subset=['sample_id']).copy()
        print(f"Deduplicated to {len(ethnicity_donor)} unique donors with ethnicity data")
        
        # Merge ethnicity data with donor results
        merged_data = donor_data.merge(
            ethnicity_donor[['sample_id', 'ethnicity', 'hancestro_term']],
            left_on='IID',
            right_on='sample_id',
            how='left',
            suffixes=('', '_ethnicity')
        )
        
        # Rename columns for clarity
        merged_data['self_reported_ethnicity'] = merged_data['ethnicity']
        merged_data['ethnicity_ontology_term'] = merged_data['hancestro_term']
        
    else:
        merged_data = donor_data.copy()
        merged_data['self_reported_ethnicity'] = None
        merged_data['ethnicity_ontology_term'] = None
        print("Warning: No ethnicity mapping found")
    
    # 3. Handle MAGE vs non-MAGE samples differently
    # MAGE samples have 1000 Genomes population labels, not self-reported ethnicity
    merged_data['mage_1kgp_population'] = None
    merged_data['mage_1kgp_ancestry'] = None
    
    # For MAGE samples, use their 1000 Genomes population labels
    mage_mask = merged_data['data_source'] == 'MAGE'
    if 'population' in merged_data.columns:
        merged_data.loc[mage_mask, 'mage_1kgp_population'] = merged_data.loc[mage_mask, 'population']
    if 'ancestry' in merged_data.columns:
        merged_data.loc[mage_mask, 'mage_1kgp_ancestry'] = merged_data.loc[mage_mask, 'ancestry']
    
    # Clear self-reported ethnicity for MAGE samples (they don't have self-reports)
    merged_data.loc[mage_mask, 'self_reported_ethnicity'] = None
    merged_data.loc[mage_mask, 'ethnicity_ontology_term'] = None
    
    # 4. Create corrected summary with essential columns
    summary_columns = [
        'IID',                           # Donor ID
        'data_source',                   # Dataset (MAGE, GTEx, ADNI, ENCODE)
        'sample_type',                   # Sample type
        'self_reported_ethnicity',       # Self-reported ethnicity (GTEx, ADNI only)
        'ethnicity_ontology_term',       # HANCESTRO term (GTEx, ADNI only)
        'mage_1kgp_population',          # 1000 Genomes population (MAGE only)
        'mage_1kgp_ancestry',            # 1000 Genomes ancestry (MAGE only)
        'inferred_ancestry',            # PLINK PCA ancestry
        'consensus_prediction',         # ML consensus prediction
        'KNN_confidence',               # ML confidence score
        'PC1', 'PC2', 'PC3'            # Population structure coordinates
    ]
    
    # Select available columns
    available_cols = [col for col in summary_columns if col in merged_data.columns]
    corrected_summary = merged_data[available_cols].copy()
    
    # Rename columns for partner clarity
    column_mapping = {
        'IID': 'donor_id',
        'consensus_prediction': 'genomic_ancestry_ml',
        'inferred_ancestry': 'genomic_ancestry_pca',
        'KNN_confidence': 'ml_confidence',
        'data_source': 'dataset'
    }
    
    corrected_summary = corrected_summary.rename(columns=column_mapping)
    
    # 4. Add concordance analysis
    if 'self_reported_ethnicity' in corrected_summary.columns and 'genomic_ancestry_ml' in corrected_summary.columns:
        corrected_summary['ethnicity_ancestry_concordance'] = corrected_summary.apply(
            lambda row: calculate_concordance(row['self_reported_ethnicity'], row['genomic_ancestry_ml']), 
            axis=1
        )
    
    # 5. Create summary statistics
    print("\n=== CORRECTED DONOR-LEVEL SUMMARY ===")
    print(f"Total unique donors: {len(corrected_summary):,}")
    
    if 'dataset' in corrected_summary.columns:
        print("\nDataset Distribution:")
        dataset_counts = corrected_summary['dataset'].value_counts()
        for dataset, count in dataset_counts.items():
            pct = count / len(corrected_summary) * 100
            print(f"  {dataset}: {count:,} ({pct:.1f}%)")
    
    if 'genomic_ancestry_ml' in corrected_summary.columns:
        print("\nGenomic Ancestry (ML) Distribution:")
        ancestry_counts = corrected_summary['genomic_ancestry_ml'].value_counts()
        for ancestry, count in ancestry_counts.items():
            pct = count / len(corrected_summary) * 100
            print(f"  {ancestry}: {count:,} ({pct:.1f}%)")
    
    if 'ethnicity_ancestry_concordance' in corrected_summary.columns:
        print("\nEthnicity-Ancestry Concordance:")
        concordance_counts = corrected_summary['ethnicity_ancestry_concordance'].value_counts()
        for concordance, count in concordance_counts.items():
            pct = count / len(corrected_summary) * 100
            print(f"  {concordance}: {count:,} ({pct:.1f}%)")
    
    return corrected_summary

def calculate_concordance(ethnicity, ancestry):
    """Calculate concordance between self-reported ethnicity and genomic ancestry"""
    if pd.isna(ethnicity) or pd.isna(ancestry):
        return 'unknown'
    
    # Define concordance mappings
    concordance_map = {
        'white': ['EUR'],
        'black or african american': ['AFR'],
        'asian': ['EAS', 'SAS'],
        'hispanic or latino': ['AMR'],
        'american indian or alaska native': ['AMR']
    }
    
    ethnicity_lower = str(ethnicity).lower()
    
    for eth_key, anc_values in concordance_map.items():
        if eth_key in ethnicity_lower:
            if ancestry in anc_values:
                return 'concordant'
            else:
                return 'discordant'
    
    return 'uncertain'

def main():
    """Main function to create corrected summary"""
    print("Creating Corrected Donor-Level Summary")
    print("=" * 50)
    
    # Create corrected summary
    corrected_summary = create_corrected_summary()
    
    # Save corrected summary
    output_file = results_dir / "corrected_donor_summary.csv"
    corrected_summary.to_csv(output_file, index=False)
    print(f"\nSaved corrected donor summary to: {output_file}")
    
    # Create partner-friendly version
    partner_columns = [
        'donor_id', 'dataset', 
        'self_reported_ethnicity', 'ethnicity_ontology_term',  # GTEx, ADNI only
        'mage_1kgp_population', 'mage_1kgp_ancestry',          # MAGE only
        'genomic_ancestry_pca', 'genomic_ancestry_ml', 'ml_confidence', 
        'ethnicity_ancestry_concordance', 'PC1', 'PC2'
    ]
    
    available_partner_cols = [col for col in partner_columns if col in corrected_summary.columns]
    partner_summary = corrected_summary[available_partner_cols].copy()
    
    partner_file = results_dir / "partner_donor_summary.csv"
    partner_summary.to_csv(partner_file, index=False)
    print(f"Saved partner summary to: {partner_file}")
    
    print(f"\n✅ Corrected summaries created!")
    print(f"   - Full summary: {len(corrected_summary)} unique donors")
    print(f"   - Partner summary: {len(partner_summary)} donors with key columns")
    print(f"\n🚫 Previous files with 20,942 rows were tissue-duplicated")
    print(f"✅ New files have {len(corrected_summary)} unique donors (correct for population genetics)")
    print(f"\n📊 Data Structure:")
    print(f"   - MAGE samples: 1000 Genomes population labels (not self-reported)")
    print(f"   - GTEx/ADNI: Self-reported ethnicity with HANCESTRO terms")
    print(f"   - All samples: Genomic ancestry from PCA + ML hybrid approach")

if __name__ == "__main__":
    main()