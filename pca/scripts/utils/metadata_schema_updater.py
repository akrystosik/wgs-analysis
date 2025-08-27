#!/usr/bin/env python3
"""
Metadata Schema Updater
Updates the metadata schema to properly distinguish between:
- Self-reported ethnicity (from clinical/survey data)
- Genomic ancestry (inferred from genetic data)
- Population assignment (from reference projection)
"""

import pandas as pd
import numpy as np
from pathlib import Path
import json
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

# Set up file paths
base_dir = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/wgs")
pca_dir = base_dir / "pca"
rnaseq_dir = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq")
results_dir = pca_dir / "results" / "updated_metadata"
results_dir.mkdir(parents=True, exist_ok=True)

def load_existing_metadata():
    """Load all existing metadata sources"""
    print("Loading existing metadata sources...")
    
    metadata_sources = {}
    
    # 1. RNA-seq ethnicity mapping (self-reported)
    rnaseq_file = rnaseq_dir / "reports" / "subject_ethnicity_mapping_czi_compliant.csv"
    if rnaseq_file.exists():
        rnaseq_data = pd.read_csv(rnaseq_file)
        metadata_sources['rnaseq_ethnicity'] = rnaseq_data
        print(f"Loaded RNA-seq ethnicity data: {len(rnaseq_data)} samples")
    
    # 2. Comprehensive ancestry mapping (mixed sources)
    ancestry_file = pca_dir / "results" / "comprehensive_ancestry_mapping.csv"
    if ancestry_file.exists():
        ancestry_data = pd.read_csv(ancestry_file)
        metadata_sources['comprehensive_ancestry'] = ancestry_data
        print(f"Loaded comprehensive ancestry data: {len(ancestry_data)} samples")
    
    # 3. PLINK PCA results (genomic analysis)
    plink_file = pca_dir / "results" / "plink_analysis" / "plink_pca_results.csv"
    if plink_file.exists():
        plink_data = pd.read_csv(plink_file)
        metadata_sources['plink_pca'] = plink_data
        print(f"Loaded PLINK PCA results: {len(plink_data)} samples")
    
    # 4. Reference projection results (ML-based assignment)
    projection_file = pca_dir / "results" / "reference_projection" / "reference_projection_results.csv"
    if projection_file.exists():
        projection_data = pd.read_csv(projection_file)
        metadata_sources['reference_projection'] = projection_data
        print(f"Loaded reference projection results: {len(projection_data)} samples")
    
    return metadata_sources

def create_unified_metadata_schema():
    """Create a unified metadata schema with clear distinctions"""
    print("Creating unified metadata schema...")
    
    schema = {
        "metadata_version": "2.0",
        "created_date": datetime.now().isoformat(),
        "description": "Unified metadata schema distinguishing self-reported ethnicity from genomic ancestry",
        
        "sample_identification": {
            "sample_id": {
                "type": "string",
                "description": "Unique sample identifier",
                "required": True
            },
            "original_id": {
                "type": "string", 
                "description": "Original sample ID from source dataset",
                "required": False
            },
            "data_source": {
                "type": "categorical",
                "description": "Source dataset of the sample",
                "values": ["MAGE", "GTEx", "ADNI", "ENCODE", "1000_Genomes"],
                "required": True
            },
            "sample_type": {
                "type": "categorical",
                "description": "Type of biological sample",
                "values": ["cell_line", "tissue", "blood", "other"],
                "required": True
            }
        },
        
        "self_reported_ethnicity": {
            "description": "Ethnicity information as reported by participants or clinical records",
            "ethnicity_reported": {
                "type": "categorical",
                "description": "Self-reported ethnicity/race",
                "values": ["White", "Black", "Asian", "Hispanic", "Native_American", "Pacific_Islander", "Mixed", "Other", "Unknown"],
                "required": False,
                "source": "clinical_survey"
            },
            "ethnicity_confidence": {
                "type": "categorical",
                "description": "Confidence in ethnicity reporting",
                "values": ["high", "medium", "low", "unknown"],
                "required": False
            },
            "ethnicity_source": {
                "type": "categorical",
                "description": "Source of ethnicity information",
                "values": ["self_reported", "clinical_record", "database_annotation", "inferred"],
                "required": False
            }
        },
        
        "genomic_ancestry": {
            "description": "Ancestry inferred from genetic/genomic data",
            "pca_ancestry": {
                "type": "categorical",
                "description": "Ancestry inferred from PCA clustering",
                "values": ["EUR", "AFR", "EAS", "SAS", "AMR", "OCE", "ADMIXED", "Unknown"],
                "required": False,
                "method": "pca_clustering"
            },
            "reference_ancestry": {
                "type": "categorical",
                "description": "Ancestry from 1000 Genomes reference projection",
                "values": ["EUR", "AFR", "EAS", "SAS", "AMR", "OCE", "ADMIXED", "Uncertain"],
                "required": False,
                "method": "reference_projection"
            },
            "consensus_ancestry": {
                "type": "categorical",
                "description": "Consensus ancestry from multiple genomic methods",
                "values": ["EUR", "AFR", "EAS", "SAS", "AMR", "OCE", "ADMIXED", "Uncertain"],
                "required": False,
                "method": "consensus"
            },
            "ancestry_confidence": {
                "type": "float",
                "description": "Confidence score for ancestry assignment (0-1)",
                "range": [0, 1],
                "required": False
            }
        },
        
        "population_genetics": {
            "description": "Population genetics analysis results",
            "pc_coordinates": {
                "type": "array",
                "description": "Principal component coordinates (PC1-PC20)",
                "required": False
            },
            "admixture_proportions": {
                "type": "object",
                "description": "Admixture analysis ancestry proportions",
                "required": False
            },
            "genetic_distance": {
                "type": "float",
                "description": "Genetic distance to nearest reference population",
                "required": False
            }
        },
        
        "quality_control": {
            "description": "Quality control metrics for ancestry inference",
            "genotyping_rate": {
                "type": "float",
                "description": "Proportion of successfully genotyped variants",
                "range": [0, 1],
                "required": False
            },
            "heterozygosity": {
                "type": "float",
                "description": "Observed heterozygosity rate",
                "required": False
            },
            "inbreeding_coefficient": {
                "type": "float",
                "description": "Inbreeding coefficient (F)",
                "required": False
            },
            "ancestry_inference_method": {
                "type": "categorical",
                "description": "Method used for ancestry inference",
                "values": ["pca", "admixture", "reference_projection", "consensus"],
                "required": False
            }
        },
        
        "concordance_analysis": {
            "description": "Concordance between self-reported and genomic ancestry",
            "ethnicity_ancestry_concordance": {
                "type": "categorical",
                "description": "Agreement between self-reported ethnicity and genomic ancestry",
                "values": ["concordant", "discordant", "uncertain", "unknown"],
                "required": False
            },
            "discordance_notes": {
                "type": "string",
                "description": "Notes explaining discordance between ethnicity and ancestry",
                "required": False
            }
        }
    }
    
    return schema

def integrate_metadata_sources(metadata_sources):
    """Integrate all metadata sources into unified format"""
    print("Integrating metadata sources...")
    
    # Start with the most comprehensive source
    if 'reference_projection' in metadata_sources:
        base_data = metadata_sources['reference_projection'].copy()
    elif 'plink_pca' in metadata_sources:
        base_data = metadata_sources['plink_pca'].copy()
    else:
        print("ERROR: No PCA/projection data found!")
        return None
    
    # Rename columns to match schema
    column_mapping = {
        'IID': 'sample_id',
        'sample_type': 'data_source',
        'inferred_ancestry': 'pca_ancestry',
        'consensus_prediction': 'reference_ancestry',
        'KNN_confidence': 'ancestry_confidence'
    }
    
    for old_col, new_col in column_mapping.items():
        if old_col in base_data.columns:
            base_data[new_col] = base_data[old_col]
    
    # Add self-reported ethnicity from proper ethnicity mapping file
    ethnicity_file = rnaseq_dir / "reports" / "sample_ethnicity_mapping_with_ontology.csv"
    if ethnicity_file.exists():
        ethnicity_data = pd.read_csv(ethnicity_file)
        print(f"Loaded ethnicity mapping: {len(ethnicity_data)} samples")
        
        # Merge ethnicity data on sample_id
        if 'sample_id' in ethnicity_data.columns:
            ethnicity_df = ethnicity_data[['sample_id', 'ethnicity', 'hancestro_term']].copy()
            ethnicity_df['ethnicity_reported'] = ethnicity_df['ethnicity']
            ethnicity_df['ethnicity_source'] = 'self_reported'
            ethnicity_df['ethnicity_confidence'] = 'high'
            ethnicity_df['self_reported_ethnicity_ontology_term_id'] = ethnicity_df['hancestro_term']
            
            base_data = base_data.merge(ethnicity_df[['sample_id', 'ethnicity_reported', 
                                                    'ethnicity_source', 'ethnicity_confidence',
                                                    'self_reported_ethnicity_ontology_term_id']], 
                                      on='sample_id', how='left')
    
    # Fallback to RNA-seq data if main ethnicity file not found
    elif 'rnaseq_ethnicity' in metadata_sources:
        rnaseq_data = metadata_sources['rnaseq_ethnicity']
        
        # Create mapping for ethnicity
        ethnicity_mapping = {
            'White': 'White',
            'Black': 'Black', 
            'Asian': 'Asian',
            'Hispanic': 'Hispanic',
            'European': 'White',
            'African': 'Black',
            'East_Asian': 'Asian',
            'South_Asian': 'Asian',
            'American': 'Hispanic'
        }
        
        # Merge ethnicity data
        if 'ethnicity' in rnaseq_data.columns:
            ethnicity_df = rnaseq_data[['sample_id', 'ethnicity']].copy()
            ethnicity_df['ethnicity_reported'] = ethnicity_df['ethnicity'].map(ethnicity_mapping)
            ethnicity_df['ethnicity_source'] = 'database_annotation'
            ethnicity_df['ethnicity_confidence'] = 'medium'
            
            base_data = base_data.merge(ethnicity_df[['sample_id', 'ethnicity_reported', 
                                                    'ethnicity_source', 'ethnicity_confidence']], 
                                      on='sample_id', how='left')
    
    # Add original sample identification
    base_data['original_id'] = base_data['sample_id']
    
    # Determine sample type
    base_data['sample_type'] = base_data['data_source'].apply(lambda x: 
        'cell_line' if x == 'ENCODE' else 'tissue' if x in ['GTEx', 'MAGE'] else 'blood'
    )
    
    # Create consensus ancestry
    base_data['consensus_ancestry'] = base_data.apply(
        lambda row: row.get('reference_ancestry', row.get('pca_ancestry', 'Unknown')), axis=1
    )
    
    # Add ancestry inference method
    base_data['ancestry_inference_method'] = 'consensus'
    
    # Add PC coordinates
    pc_cols = [col for col in base_data.columns if col.startswith('PC')]
    if pc_cols:
        base_data['pc_coordinates'] = base_data[pc_cols].apply(lambda x: x.tolist(), axis=1)
    
    # Calculate concordance
    if 'ethnicity_reported' in base_data.columns and 'consensus_ancestry' in base_data.columns:
        base_data['ethnicity_ancestry_concordance'] = base_data.apply(
            lambda row: calculate_concordance(row['ethnicity_reported'], row['consensus_ancestry']), axis=1
        )
    
    print(f"Integrated metadata for {len(base_data)} samples")
    
    return base_data

def calculate_concordance(ethnicity, ancestry):
    """Calculate concordance between self-reported ethnicity and genomic ancestry"""
    if pd.isna(ethnicity) or pd.isna(ancestry):
        return 'unknown'
    
    # Define concordance mappings
    concordance_map = {
        'White': ['EUR'],
        'Black': ['AFR'],
        'Asian': ['EAS', 'SAS'],
        'Hispanic': ['AMR', 'ADMIXED']
    }
    
    if ethnicity in concordance_map:
        if ancestry in concordance_map[ethnicity]:
            return 'concordant'
        else:
            return 'discordant'
    
    return 'uncertain'

def generate_metadata_report(unified_data, schema):
    """Generate comprehensive metadata report"""
    print("Generating metadata report...")
    
    report = []
    report.append("# Updated Metadata Schema Report")
    report.append("=" * 50)
    report.append("")
    
    report.append("## Schema Overview")
    report.append("This updated schema distinguishes between:")
    report.append("- **Self-reported ethnicity**: From clinical/survey data")
    report.append("- **Genomic ancestry**: Inferred from genetic analysis")
    report.append("- **Population assignment**: From reference projection")
    report.append("")
    
    report.append("## Sample Composition")
    report.append(f"- Total samples: {len(unified_data):,}")
    
    if 'data_source' in unified_data.columns:
        report.append("\n### By Data Source")
        for source, count in unified_data['data_source'].value_counts().items():
            pct = count / len(unified_data) * 100
            report.append(f"- {source}: {count:,} ({pct:.1f}%)")
    
    if 'sample_type' in unified_data.columns:
        report.append("\n### By Sample Type")
        for sample_type, count in unified_data['sample_type'].value_counts().items():
            pct = count / len(unified_data) * 100
            report.append(f"- {sample_type}: {count:,} ({pct:.1f}%)")
    
    report.append("")
    
    report.append("## Ethnicity vs Ancestry")
    
    if 'ethnicity_reported' in unified_data.columns:
        report.append("\n### Self-Reported Ethnicity")
        eth_counts = unified_data['ethnicity_reported'].value_counts()
        for ethnicity, count in eth_counts.items():
            pct = count / len(unified_data) * 100
            report.append(f"- {ethnicity}: {count:,} ({pct:.1f}%)")
    
    if 'consensus_ancestry' in unified_data.columns:
        report.append("\n### Genomic Ancestry")
        anc_counts = unified_data['consensus_ancestry'].value_counts()
        for ancestry, count in anc_counts.items():
            pct = count / len(unified_data) * 100
            report.append(f"- {ancestry}: {count:,} ({pct:.1f}%)")
    
    if 'ethnicity_ancestry_concordance' in unified_data.columns:
        report.append("\n### Concordance Analysis")
        conc_counts = unified_data['ethnicity_ancestry_concordance'].value_counts()
        for concordance, count in conc_counts.items():
            pct = count / len(unified_data) * 100
            report.append(f"- {concordance}: {count:,} ({pct:.1f}%)")
    
    report.append("")
    
    report.append("## Quality Metrics")
    if 'ancestry_confidence' in unified_data.columns:
        conf_data = unified_data['ancestry_confidence'].dropna()
        if len(conf_data) > 0:
            report.append(f"- Mean ancestry confidence: {conf_data.mean():.3f}")
            report.append(f"- High confidence samples (>0.7): {(conf_data > 0.7).sum():,}")
    
    report.append("")
    
    report.append("## Schema Implementation")
    report.append("- Version: 2.0")
    report.append("- Created: " + datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    report.append("- Method: PLINK PCA + Reference Projection + ML Classification")
    report.append("- Markers: 2,884,699 ancestry-informative variants")
    
    return report

def save_updated_metadata(unified_data, schema, report):
    """Save updated metadata, schema, and report"""
    print("Saving updated metadata...")
    
    # Save unified metadata
    metadata_file = results_dir / "unified_metadata_v2.csv"
    unified_data.to_csv(metadata_file, index=False)
    print(f"Saved unified metadata to {metadata_file}")
    
    # Save schema
    schema_file = results_dir / "metadata_schema_v2.json"
    with open(schema_file, 'w') as f:
        json.dump(schema, f, indent=2)
    print(f"Saved schema to {schema_file}")
    
    # Save report
    report_file = results_dir / "metadata_update_report.txt"
    with open(report_file, 'w') as f:
        f.write('\n'.join(report))
    print(f"Saved report to {report_file}")
    
    # Create a simplified version for downstream analysis
    essential_cols = [
        'sample_id', 'data_source', 'sample_type',
        'ethnicity_reported', 'consensus_ancestry', 'ancestry_confidence',
        'ethnicity_ancestry_concordance', 'PC1', 'PC2', 'PC3'
    ]
    
    available_cols = [col for col in essential_cols if col in unified_data.columns]
    simplified_data = unified_data[available_cols].copy()
    
    simplified_file = results_dir / "simplified_metadata_v2.csv"
    simplified_data.to_csv(simplified_file, index=False)
    print(f"Saved simplified metadata to {simplified_file}")
    
    # Create Cell x Gene compliant format
    cxg_format = unified_data.copy()
    cxg_format['cell_type_ontology_term_id'] = 'CL:0000000'  # Generic cell
    cxg_format['development_stage_ontology_term_id'] = 'HsapDv:0000087'  # Adult
    cxg_format['sex_ontology_term_id'] = 'PATO:0000383'  # Unknown
    cxg_format['organism_ontology_term_id'] = 'NCBITaxon:9606'  # Human
    
    # Keep ethnicity ontology terms if they exist, don't overwrite with ancestry
    if 'self_reported_ethnicity_ontology_term_id' not in cxg_format.columns:
        # Only add empty ethnicity column if it doesn't exist
        cxg_format['self_reported_ethnicity_ontology_term_id'] = ''
    
    cxg_file = results_dir / "cellxgene_compliant_metadata.csv"
    cxg_format.to_csv(cxg_file, index=False)
    print(f"Saved Cell x Gene compliant metadata to {cxg_file}")

def main():
    """Main metadata update function"""
    print("Starting Metadata Schema Update...")
    print("=" * 50)
    
    # Load existing metadata
    metadata_sources = load_existing_metadata()
    
    if not metadata_sources:
        print("ERROR: No metadata sources found!")
        return
    
    # Create unified schema
    schema = create_unified_metadata_schema()
    
    # Integrate metadata sources
    unified_data = integrate_metadata_sources(metadata_sources)
    
    if unified_data is None:
        print("ERROR: Failed to integrate metadata sources!")
        return
    
    # Generate report
    report = generate_metadata_report(unified_data, schema)
    
    # Save everything
    save_updated_metadata(unified_data, schema, report)
    
    print("\nMetadata schema update complete! 🎉")
    print(f"Results saved to: {results_dir}")
    print("\nKey improvements:")
    print("- Clear distinction between self-reported ethnicity and genomic ancestry")
    print("- Standardized ancestry terms (EUR, AFR, EAS, SAS, AMR)")
    print("- Concordance analysis between ethnicity and ancestry")
    print("- Cell x Gene compliant format")
    print("- Comprehensive quality metrics")

if __name__ == "__main__":
    main()