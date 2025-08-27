#!/usr/bin/env python3
"""
Reference Population Projection Analysis
Projects ADNI samples onto reference population PC space using MAGE samples
with known 1000 Genomes ancestry as reference populations.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, confusion_matrix
import warnings
warnings.filterwarnings('ignore')

# Set up file paths
base_dir = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/wgs")
pca_dir = base_dir / "pca"
results_dir = pca_dir / "results" / "reference_projection"
results_dir.mkdir(parents=True, exist_ok=True)

def load_data():
    """Load PCA results and ancestry metadata"""
    print("Loading PCA results and ancestry metadata...")
    
    # Load PCA results from Stage 1 (updated with reference-guided clustering)
    pca_file = pca_dir / "results" / "plink_analysis" / "plink_pca_results.csv"
    pca_data = pd.read_csv(pca_file)
    
    # Load donor summary with comprehensive ancestry information (only metadata columns)
    donor_file = pca_dir / "results" / "updated_metadata" / "donor_summary.csv"
    donor_data = pd.read_csv(donor_file)
    
    # Select only metadata columns from donor_data (avoid PC column conflicts)
    metadata_cols = ['donor_id', 'dataset', 'self_reported_ethnicity', 'ethnicity_ontology_term', 
                    'mage_1kgp_population', 'mage_1kgp_ancestry', 'genomic_ancestry_pca', 
                    'genomic_ancestry_ml', 'ml_confidence', 'ethnicity_ancestry_concordance']
    donor_metadata = donor_data[metadata_cols].copy()
    
    # Merge PCA with donor metadata using IID/donor_id
    merged_data = pca_data.merge(donor_metadata, left_on='IID', right_on='donor_id', how='left')
    
    print(f"Loaded {len(merged_data)} samples with PCA coordinates (PC1-PC20)")
    print(f"MAGE reference samples with 1KGP ancestry: {(merged_data['dataset'] == 'MAGE').sum()}")
    print(f"Samples with Stage 1 ancestry assignments: {(merged_data['inferred_ancestry'] != 'Unknown').sum()}")
    print(f"Unknown samples for Stage 2 ML: {(merged_data['inferred_ancestry'] == 'Unknown').sum()}")
    
    return merged_data

def create_reference_populations(data):
    """Create reference populations from MAGE samples with known ancestry"""
    print("Creating reference populations...")
    
    # Filter for MAGE samples with known 1000 Genomes ancestry
    reference_samples = data[
        (data['dataset'] == 'MAGE') & 
        (data['mage_1kgp_ancestry'].notna()) & 
        (data['mage_1kgp_ancestry'] != 'Unknown')
    ].copy()
    
    # Use the standardized 1000 Genomes ancestry labels
    reference_samples['ancestry_standard'] = reference_samples['mage_1kgp_ancestry']
    
    # Filter for major populations with sufficient samples (>= 20 for robust ML training)
    ancestry_counts = reference_samples['ancestry_standard'].value_counts()
    major_ancestries = ancestry_counts[ancestry_counts >= 20].index.tolist()
    
    reference_samples = reference_samples[reference_samples['ancestry_standard'].isin(major_ancestries)]
    
    print(f"Reference population composition (1000 Genomes ancestry):")
    for ancestry, count in reference_samples['ancestry_standard'].value_counts().items():
        print(f"  {ancestry}: {count} samples")
    
    return reference_samples

def train_ancestry_classifiers(reference_samples):
    """Train ML classifiers on reference populations"""
    print("Training ancestry classifiers...")
    
    # Prepare features (PC coordinates)
    pc_cols = [col for col in reference_samples.columns if col.startswith('PC')]
    X_ref = reference_samples[pc_cols].values
    y_ref = reference_samples['ancestry_standard'].values
    
    # Train classifiers
    classifiers = {}
    
    # K-Nearest Neighbors
    knn = KNeighborsClassifier(n_neighbors=5, weights='distance')
    knn.fit(X_ref, y_ref)
    classifiers['KNN'] = knn
    
    # Random Forest
    rf = RandomForestClassifier(n_estimators=100, random_state=42)
    rf.fit(X_ref, y_ref)
    classifiers['RandomForest'] = rf
    
    # Cross-validation on reference samples
    from sklearn.model_selection import cross_val_score
    print("Classifier performance on reference samples:")
    for name, clf in classifiers.items():
        cv_scores = cross_val_score(clf, X_ref, y_ref, cv=5)
        print(f"  {name}: {cv_scores.mean():.3f} ± {cv_scores.std():.3f}")
    
    return classifiers, pc_cols

def project_unknown_samples(data, reference_samples, classifiers, pc_cols):
    """Project unknown samples (Stage 1 'Unknown') onto reference population space"""
    print("Projecting unknown samples from Stage 1...")
    
    # Identify samples that Stage 1 couldn't classify (inferred_ancestry == 'Unknown')
    unknown_samples = data[data['inferred_ancestry'] == 'Unknown'].copy()
    
    if len(unknown_samples) == 0:
        print("No unknown samples to project - Stage 1 classified all samples")
        return data
    
    # Get PC coordinates for unknown samples
    X_unknown = unknown_samples[pc_cols].values
    
    # Project using each classifier
    projections = {}
    for name, clf in classifiers.items():
        predictions = clf.predict(X_unknown)
        probabilities = clf.predict_proba(X_unknown)
        
        # Store predictions
        projections[f'{name}_prediction'] = predictions
        projections[f'{name}_confidence'] = probabilities.max(axis=1)
    
    # Add projections to unknown samples
    for col, values in projections.items():
        unknown_samples[col] = values
    
    # Create consensus prediction with confidence threshold
    unknown_samples['ml_consensus_prediction'] = unknown_samples.apply(
        lambda row: row['KNN_prediction'] if row['KNN_confidence'] > 0.7 else 
                   row['RandomForest_prediction'] if row['RandomForest_confidence'] > 0.7 else 
                   'Uncertain', axis=1
    )
    
    # Calculate overall ML confidence (average of KNN and RF)
    unknown_samples['ml_consensus_confidence'] = (
        unknown_samples['KNN_confidence'] + unknown_samples['RandomForest_confidence']
    ) / 2
    
    # Update the main dataset with ML projections
    result_data = data.copy()
    for idx, row in unknown_samples.iterrows():
        result_data.loc[result_data['IID'] == row['IID'], 'KNN_prediction'] = row['KNN_prediction']
        result_data.loc[result_data['IID'] == row['IID'], 'RandomForest_prediction'] = row['RandomForest_prediction']
        result_data.loc[result_data['IID'] == row['IID'], 'ml_consensus_prediction'] = row['ml_consensus_prediction']
        result_data.loc[result_data['IID'] == row['IID'], 'KNN_confidence'] = row['KNN_confidence']
        result_data.loc[result_data['IID'] == row['IID'], 'RandomForest_confidence'] = row['RandomForest_confidence']
        result_data.loc[result_data['IID'] == row['IID'], 'ml_consensus_confidence'] = row['ml_consensus_confidence']
    
    # Add projection status
    result_data['projection_status'] = 'Stage1_Assigned'
    result_data.loc[result_data['dataset'] == 'MAGE', 'projection_status'] = 'Reference'
    result_data.loc[result_data['inferred_ancestry'] == 'Unknown', 'projection_status'] = 'Stage2_ML_Projected'
    
    print(f"Projected {len(unknown_samples)} unknown samples from Stage 1")
    print(f"Stage 2 ML consensus predictions:")
    if 'ml_consensus_prediction' in unknown_samples.columns:
        for pred, count in unknown_samples['ml_consensus_prediction'].value_counts().items():
            print(f"  {pred}: {count} samples")
    
    return result_data

def create_projection_plots(data, reference_samples):
    """Create comprehensive projection visualization"""
    print("Creating projection plots...")
    
    # Set up plotting
    plt.style.use('default')
    sns.set_palette("husl")
    
    # Create figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Reference Population Projection Analysis\nMAGE Samples as Reference Populations', 
                 fontsize=16, fontweight='bold')
    
    # Define colors
    ancestry_colors = {
        'EUR': '#1f77b4',
        'AFR': '#ff7f0e',
        'EAS': '#2ca02c', 
        'SAS': '#d62728',
        'AMR': '#9467bd',
        'ADMIXED': '#8c564b',
        'Other': '#e377c2',
        'Uncertain': '#7f7f7f'
    }
    
    # Plot 1: Reference populations in PC space
    ax1 = axes[0, 0]
    for ancestry in reference_samples['ancestry_standard'].unique():
        mask = reference_samples['ancestry_standard'] == ancestry
        if mask.sum() > 0:
            ax1.scatter(reference_samples.loc[mask, 'PC1'], 
                       reference_samples.loc[mask, 'PC2'],
                       c=ancestry_colors.get(ancestry, '#7f7f7f'),
                       label=f'{ancestry} (n={mask.sum()})', 
                       alpha=0.8, s=40, edgecolor='black', linewidth=0.5)
    
    ax1.set_xlabel('PC1 (45.26%)')
    ax1.set_ylabel('PC2 (13.79%)')
    ax1.set_title('Reference Populations (MAGE Samples)')
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: All samples with projection status
    ax2 = axes[0, 1]
    
    # Plot reference samples
    ref_mask = data['projection_status'] == 'Reference'
    if ref_mask.sum() > 0:
        ax2.scatter(data.loc[ref_mask, 'PC1'], data.loc[ref_mask, 'PC2'],
                   c='blue', alpha=0.6, s=20, label=f'Reference (n={ref_mask.sum()})')
    
    # Plot Stage 1 assigned samples
    stage1_mask = data['projection_status'] == 'Stage1_Assigned'
    if stage1_mask.sum() > 0:
        ax2.scatter(data.loc[stage1_mask, 'PC1'], data.loc[stage1_mask, 'PC2'],
                   c='green', alpha=0.6, s=20, label=f'Stage1 Assigned (n={stage1_mask.sum()})')
    
    # Plot projected samples (Stage 2 ML)
    proj_mask = data['projection_status'] == 'Stage2_ML_Projected'
    if proj_mask.sum() > 0:
        ax2.scatter(data.loc[proj_mask, 'PC1'], data.loc[proj_mask, 'PC2'],
                   c='red', alpha=0.6, s=20, label=f'Stage2 ML Projected (n={proj_mask.sum()})')
    
    ax2.set_xlabel('PC1 (45.26%)')
    ax2.set_ylabel('PC2 (13.79%)')
    ax2.set_title('Sample Assignment Method')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: ML Consensus predictions (only projected samples)
    ax3 = axes[1, 0]
    projected_data = data[data['projection_status'] == 'Stage2_ML_Projected']
    
    if len(projected_data) > 0 and 'ml_consensus_prediction' in projected_data.columns:
        for pred in projected_data['ml_consensus_prediction'].unique():
            if pd.notna(pred):
                mask = projected_data['ml_consensus_prediction'] == pred
                if mask.sum() > 0:
                    ax3.scatter(projected_data.loc[mask, 'PC1'], projected_data.loc[mask, 'PC2'],
                               c=ancestry_colors.get(pred, '#7f7f7f'),
                               label=f'{pred} (n={mask.sum()})', 
                               alpha=0.8, s=30)
    
    # Also show reference samples for context (smaller, lighter)
    for ancestry in reference_samples['ancestry_standard'].unique():
        mask = reference_samples['ancestry_standard'] == ancestry
        if mask.sum() > 0:
            ax3.scatter(reference_samples.loc[mask, 'PC1'], 
                       reference_samples.loc[mask, 'PC2'],
                       c=ancestry_colors.get(ancestry, '#7f7f7f'),
                       alpha=0.3, s=10, edgecolor='none')
    
    ax3.set_xlabel('PC1 (45.26%)')
    ax3.set_ylabel('PC2 (13.79%)')
    ax3.set_title('Stage 2 ML Consensus Predictions')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Confidence distribution
    ax4 = axes[1, 1]
    projected_data = data[data['projection_status'] == 'Stage2_ML_Projected']
    
    if len(projected_data) > 0 and 'KNN_confidence' in projected_data.columns:
        # Filter out NaN values
        knn_conf = projected_data['KNN_confidence'].dropna()
        rf_conf = projected_data['RandomForest_confidence'].dropna()
        
        if len(knn_conf) > 0:
            ax4.hist(knn_conf, bins=15, alpha=0.7, 
                    label=f'KNN Confidence (n={len(knn_conf)})', color='skyblue', density=True)
        
        if len(rf_conf) > 0:
            ax4.hist(rf_conf, bins=15, alpha=0.7, 
                    label=f'RF Confidence (n={len(rf_conf)})', color='lightcoral', density=True)
        
        ax4.axvline(x=0.7, color='red', linestyle='--', linewidth=2, label='High Confidence Threshold')
        ax4.axvline(x=0.5, color='orange', linestyle='--', linewidth=2, label='Medium Confidence Threshold')
        
        # Add statistics text
        if len(knn_conf) > 0 and len(rf_conf) > 0:
            stats_text = f'Mean KNN: {knn_conf.mean():.3f}\nMean RF: {rf_conf.mean():.3f}'
            ax4.text(0.02, 0.98, stats_text, transform=ax4.transAxes, 
                    verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    else:
        ax4.text(0.5, 0.5, 'No projected samples\nwith confidence scores', 
                ha='center', va='center', transform=ax4.transAxes, fontsize=12)
    
    ax4.set_xlabel('Prediction Confidence')
    ax4.set_ylabel('Density')
    ax4.set_title('Stage 2 ML Confidence Distribution')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save the plot
    plot_file = results_dir / "reference_projection_analysis.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"Saved projection plot to {plot_file}")
    
    plt.show()

def generate_projection_report(data, reference_samples):
    """Generate comprehensive projection analysis report"""
    print("Generating projection analysis report...")
    
    report = []
    report.append("# Reference Population Projection Analysis Report")
    report.append("=" * 60)
    report.append("")
    
    report.append("## Overview")
    report.append("Analysis projects ADNI samples onto reference population PC space")
    report.append("using MAGE samples with known 1000 Genomes ancestry as references.")
    report.append("")
    
    report.append("## Reference Population Composition")
    ref_counts = reference_samples['ancestry_standard'].value_counts()
    for ancestry, count in ref_counts.items():
        pct = count / len(reference_samples) * 100
        report.append(f"- {ancestry}: {count} samples ({pct:.1f}%)")
    report.append("")
    
    report.append("## Projection Results")
    total_samples = len(data)
    reference_count = (data['projection_status'] == 'Reference').sum()
    projected_count = (data['projection_status'] == 'Projected').sum()
    
    report.append(f"- Total samples: {total_samples:,}")
    report.append(f"- Reference samples: {reference_count:,} ({reference_count/total_samples*100:.1f}%)")
    report.append(f"- Projected samples: {projected_count:,} ({projected_count/total_samples*100:.1f}%)")
    report.append("")
    
    if 'consensus_prediction' in data.columns:
        report.append("## Consensus Ancestry Predictions")
        projected_data = data[data['projection_status'] == 'Projected']
        if len(projected_data) > 0:
            pred_counts = projected_data['consensus_prediction'].value_counts()
            for pred, count in pred_counts.items():
                pct = count / len(projected_data) * 100
                report.append(f"- {pred}: {count} samples ({pct:.1f}%)")
        report.append("")
    
    report.append("## Confidence Metrics")
    if 'KNN_confidence' in data.columns:
        projected_data = data[data['projection_status'] == 'Projected']
        if len(projected_data) > 0:
            high_conf_knn = (projected_data['KNN_confidence'] > 0.7).sum()
            high_conf_rf = (projected_data['RandomForest_confidence'] > 0.7).sum()
            report.append(f"- High confidence KNN predictions (>70%): {high_conf_knn} samples")
            report.append(f"- High confidence RF predictions (>70%): {high_conf_rf} samples")
            report.append(f"- Mean KNN confidence: {projected_data['KNN_confidence'].mean():.3f}")
            report.append(f"- Mean RF confidence: {projected_data['RandomForest_confidence'].mean():.3f}")
    report.append("")
    
    report.append("## Methodology")
    report.append("- Reference populations: MAGE samples with 1000 Genomes ancestry")
    report.append("- Projection method: K-Nearest Neighbors + Random Forest ensemble")
    report.append("- Feature space: Top 20 principal components")
    report.append("- Confidence threshold: 70% for high-confidence predictions")
    report.append("- Consensus rule: Use KNN if confidence >70%, else RF if confidence >70%")
    
    # Save report
    report_file = results_dir / "reference_projection_report.txt"
    with open(report_file, 'w') as f:
        f.write('\n'.join(report))
    
    print(f"Saved projection report to {report_file}")
    
    return report

def save_projection_results(data):
    """Save projection results to CSV files"""
    print("Saving projection results...")
    
    # Save complete results
    results_file = results_dir / "reference_projection_results.csv"
    data.to_csv(results_file, index=False)
    print(f"Saved complete results to {results_file}")
    
    # Save just the projected samples
    projected_samples = data[data['projection_status'] == 'Projected']
    if len(projected_samples) > 0:
        projected_file = results_dir / "projected_ancestry_assignments.csv"
        projected_samples[['IID', 'sample_type', 'KNN_prediction', 'RandomForest_prediction', 
                          'consensus_prediction', 'KNN_confidence', 'RandomForest_confidence']].to_csv(
            projected_file, index=False)
        print(f"Saved projected ancestry assignments to {projected_file}")
    
    # Save ancestry summary
    summary_data = []
    for sample_type in data['sample_type'].unique():
        subset = data[data['sample_type'] == sample_type]
        if 'consensus_prediction' in subset.columns:
            for ancestry in subset['consensus_prediction'].unique():
                count = (subset['consensus_prediction'] == ancestry).sum()
                summary_data.append({
                    'sample_type': sample_type,
                    'ancestry': ancestry,
                    'count': count
                })
    
    summary_df = pd.DataFrame(summary_data)
    summary_file = results_dir / "ancestry_by_sample_type.csv"
    summary_df.to_csv(summary_file, index=False)
    print(f"Saved ancestry summary to {summary_file}")

def main():
    """Main projection analysis function"""
    print("Starting Reference Population Projection Analysis...")
    print("=" * 60)
    
    # Load data
    data = load_data()
    
    # Create reference populations
    reference_samples = create_reference_populations(data)
    
    if len(reference_samples) == 0:
        print("ERROR: No reference samples found!")
        return
    
    # Train classifiers
    classifiers, pc_cols = train_ancestry_classifiers(reference_samples)
    
    # Project ADNI samples
    projected_data = project_unknown_samples(data, reference_samples, classifiers, pc_cols)
    
    # Create visualizations
    create_projection_plots(projected_data, reference_samples)
    
    # Generate report
    generate_projection_report(projected_data, reference_samples)
    
    # Save results
    save_projection_results(projected_data)
    
    print("\nReference population projection analysis complete! 🎉")
    print(f"Results saved to: {results_dir}")

if __name__ == "__main__":
    main()