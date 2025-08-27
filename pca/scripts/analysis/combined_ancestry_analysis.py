#!/usr/bin/env python3
"""
Combined Ancestry Analysis
Creates comprehensive comparison of all ancestry sources:
- Self-reported ethnicity
- MAGE 1000 Genomes reference populations  
- Stage 1 reference-guided PCA clustering
- Stage 2 ML inference
- Final hybrid consensus
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Set up file paths
base_dir = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/wgs")
pca_dir = base_dir / "pca"
results_dir = pca_dir / "results" / "combined_analysis"
results_dir.mkdir(parents=True, exist_ok=True)

def load_all_data():
    """Load and combine all ancestry data sources"""
    print("Loading comprehensive ancestry data...")
    
    # Load Stage 1 PCA results
    pca_file = pca_dir / "results" / "plink_analysis" / "plink_pca_results.csv"
    pca_data = pd.read_csv(pca_file)
    
    # Load Stage 2 ML results  
    ml_file = pca_dir / "results" / "reference_projection" / "reference_projection_results.csv"
    ml_data = pd.read_csv(ml_file)
    
    # Load donor metadata
    donor_file = pca_dir / "results" / "updated_metadata" / "donor_summary.csv"
    donor_data = pd.read_csv(donor_file)
    
    # Start with PCA data as base
    combined_data = pca_data.copy()
    
    # Add donor metadata (avoiding PC column conflicts)
    metadata_cols = ['donor_id', 'dataset', 'self_reported_ethnicity', 'ethnicity_ontology_term', 
                    'mage_1kgp_population', 'mage_1kgp_ancestry', 'genomic_ancestry_pca', 
                    'genomic_ancestry_ml', 'ml_confidence', 'ethnicity_ancestry_concordance']
    
    donor_metadata = donor_data[metadata_cols].copy()
    combined_data = combined_data.merge(donor_metadata, left_on='IID', right_on='donor_id', how='left')
    
    # Add Stage 2 ML results for samples that were projected
    ml_cols = ['IID', 'KNN_prediction', 'RandomForest_prediction', 'ml_consensus_prediction', 
               'KNN_confidence', 'RandomForest_confidence', 'ml_consensus_confidence', 'projection_status']
    
    ml_subset = ml_data[ml_cols].copy()
    combined_data = combined_data.merge(ml_subset, on='IID', how='left')
    
    print(f"Combined dataset: {len(combined_data)} samples")
    return combined_data

def create_final_consensus_ancestry(data):
    """Create final hybrid consensus ancestry assignments"""
    print("Creating final hybrid consensus ancestry...")
    
    # Create final consensus column
    data['final_consensus_ancestry'] = data['inferred_ancestry'].copy()
    data['final_consensus_confidence'] = data['assignment_confidence'].copy()
    data['final_consensus_method'] = 'Stage1_PCA'
    
    # For samples that were uncertain in Stage 1 but classified by Stage 2 ML
    ml_mask = (data['inferred_ancestry'] == 'Unknown') & (data['ml_consensus_prediction'].notna())
    
    # Use ML prediction if confidence > 0.7, otherwise keep as uncertain
    high_conf_ml = ml_mask & (data['ml_consensus_confidence'] > 0.7)
    med_conf_ml = ml_mask & (data['ml_consensus_confidence'] > 0.5) & (data['ml_consensus_confidence'] <= 0.7)
    
    # High confidence ML predictions
    data.loc[high_conf_ml, 'final_consensus_ancestry'] = data.loc[high_conf_ml, 'ml_consensus_prediction']
    data.loc[high_conf_ml, 'final_consensus_confidence'] = data.loc[high_conf_ml, 'ml_consensus_confidence']
    data.loc[high_conf_ml, 'final_consensus_method'] = 'Stage2_ML_High'
    
    # Medium confidence ML predictions (lower priority)
    data.loc[med_conf_ml, 'final_consensus_ancestry'] = data.loc[med_conf_ml, 'ml_consensus_prediction']
    data.loc[med_conf_ml, 'final_consensus_confidence'] = data.loc[med_conf_ml, 'ml_consensus_confidence']
    data.loc[med_conf_ml, 'final_consensus_method'] = 'Stage2_ML_Med'
    
    # Still uncertain after both stages
    still_uncertain = (data['final_consensus_ancestry'] == 'Unknown') | (data['final_consensus_ancestry'] == 'Uncertain')
    data.loc[still_uncertain, 'final_consensus_method'] = 'Still_Uncertain'
    
    print("Final consensus ancestry distribution:")
    final_counts = data['final_consensus_ancestry'].value_counts()
    for ancestry, count in final_counts.items():
        pct = count / len(data) * 100
        print(f"  {ancestry}: {count:,} ({pct:.1f}%)")
    
    print("\nFinal consensus method distribution:")
    method_counts = data['final_consensus_method'].value_counts()
    for method, count in method_counts.items():
        pct = count / len(data) * 100
        print(f"  {method}: {count:,} ({pct:.1f}%)")
    
    return data

def create_concordance_analysis(data):
    """Analyze concordance between different ancestry sources"""
    print("Performing concordance analysis...")
    
    # Focus on samples with both self-reported ethnicity and genomic ancestry
    analysis_data = data[
        (data['self_reported_ethnicity'].notna()) & 
        (data['final_consensus_ancestry'] != 'Unknown') &
        (data['final_consensus_ancestry'] != 'Uncertain')
    ].copy()
    
    if len(analysis_data) == 0:
        print("No samples with both self-reported and inferred ancestry for concordance analysis")
        return pd.DataFrame()
    
    # Standardize self-reported ethnicity labels
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
        'other': 'OTHER'
    }
    
    analysis_data['self_reported_ancestry'] = analysis_data['self_reported_ethnicity'].str.lower().map(ethnicity_mapping)
    
    # Calculate concordance rates
    concordance_results = []
    
    for dataset in analysis_data['dataset'].unique():
        if pd.isna(dataset):
            continue
            
        subset = analysis_data[analysis_data['dataset'] == dataset]
        if len(subset) == 0:
            continue
            
        # Overall concordance
        concordant = subset['self_reported_ancestry'] == subset['final_consensus_ancestry']
        concordance_rate = concordant.sum() / len(subset) * 100
        
        concordance_results.append({
            'dataset': dataset,
            'total_samples': len(subset),
            'concordant_samples': concordant.sum(),
            'concordance_rate': concordance_rate
        })
        
        print(f"{dataset} concordance: {concordant.sum()}/{len(subset)} ({concordance_rate:.1f}%)")
    
    concordance_df = pd.DataFrame(concordance_results)
    return concordance_df, analysis_data

def create_comprehensive_visualization(data):
    """Create comprehensive ancestry comparison visualization"""
    print("Creating comprehensive visualization...")
    
    # Set up the plotting style
    plt.style.use('default')
    sns.set_palette("husl")
    
    # Create figure with subplots
    fig, axes = plt.subplots(2, 3, figsize=(20, 12))
    fig.suptitle('Comprehensive Ancestry Inference Analysis\nComparison of All Methods and Sources', 
                 fontsize=16, fontweight='bold')
    
    # Define colors
    ancestry_colors = {
        'EUR': '#1f77b4',
        'AFR': '#ff7f0e',
        'EAS': '#2ca02c', 
        'SAS': '#d62728',
        'AMR': '#9467bd',
        'Unknown': '#7f7f7f',
        'Uncertain': '#7f7f7f',
        'MIXED': '#8c564b',
        'OTHER': '#e377c2'
    }
    
    # Plot 1: Stage 1 PCA Results
    ax1 = axes[0, 0]
    for ancestry in data['inferred_ancestry'].unique():
        mask = data['inferred_ancestry'] == ancestry
        if mask.sum() > 0:
            ax1.scatter(data.loc[mask, 'PC1'], data.loc[mask, 'PC2'],
                       c=ancestry_colors.get(ancestry, '#7f7f7f'),
                       label=f'{ancestry} (n={mask.sum()})', alpha=0.7, s=15)
    
    ax1.set_xlabel('PC1 (45.26%)')
    ax1.set_ylabel('PC2 (13.79%)')
    ax1.set_title('Stage 1: Reference-Guided PCA')
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Final Consensus Results
    ax2 = axes[0, 1]
    for ancestry in data['final_consensus_ancestry'].unique():
        if pd.notna(ancestry):
            mask = data['final_consensus_ancestry'] == ancestry
            if mask.sum() > 0:
                ax2.scatter(data.loc[mask, 'PC1'], data.loc[mask, 'PC2'],
                           c=ancestry_colors.get(ancestry, '#7f7f7f'),
                           label=f'{ancestry} (n={mask.sum()})', alpha=0.7, s=15)
    
    ax2.set_xlabel('PC1 (45.26%)')
    ax2.set_ylabel('PC2 (13.79%)')
    ax2.set_title('Final Hybrid Consensus')
    ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Sample Type Distribution
    ax3 = axes[0, 2]
    sample_counts = data['dataset'].value_counts()
    ax3.pie(sample_counts.values, labels=sample_counts.index, autopct='%1.1f%%', startangle=90)
    ax3.set_title('Sample Composition by Dataset')
    
    # Plot 4: Method Comparison
    ax4 = axes[1, 0]
    method_counts = data['final_consensus_method'].value_counts()
    bars = ax4.bar(range(len(method_counts)), method_counts.values, color='skyblue')
    ax4.set_xticks(range(len(method_counts)))
    ax4.set_xticklabels(method_counts.index, rotation=45, ha='right')
    ax4.set_ylabel('Number of Samples')
    ax4.set_title('Final Assignment Method')
    ax4.grid(True, alpha=0.3)
    
    # Add value labels on bars
    for bar, count in zip(bars, method_counts.values):
        ax4.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 10,
                f'{count}', ha='center', va='bottom')
    
    # Plot 5: Confidence Distribution
    ax5 = axes[1, 1] 
    assigned_samples = data[data['final_consensus_ancestry'] != 'Unknown']
    if len(assigned_samples) > 0:
        ax5.hist(assigned_samples['final_consensus_confidence'], bins=20, alpha=0.7, 
                color='lightcoral', edgecolor='black')
        ax5.axvline(x=0.7, color='red', linestyle='--', label='High Confidence Threshold')
        ax5.axvline(x=0.5, color='orange', linestyle='--', label='Medium Confidence Threshold')
    
    ax5.set_xlabel('Confidence Score')
    ax5.set_ylabel('Number of Samples')
    ax5.set_title('Final Assignment Confidence Distribution')
    ax5.legend()
    ax5.grid(True, alpha=0.3)
    
    # Plot 6: Performance Comparison
    ax6 = axes[1, 2]
    
    # Calculate assignment rates
    stage1_assigned = (data['inferred_ancestry'] != 'Unknown').sum()
    stage1_rate = stage1_assigned / len(data) * 100
    
    final_assigned = (data['final_consensus_ancestry'] != 'Unknown').sum()
    final_rate = final_assigned / len(data) * 100
    
    categories = ['Stage 1\nPCA Only', 'Final\nHybrid']
    rates = [stage1_rate, final_rate]
    
    bars = ax6.bar(categories, rates, color=['lightblue', 'lightgreen'])
    ax6.set_ylabel('Assignment Rate (%)')
    ax6.set_title('Pipeline Performance Comparison')
    ax6.set_ylim(0, 100)
    ax6.grid(True, alpha=0.3)
    
    # Add value labels on bars
    for bar, rate in zip(bars, rates):
        ax6.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                f'{rate:.1f}%', ha='center', va='bottom', fontweight='bold')
    
    plt.tight_layout()
    
    # Save the plot
    plot_file = results_dir / "comprehensive_ancestry_analysis.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"Saved comprehensive plot to {plot_file}")
    
    plt.show()

def generate_comprehensive_report(data, concordance_df):
    """Generate comprehensive analysis report"""
    print("Generating comprehensive analysis report...")
    
    report = []
    report.append("# Comprehensive Ancestry Inference Analysis Report")
    report.append("=" * 65)
    report.append("")
    
    report.append("## Executive Summary")
    report.append("This report presents the complete results of the hybrid ancestry inference")
    report.append("pipeline, comparing all methods and data sources.")
    report.append("")
    
    report.append("## Dataset Overview")
    report.append(f"- Total samples: {len(data):,}")
    dataset_counts = data['dataset'].value_counts()
    for dataset, count in dataset_counts.items():
        if pd.notna(dataset):
            pct = count / len(data) * 100
            report.append(f"- {dataset}: {count:,} samples ({pct:.1f}%)")
    report.append("")
    
    report.append("## Performance Summary")
    
    # Stage 1 performance
    stage1_assigned = (data['inferred_ancestry'] != 'Unknown').sum()
    stage1_rate = stage1_assigned / len(data) * 100
    report.append(f"### Stage 1 (Reference-Guided PCA)")
    report.append(f"- Assignment rate: {stage1_rate:.1f}% ({stage1_assigned:,}/{len(data):,} samples)")
    
    stage1_counts = data['inferred_ancestry'].value_counts()
    for ancestry, count in stage1_counts.items():
        pct = count / len(data) * 100
        report.append(f"- {ancestry}: {count:,} ({pct:.1f}%)")
    report.append("")
    
    # Stage 2 performance  
    ml_projected = (data['projection_status'] == 'Stage2_ML_Projected').sum()
    if ml_projected > 0:
        report.append(f"### Stage 2 (ML Projection)")
        report.append(f"- Samples projected: {ml_projected:,}")
        
        stage2_data = data[data['projection_status'] == 'Stage2_ML_Projected']
        ml_counts = stage2_data['ml_consensus_prediction'].value_counts()
        for pred, count in ml_counts.items():
            if pd.notna(pred):
                pct = count / ml_projected * 100
                report.append(f"- {pred}: {count:,} ({pct:.1f}%)")
        report.append("")
    
    # Final consensus performance
    final_assigned = (data['final_consensus_ancestry'] != 'Unknown').sum()
    final_rate = final_assigned / len(data) * 100
    report.append(f"### Final Hybrid Consensus")
    report.append(f"- Overall assignment rate: {final_rate:.1f}% ({final_assigned:,}/{len(data):,} samples)")
    
    final_counts = data['final_consensus_ancestry'].value_counts()
    for ancestry, count in final_counts.items():
        pct = count / len(data) * 100
        report.append(f"- {ancestry}: {count:,} ({pct:.1f}%)")
    report.append("")
    
    # Method breakdown
    report.append("### Assignment Method Breakdown")
    method_counts = data['final_consensus_method'].value_counts()
    for method, count in method_counts.items():
        pct = count / len(data) * 100
        report.append(f"- {method}: {count:,} ({pct:.1f}%)")
    report.append("")
    
    # Concordance analysis
    if len(concordance_df) > 0:
        report.append("## Concordance Analysis")
        report.append("Agreement between self-reported ethnicity and genomic ancestry:")
        for _, row in concordance_df.iterrows():
            report.append(f"- {row['dataset']}: {row['concordance_rate']:.1f}% ({row['concordant_samples']}/{row['total_samples']})")
        report.append("")
    
    # Confidence metrics
    assigned_samples = data[(data['final_consensus_ancestry'] != 'Unknown') & 
                           (data['final_consensus_confidence'].notna())]
    if len(assigned_samples) > 0:
        report.append("## Confidence Metrics")
        mean_conf = assigned_samples['final_consensus_confidence'].mean()
        high_conf = (assigned_samples['final_consensus_confidence'] > 0.7).sum()
        med_conf = ((assigned_samples['final_consensus_confidence'] >= 0.5) & 
                   (assigned_samples['final_consensus_confidence'] <= 0.7)).sum()
        low_conf = (assigned_samples['final_consensus_confidence'] < 0.5).sum()
        
        report.append(f"- Mean confidence: {mean_conf:.3f}")
        report.append(f"- High confidence (>0.7): {high_conf:,} ({high_conf/len(assigned_samples)*100:.1f}%)")
        report.append(f"- Medium confidence (0.5-0.7): {med_conf:,} ({med_conf/len(assigned_samples)*100:.1f}%)")
        report.append(f"- Low confidence (<0.5): {low_conf:,} ({low_conf/len(assigned_samples)*100:.1f}%)")
        report.append("")
    
    report.append("## Methodology")
    report.append("### Stage 1: Reference-Guided PCA Clustering")
    report.append("- Uses MAGE 1000 Genomes population centroids")
    report.append("- Assignment threshold: 2.5 standard deviations")
    report.append("- Minimum confidence: 0.3")
    report.append("")
    
    report.append("### Stage 2: Machine Learning Projection")
    report.append("- Ensemble: K-Nearest Neighbors + Random Forest")
    report.append("- Feature space: 20 principal components")
    report.append("- Confidence threshold: 0.7 for high confidence")
    report.append("")
    
    report.append("### Final Consensus")
    report.append("- Prioritizes Stage 1 assignments")
    report.append("- Uses Stage 2 for high-confidence ML predictions (>0.7)")
    report.append("- Maintains uncertainty for low-confidence predictions")
    
    # Save report
    report_file = results_dir / "comprehensive_ancestry_analysis_report.txt"
    with open(report_file, 'w') as f:
        f.write('\n'.join(report))
    
    print(f"Saved comprehensive report to {report_file}")
    
    return report

def save_combined_results(data):
    """Save all combined results"""
    print("Saving combined ancestry results...")
    
    # Save complete dataset
    complete_file = results_dir / "complete_ancestry_results.csv"
    data.to_csv(complete_file, index=False)
    print(f"Saved complete results to {complete_file}")
    
    # Save key comparisons
    comparison_cols = ['IID', 'dataset', 'sample_type', 'self_reported_ethnicity', 
                      'mage_1kgp_population', 'mage_1kgp_ancestry',
                      'inferred_ancestry', 'assignment_confidence',
                      'ml_consensus_prediction', 'ml_consensus_confidence',
                      'final_consensus_ancestry', 'final_consensus_confidence', 'final_consensus_method',
                      'PC1', 'PC2']
    
    comparison_data = data[comparison_cols].copy()
    comparison_file = results_dir / "ancestry_method_comparison.csv"
    comparison_data.to_csv(comparison_file, index=False)
    print(f"Saved method comparison to {comparison_file}")
    
    # Save summary statistics
    summary_stats = []
    
    # Overall statistics
    summary_stats.append({
        'metric': 'total_samples',
        'value': len(data),
        'percentage': 100.0
    })
    
    # Stage 1 statistics
    stage1_assigned = (data['inferred_ancestry'] != 'Unknown').sum()
    summary_stats.append({
        'metric': 'stage1_assigned',
        'value': stage1_assigned,
        'percentage': stage1_assigned / len(data) * 100
    })
    
    # Final consensus statistics
    final_assigned = (data['final_consensus_ancestry'] != 'Unknown').sum()
    summary_stats.append({
        'metric': 'final_assigned', 
        'value': final_assigned,
        'percentage': final_assigned / len(data) * 100
    })
    
    # Improvement
    improvement = final_assigned - stage1_assigned
    summary_stats.append({
        'metric': 'stage2_improvement',
        'value': improvement,
        'percentage': improvement / len(data) * 100
    })
    
    summary_df = pd.DataFrame(summary_stats)
    summary_file = results_dir / "pipeline_performance_summary.csv"
    summary_df.to_csv(summary_file, index=False)
    print(f"Saved performance summary to {summary_file}")

def main():
    """Main comprehensive analysis function"""
    print("Starting Comprehensive Ancestry Analysis...")
    print("=" * 60)
    
    # Load all data
    data = load_all_data()
    
    # Create final consensus
    data = create_final_consensus_ancestry(data)
    
    # Concordance analysis
    concordance_df, analysis_data = create_concordance_analysis(data)
    
    # Create visualizations
    create_comprehensive_visualization(data)
    
    # Generate report
    generate_comprehensive_report(data, concordance_df)
    
    # Save results
    save_combined_results(data)
    
    print("\nComprehensive ancestry analysis complete! 🎉")
    print(f"Results saved to: {results_dir}")

if __name__ == "__main__":
    main()