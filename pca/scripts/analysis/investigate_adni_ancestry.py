#!/usr/bin/env python3
"""
Investigate ADNI Ancestry Discordance
Deep dive into why ADNI samples show low concordance between self-reported and genomic ancestry
Testing alternative hypotheses:
1. AMR individuals self-identifying as "white"
2. Errors in ADNI self-reported ethnicity data
3. Population structure differences
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Set up file paths
base_dir = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/wgs")
pca_dir = base_dir / "pca"
results_dir = pca_dir / "results" / "adni_investigation"
results_dir.mkdir(parents=True, exist_ok=True)

def investigate_adni_ancestry():
    """Investigate ADNI ancestry assignment in detail"""
    print("Investigating ADNI Ancestry Assignments...")
    print("=" * 60)
    
    # Load complete results
    complete_file = pca_dir / "results" / "combined_analysis" / "complete_ancestry_results.csv"
    data = pd.read_csv(complete_file)
    
    # Focus on ADNI samples with self-reported ethnicity
    adni_data = data[
        (data['dataset'] == 'ADNI') & 
        (data['self_reported_ethnicity'].notna())
    ].copy()
    
    print(f"ADNI samples with self-reported ethnicity: {len(adni_data)}")
    
    # Check for missing genomic ancestry assignments
    adni_assigned = adni_data[
        (adni_data['final_consensus_ancestry'] != 'Unknown') &
        (adni_data['final_consensus_ancestry'] != 'Uncertain')
    ]
    
    print(f"ADNI samples with genomic ancestry assigned: {len(adni_assigned)}")
    
    # Hypothesis 1: AMR individuals self-identifying as "white"
    print("\n" + "="*60)
    print("HYPOTHESIS 1: AMR INDIVIDUALS SELF-IDENTIFYING AS 'WHITE'")
    print("="*60)
    
    white_identified = adni_assigned[adni_assigned['self_reported_ethnicity'] == 'white']
    print(f"ADNI samples self-identifying as 'white': {len(white_identified)}")
    
    white_genomic_breakdown = white_identified['final_consensus_ancestry'].value_counts()
    print("\nGenomic ancestry of self-identified 'white' ADNI participants:")
    for ancestry, count in white_genomic_breakdown.items():
        pct = count / len(white_identified) * 100
        print(f"  {ancestry}: {count} ({pct:.1f}%)")
    
    # Check PC coordinates of white-identified samples
    print("\n--- PC Coordinate Analysis of 'White' ADNI Samples ---")
    white_amr = white_identified[white_identified['final_consensus_ancestry'] == 'AMR']
    white_eur = white_identified[white_identified['final_consensus_ancestry'] == 'EUR']
    white_afr = white_identified[white_identified['final_consensus_ancestry'] == 'AFR']
    
    if len(white_amr) > 0:
        print(f"\n'White' → AMR samples (n={len(white_amr)}):")
        print(f"  PC1 range: {white_amr['PC1'].min():.4f} to {white_amr['PC1'].max():.4f}")
        print(f"  PC2 range: {white_amr['PC2'].min():.4f} to {white_amr['PC2'].max():.4f}")
        print(f"  PC1 mean: {white_amr['PC1'].mean():.4f} ± {white_amr['PC1'].std():.4f}")
        print(f"  PC2 mean: {white_amr['PC2'].mean():.4f} ± {white_amr['PC2'].std():.4f}")
    
    if len(white_eur) > 0:
        print(f"\n'White' → EUR samples (n={len(white_eur)}):")
        print(f"  PC1 range: {white_eur['PC1'].min():.4f} to {white_eur['PC1'].max():.4f}")
        print(f"  PC2 range: {white_eur['PC2'].min():.4f} to {white_eur['PC2'].max():.4f}")
        print(f"  PC1 mean: {white_eur['PC1'].mean():.4f} ± {white_eur['PC1'].std():.4f}")
        print(f"  PC2 mean: {white_eur['PC2'].mean():.4f} ± {white_eur['PC2'].std():.4f}")
    
    if len(white_afr) > 0:
        print(f"\n'White' → AFR samples (n={len(white_afr)}):")
        print(f"  PC1 range: {white_afr['PC1'].min():.4f} to {white_afr['PC1'].max():.4f}")
        print(f"  PC2 range: {white_afr['PC2'].min():.4f} to {white_afr['PC2'].max():.4f}")
        print(f"  PC1 mean: {white_afr['PC1'].mean():.4f} ± {white_afr['PC1'].std():.4f}")
        print(f"  PC2 mean: {white_afr['PC2'].mean():.4f} ± {white_afr['PC2'].std():.4f}")
    
    # Compare with reference populations
    print("\n--- Reference Population Comparison ---")
    
    # Load reference populations from MAGE
    mage_data = data[data['dataset'] == 'MAGE'].copy()
    
    if len(mage_data) > 0:
        print("\nMAGE Reference Population Centers:")
        for ancestry in ['EUR', 'AMR', 'AFR', 'EAS', 'SAS']:
            ref_samples = mage_data[mage_data['mage_1kgp_ancestry'] == ancestry]
            if len(ref_samples) > 0:
                print(f"  {ancestry} (n={len(ref_samples)}): PC1={ref_samples['PC1'].mean():.4f}, PC2={ref_samples['PC2'].mean():.4f}")
    
    # Hypothesis 2: Check for systematic biases in ADNI data
    print("\n" + "="*60)
    print("HYPOTHESIS 2: SYSTEMATIC ISSUES IN ADNI DATA")
    print("="*60)
    
    # Check confidence scores
    print("\nConfidence Score Analysis:")
    print("ADNI genomic ancestry assignment confidence:")
    
    adni_confidence = adni_assigned['final_consensus_confidence'].dropna()
    if len(adni_confidence) > 0:
        print(f"  Mean confidence: {adni_confidence.mean():.3f}")
        print(f"  Median confidence: {adni_confidence.median():.3f}")
        print(f"  High confidence (>0.7): {(adni_confidence > 0.7).sum()} ({(adni_confidence > 0.7).sum()/len(adni_confidence)*100:.1f}%)")
        print(f"  Low confidence (<0.5): {(adni_confidence < 0.5).sum()} ({(adni_confidence < 0.5).sum()/len(adni_confidence)*100:.1f}%)")
    
    # Compare confidence by assignment method
    print("\nAssignment Method Breakdown for ADNI:")
    method_counts = adni_assigned['final_consensus_method'].value_counts()
    for method, count in method_counts.items():
        pct = count / len(adni_assigned) * 100
        subset_conf = adni_assigned[adni_assigned['final_consensus_method'] == method]['final_consensus_confidence'].mean()
        print(f"  {method}: {count} ({pct:.1f}%) - Mean confidence: {subset_conf:.3f}")
    
    # Check specific discordant cases
    print("\n--- Detailed Discordant Case Analysis ---")
    
    # Focus on the most problematic: white → AMR
    white_to_amr = adni_assigned[
        (adni_assigned['self_reported_ethnicity'] == 'white') &
        (adni_assigned['final_consensus_ancestry'] == 'AMR')
    ]
    
    print(f"\n'White' → AMR cases (n={len(white_to_amr)}):")
    if len(white_to_amr) > 0:
        print("Sample details (first 10):")
        for i, (_, row) in enumerate(white_to_amr.head(10).iterrows()):
            print(f"  {row['IID']}: PC1={row['PC1']:.4f}, PC2={row['PC2']:.4f}, "
                  f"Conf={row['final_consensus_confidence']:.3f}, Method={row['final_consensus_method']}")
    
    # Hypothesis 3: Check for batch effects or population substructure
    print("\n" + "="*60) 
    print("HYPOTHESIS 3: POPULATION SUBSTRUCTURE ANALYSIS")
    print("="*60)
    
    # Compare PC coordinate distributions
    print("\nPC Coordinate Distribution Comparison:")
    
    # ADNI vs GTEx (both should have similar EUR populations)
    gtex_white = data[
        (data['dataset'] == 'GTEx') & 
        (data['self_reported_ethnicity'] == 'white') &
        (data['final_consensus_ancestry'] == 'EUR')
    ]
    
    adni_white = adni_assigned[adni_assigned['self_reported_ethnicity'] == 'white']
    
    if len(gtex_white) > 0 and len(adni_white) > 0:
        print(f"\nGTEx 'white' → EUR samples (n={len(gtex_white)}):")
        print(f"  PC1: {gtex_white['PC1'].mean():.4f} ± {gtex_white['PC1'].std():.4f}")
        print(f"  PC2: {gtex_white['PC2'].mean():.4f} ± {gtex_white['PC2'].std():.4f}")
        
        print(f"\nADNI 'white' samples (all genomic ancestries, n={len(adni_white)}):")
        print(f"  PC1: {adni_white['PC1'].mean():.4f} ± {adni_white['PC1'].std():.4f}")
        print(f"  PC2: {adni_white['PC2'].mean():.4f} ± {adni_white['PC2'].std():.4f}")
        
        # Statistical test
        from scipy import stats
        pc1_ttest = stats.ttest_ind(gtex_white['PC1'], adni_white['PC1'])
        pc2_ttest = stats.ttest_ind(gtex_white['PC2'], adni_white['PC2'])
        
        print(f"\nPC1 difference test: t={pc1_ttest.statistic:.3f}, p={pc1_ttest.pvalue:.2e}")
        print(f"PC2 difference test: t={pc2_ttest.statistic:.3f}, p={pc2_ttest.pvalue:.2e}")
    
    return adni_data, white_to_amr

def create_adni_visualization(adni_data, white_to_amr):
    """Create visualization of ADNI ancestry assignment patterns"""
    print("\nCreating ADNI ancestry visualization...")
    
    # Load all data for context
    complete_file = pca_dir / "results" / "combined_analysis" / "complete_ancestry_results.csv"
    all_data = pd.read_csv(complete_file)
    
    # Create figure
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('ADNI Ancestry Assignment Investigation', fontsize=16, fontweight='bold')
    
    # Define colors
    ancestry_colors = {
        'EUR': '#1f77b4',
        'AFR': '#ff7f0e',
        'EAS': '#2ca02c',
        'SAS': '#d62728',
        'AMR': '#9467bd',
        'Unknown': '#7f7f7f',
        'Uncertain': '#7f7f7f'
    }
    
    # Plot 1: All samples by dataset and ancestry
    ax1 = axes[0, 0]
    
    # Plot MAGE references lightly
    mage_data = all_data[all_data['dataset'] == 'MAGE']
    for ancestry in ['EUR', 'AFR', 'EAS', 'SAS', 'AMR']:
        mage_subset = mage_data[mage_data['mage_1kgp_ancestry'] == ancestry]
        if len(mage_subset) > 0:
            ax1.scatter(mage_subset['PC1'], mage_subset['PC2'],
                       c=ancestry_colors.get(ancestry, '#7f7f7f'),
                       alpha=0.3, s=10, label=f'MAGE {ancestry}')
    
    # Plot ADNI samples prominently
    adni_assigned = adni_data[
        (adni_data['final_consensus_ancestry'] != 'Unknown') &
        (adni_data['final_consensus_ancestry'] != 'Uncertain')
    ]
    
    for ancestry in adni_assigned['final_consensus_ancestry'].unique():
        adni_subset = adni_assigned[adni_assigned['final_consensus_ancestry'] == ancestry]
        if len(adni_subset) > 0:
            ax1.scatter(adni_subset['PC1'], adni_subset['PC2'],
                       c=ancestry_colors.get(ancestry, '#7f7f7f'),
                       alpha=0.8, s=30, edgecolor='black', linewidth=0.5,
                       label=f'ADNI {ancestry} (n={len(adni_subset)})')
    
    ax1.set_xlabel('PC1 (45.26%)')
    ax1.set_ylabel('PC2 (13.79%)')
    ax1.set_title('ADNI vs MAGE Reference Populations')
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    ax1.grid(True, alpha=0.3)
    
    # Plot 2: Focus on discordant cases
    ax2 = axes[0, 1]
    
    # Plot reference populations
    mage_eur = mage_data[mage_data['mage_1kgp_ancestry'] == 'EUR']
    mage_amr = mage_data[mage_data['mage_1kgp_ancestry'] == 'AMR']
    
    if len(mage_eur) > 0:
        ax2.scatter(mage_eur['PC1'], mage_eur['PC2'], c='blue', alpha=0.3, s=15, label='MAGE EUR')
    if len(mage_amr) > 0:
        ax2.scatter(mage_amr['PC1'], mage_amr['PC2'], c='purple', alpha=0.3, s=15, label='MAGE AMR')
    
    # Highlight white → AMR discordant cases
    if len(white_to_amr) > 0:
        ax2.scatter(white_to_amr['PC1'], white_to_amr['PC2'],
                   c='red', s=50, alpha=0.8, edgecolor='black', linewidth=1,
                   label=f'ADNI "White" → AMR (n={len(white_to_amr)})')
    
    # Add white → EUR cases for comparison
    white_eur = adni_assigned[
        (adni_assigned['self_reported_ethnicity'] == 'white') &
        (adni_assigned['final_consensus_ancestry'] == 'EUR')
    ]
    if len(white_eur) > 0:
        ax2.scatter(white_eur['PC1'], white_eur['PC2'],
                   c='green', s=50, alpha=0.8, edgecolor='black', linewidth=1,
                   label=f'ADNI "White" → EUR (n={len(white_eur)})')
    
    ax2.set_xlabel('PC1 (45.26%)')
    ax2.set_ylabel('PC2 (13.79%)')
    ax2.set_title('Discordant Cases: "White" ADNI Participants')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # Plot 3: Confidence distribution
    ax3 = axes[1, 0]
    
    datasets_to_compare = ['ADNI', 'GTEx']
    colors = ['red', 'blue']
    
    for dataset, color in zip(datasets_to_compare, colors):
        dataset_samples = all_data[
            (all_data['dataset'] == dataset) &
            (all_data['final_consensus_confidence'].notna())
        ]
        
        if len(dataset_samples) > 0:
            confidence_scores = dataset_samples['final_consensus_confidence']
            ax3.hist(confidence_scores, bins=20, alpha=0.6, color=color,
                    label=f'{dataset} (n={len(dataset_samples)}, mean={confidence_scores.mean():.3f})')
    
    ax3.axvline(x=0.7, color='black', linestyle='--', label='High Conf Threshold')
    ax3.set_xlabel('Assignment Confidence')
    ax3.set_ylabel('Number of Samples')
    ax3.set_title('Confidence Score Distribution by Dataset')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Plot 4: Method breakdown
    ax4 = axes[1, 1]
    
    method_data = []
    for dataset in ['ADNI', 'GTEx']:
        dataset_samples = all_data[all_data['dataset'] == dataset]
        method_counts = dataset_samples['final_consensus_method'].value_counts()
        for method, count in method_counts.items():
            method_data.append({
                'Dataset': dataset,
                'Method': method,
                'Count': count,
                'Percentage': count / len(dataset_samples) * 100
            })
    
    method_df = pd.DataFrame(method_data)
    
    # Create grouped bar chart
    methods = method_df['Method'].unique()
    x = np.arange(len(methods))
    width = 0.35
    
    adni_counts = [method_df[(method_df['Dataset'] == 'ADNI') & (method_df['Method'] == m)]['Count'].sum() 
                   for m in methods]
    gtex_counts = [method_df[(method_df['Dataset'] == 'GTEx') & (method_df['Method'] == m)]['Count'].sum() 
                   for m in methods]
    
    ax4.bar(x - width/2, adni_counts, width, label='ADNI', color='red', alpha=0.7)
    ax4.bar(x + width/2, gtex_counts, width, label='GTEx', color='blue', alpha=0.7)
    
    ax4.set_xlabel('Assignment Method')
    ax4.set_ylabel('Number of Samples')
    ax4.set_title('Assignment Method Distribution')
    ax4.set_xticks(x)
    ax4.set_xticklabels(methods, rotation=45, ha='right')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save plot
    plot_file = results_dir / "adni_ancestry_investigation.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"Saved visualization to {plot_file}")
    
    plt.show()

def main():
    """Main investigation function"""
    print("Starting ADNI Ancestry Investigation...")
    
    # Run investigation
    adni_data, white_to_amr = investigate_adni_ancestry()
    
    # Create visualization
    create_adni_visualization(adni_data, white_to_amr)
    
    print("\nInvestigation complete!")
    print(f"Results saved to: {results_dir}")

if __name__ == "__main__":
    main()