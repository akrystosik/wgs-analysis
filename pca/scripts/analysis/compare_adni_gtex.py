#!/usr/bin/env python3
"""
Compare ADNI vs GTEx Population Structure
Systematic comparison of demographic and genomic ancestry patterns
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from scipy import stats

# Set up file paths
base_dir = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/wgs")
pca_dir = base_dir / "pca"
results_dir = pca_dir / "results" / "adni_gtex_comparison"
results_dir.mkdir(parents=True, exist_ok=True)

def compare_adni_gtex():
    """Compare ADNI vs GTEx population structure and demographics"""
    print("Comparing ADNI vs GTEx Population Structure")
    print("=" * 60)
    
    # Load complete results
    complete_file = pca_dir / "results" / "combined_analysis" / "complete_ancestry_results.csv"
    data = pd.read_csv(complete_file)
    
    # Extract ADNI and GTEx data
    adni_data = data[data['dataset'] == 'ADNI'].copy()
    gtex_data = data[data['dataset'] == 'GTEx'].copy()
    
    print(f"ADNI samples: {len(adni_data)}")
    print(f"GTEx samples: {len(gtex_data)}")
    
    # 1. SELF-REPORTED ETHNICITY COMPARISON
    print("\n" + "="*60)
    print("SELF-REPORTED ETHNICITY COMPARISON")
    print("="*60)
    
    print("\nADNI Self-Reported Ethnicity:")
    adni_ethnicity = adni_data['self_reported_ethnicity'].value_counts()
    adni_total = adni_ethnicity.sum()
    for ethnicity, count in adni_ethnicity.items():
        pct = count / adni_total * 100
        print(f"  {ethnicity}: {count} ({pct:.1f}%)")
    
    print("\nGTEx Self-Reported Ethnicity:")
    gtex_ethnicity = gtex_data['self_reported_ethnicity'].value_counts()
    gtex_total = gtex_ethnicity.sum()
    for ethnicity, count in gtex_ethnicity.items():
        pct = count / gtex_total * 100
        print(f"  {ethnicity}: {count} ({pct:.1f}%)")
    
    # Compare white populations specifically
    print("\n--- 'White' Population Analysis ---")
    
    adni_white = adni_data[adni_data['self_reported_ethnicity'] == 'white']
    gtex_white = gtex_data[gtex_data['self_reported_ethnicity'] == 'white']
    
    print(f"\nADNI 'white' participants: {len(adni_white)} ({len(adni_white)/len(adni_data)*100:.1f}%)")
    print(f"GTEx 'white' participants: {len(gtex_white)} ({len(gtex_white)/len(gtex_data)*100:.1f}%)")
    
    # 2. GENOMIC ANCESTRY COMPARISON
    print("\n" + "="*60)
    print("GENOMIC ANCESTRY COMPARISON")
    print("="*60)
    
    # Compare final genomic ancestry assignments
    adni_assigned = adni_data[
        (adni_data['final_consensus_ancestry'] != 'Unknown') &
        (adni_data['final_consensus_ancestry'] != 'Uncertain')
    ]
    gtex_assigned = gtex_data[
        (gtex_data['final_consensus_ancestry'] != 'Unknown') &
        (gtex_data['final_consensus_ancestry'] != 'Uncertain')
    ]
    
    print(f"\nADNI Genomic Ancestry (assigned samples: {len(adni_assigned)}):")
    adni_genomic = adni_assigned['final_consensus_ancestry'].value_counts()
    for ancestry, count in adni_genomic.items():
        pct = count / len(adni_assigned) * 100
        print(f"  {ancestry}: {count} ({pct:.1f}%)")
    
    print(f"\nGTEx Genomic Ancestry (assigned samples: {len(gtex_assigned)}):")
    gtex_genomic = gtex_assigned['final_consensus_ancestry'].value_counts()
    for ancestry, count in gtex_genomic.items():
        pct = count / len(gtex_assigned) * 100
        print(f"  {ancestry}: {count} ({pct:.1f}%)")
    
    # 3. "WHITE" POPULATION GENOMIC BREAKDOWN
    print("\n" + "="*60)
    print("'WHITE' POPULATION GENOMIC ANCESTRY BREAKDOWN")
    print("="*60)
    
    # Filter for assigned samples only
    adni_white_assigned = adni_white[
        (adni_white['final_consensus_ancestry'] != 'Unknown') &
        (adni_white['final_consensus_ancestry'] != 'Uncertain')
    ]
    gtex_white_assigned = gtex_white[
        (gtex_white['final_consensus_ancestry'] != 'Unknown') &
        (gtex_white['final_consensus_ancestry'] != 'Uncertain')
    ]
    
    print(f"\nADNI 'white' → Genomic ancestry breakdown (n={len(adni_white_assigned)}):")
    adni_white_genomic = adni_white_assigned['final_consensus_ancestry'].value_counts()
    for ancestry, count in adni_white_genomic.items():
        pct = count / len(adni_white_assigned) * 100
        print(f"  {ancestry}: {count} ({pct:.1f}%)")
    
    print(f"\nGTEx 'white' → Genomic ancestry breakdown (n={len(gtex_white_assigned)}):")
    gtex_white_genomic = gtex_white_assigned['final_consensus_ancestry'].value_counts()
    for ancestry, count in gtex_white_genomic.items():
        pct = count / len(gtex_white_assigned) * 100
        print(f"  {ancestry}: {count} ({pct:.1f}%)")
    
    # 4. PC COORDINATE ANALYSIS
    print("\n" + "="*60)
    print("PC COORDINATE DISTRIBUTION ANALYSIS")
    print("="*60)
    
    # Compare PC coordinates of white populations
    print("\nPC Coordinate Statistics for 'White' Populations:")
    
    print(f"\nADNI 'white' samples (n={len(adni_white_assigned)}):")
    print(f"  PC1: {adni_white_assigned['PC1'].mean():.5f} ± {adni_white_assigned['PC1'].std():.5f}")
    print(f"  PC2: {adni_white_assigned['PC2'].mean():.5f} ± {adni_white_assigned['PC2'].std():.5f}")
    print(f"  PC1 range: {adni_white_assigned['PC1'].min():.5f} to {adni_white_assigned['PC1'].max():.5f}")
    print(f"  PC2 range: {adni_white_assigned['PC2'].min():.5f} to {adni_white_assigned['PC2'].max():.5f}")
    
    print(f"\nGTEx 'white' samples (n={len(gtex_white_assigned)}):")
    print(f"  PC1: {gtex_white_assigned['PC1'].mean():.5f} ± {gtex_white_assigned['PC1'].std():.5f}")
    print(f"  PC2: {gtex_white_assigned['PC2'].mean():.5f} ± {gtex_white_assigned['PC2'].std():.5f}")
    print(f"  PC1 range: {gtex_white_assigned['PC1'].min():.5f} to {gtex_white_assigned['PC1'].max():.5f}")
    print(f"  PC2 range: {gtex_white_assigned['PC2'].min():.5f} to {gtex_white_assigned['PC2'].max():.5f}")
    
    # Statistical comparison
    pc1_ttest = stats.ttest_ind(adni_white_assigned['PC1'], gtex_white_assigned['PC1'])
    pc2_ttest = stats.ttest_ind(adni_white_assigned['PC2'], gtex_white_assigned['PC2'])
    
    print(f"\nStatistical Tests:")
    print(f"  PC1 difference: t={pc1_ttest.statistic:.3f}, p={pc1_ttest.pvalue:.2e}")
    print(f"  PC2 difference: t={pc2_ttest.statistic:.3f}, p={pc2_ttest.pvalue:.2e}")
    
    # Effect sizes (Cohen's d)
    pc1_cohens_d = (adni_white_assigned['PC1'].mean() - gtex_white_assigned['PC1'].mean()) / \
                   np.sqrt((adni_white_assigned['PC1'].var() + gtex_white_assigned['PC1'].var()) / 2)
    pc2_cohens_d = (adni_white_assigned['PC2'].mean() - gtex_white_assigned['PC2'].mean()) / \
                   np.sqrt((adni_white_assigned['PC2'].var() + gtex_white_assigned['PC2'].var()) / 2)
    
    print(f"  PC1 effect size (Cohen's d): {pc1_cohens_d:.3f}")
    print(f"  PC2 effect size (Cohen's d): {pc2_cohens_d:.3f}")
    
    # 5. CONFIDENCE AND ASSIGNMENT METHOD COMPARISON
    print("\n" + "="*60)
    print("ASSIGNMENT QUALITY COMPARISON")
    print("="*60)
    
    print("\nAssignment Success Rates:")
    adni_success_rate = len(adni_assigned) / len(adni_data) * 100
    gtex_success_rate = len(gtex_assigned) / len(gtex_data) * 100
    
    print(f"  ADNI: {len(adni_assigned)}/{len(adni_data)} = {adni_success_rate:.1f}%")
    print(f"  GTEx: {len(gtex_assigned)}/{len(gtex_data)} = {gtex_success_rate:.1f}%")
    
    print("\nConfidence Score Comparison:")
    adni_confidence = adni_assigned['final_consensus_confidence'].dropna()
    gtex_confidence = gtex_assigned['final_consensus_confidence'].dropna()
    
    if len(adni_confidence) > 0 and len(gtex_confidence) > 0:
        print(f"  ADNI mean confidence: {adni_confidence.mean():.3f} ± {adni_confidence.std():.3f}")
        print(f"  GTEx mean confidence: {gtex_confidence.mean():.3f} ± {gtex_confidence.std():.3f}")
        
        conf_ttest = stats.ttest_ind(adni_confidence, gtex_confidence)
        print(f"  Confidence difference: t={conf_ttest.statistic:.3f}, p={conf_ttest.pvalue:.2e}")
    
    # 6. DATASET RECRUITMENT/DEMOGRAPHIC CONTEXT
    print("\n" + "="*60)
    print("DATASET CONTEXT AND RECRUITMENT PATTERNS")
    print("="*60)
    
    print("\nDataset Characteristics:")
    print("\nADNI (Alzheimer's Disease Neuroimaging Initiative):")
    print("  - Clinical study focused on Alzheimer's disease")
    print("  - US-based recruitment from multiple sites")
    print("  - Likely includes diverse US population including Latino/Hispanic")
    print("  - Self-reported ethnicity may not capture genetic admixture")
    
    print("\nGTEx (Genotype-Tissue Expression):")
    print("  - Post-mortem tissue donor study") 
    print("  - US-based with some geographic bias")
    print("  - Different demographic recruitment patterns")
    print("  - May have different population representation")
    
    # Return data for visualization
    return {
        'adni_data': adni_data,
        'gtex_data': gtex_data,
        'adni_white_assigned': adni_white_assigned,
        'gtex_white_assigned': gtex_white_assigned,
        'adni_assigned': adni_assigned,
        'gtex_assigned': gtex_assigned
    }

def create_comparison_visualization(comparison_data):
    """Create comprehensive ADNI vs GTEx comparison visualization"""
    print("\nCreating ADNI vs GTEx comparison visualization...")
    
    # Unpack data
    adni_data = comparison_data['adni_data']
    gtex_data = comparison_data['gtex_data']
    adni_white = comparison_data['adni_white_assigned']
    gtex_white = comparison_data['gtex_white_assigned']
    adni_assigned = comparison_data['adni_assigned']
    gtex_assigned = comparison_data['gtex_assigned']
    
    # Create figure
    fig, axes = plt.subplots(2, 3, figsize=(18, 12))
    fig.suptitle('ADNI vs GTEx Population Structure Comparison', fontsize=16, fontweight='bold')
    
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
    
    # Plot 1: Self-reported ethnicity comparison
    ax1 = axes[0, 0]
    
    # ADNI ethnicity pie chart
    adni_ethnicity_counts = adni_data['self_reported_ethnicity'].value_counts()
    colors1 = plt.cm.Set3(np.linspace(0, 1, len(adni_ethnicity_counts)))
    
    wedges, texts, autotexts = ax1.pie(adni_ethnicity_counts.values, 
                                       labels=[f"{label}\n({count})" for label, count in adni_ethnicity_counts.items()],
                                       autopct='%1.1f%%', colors=colors1, startangle=90)
    ax1.set_title('ADNI Self-Reported Ethnicity')
    
    # Plot 2: GTEx ethnicity 
    ax2 = axes[0, 1]
    gtex_ethnicity_counts = gtex_data['self_reported_ethnicity'].value_counts()
    colors2 = plt.cm.Set3(np.linspace(0, 1, len(gtex_ethnicity_counts)))
    
    wedges, texts, autotexts = ax2.pie(gtex_ethnicity_counts.values,
                                       labels=[f"{label}\n({count})" for label, count in gtex_ethnicity_counts.items()],
                                       autopct='%1.1f%%', colors=colors2, startangle=90)
    ax2.set_title('GTEx Self-Reported Ethnicity')
    
    # Plot 3: Genomic ancestry comparison
    ax3 = axes[0, 2]
    
    # Prepare data for side-by-side comparison
    all_ancestries = set(adni_assigned['final_consensus_ancestry'].unique()) | set(gtex_assigned['final_consensus_ancestry'].unique())
    
    adni_counts = []
    gtex_counts = []
    labels = []
    
    for ancestry in sorted(all_ancestries):
        adni_count = (adni_assigned['final_consensus_ancestry'] == ancestry).sum()
        gtex_count = (gtex_assigned['final_consensus_ancestry'] == ancestry).sum()
        
        adni_counts.append(adni_count)
        gtex_counts.append(gtex_count)
        labels.append(ancestry)
    
    x = np.arange(len(labels))
    width = 0.35
    
    bars1 = ax3.bar(x - width/2, adni_counts, width, label='ADNI', color='red', alpha=0.7)
    bars2 = ax3.bar(x + width/2, gtex_counts, width, label='GTEx', color='blue', alpha=0.7)
    
    ax3.set_xlabel('Genomic Ancestry')
    ax3.set_ylabel('Number of Samples')
    ax3.set_title('Genomic Ancestry Distribution')
    ax3.set_xticks(x)
    ax3.set_xticklabels(labels, rotation=45)
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Add value labels on bars
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax3.annotate(f'{int(height)}',
                           xy=(bar.get_x() + bar.get_width() / 2, height),
                           xytext=(0, 3),
                           textcoords="offset points",
                           ha='center', va='bottom', fontsize=8)
    
    # Plot 4: PC coordinate distribution of "white" populations
    ax4 = axes[1, 0]
    
    ax4.scatter(adni_white['PC1'], adni_white['PC2'], 
               c='red', alpha=0.6, s=30, label=f'ADNI "white" (n={len(adni_white)})')
    ax4.scatter(gtex_white['PC1'], gtex_white['PC2'], 
               c='blue', alpha=0.6, s=30, label=f'GTEx "white" (n={len(gtex_white)})')
    
    ax4.set_xlabel('PC1 (45.26%)')
    ax4.set_ylabel('PC2 (13.79%)')
    ax4.set_title('"White" Population PC Coordinates')
    ax4.legend()
    ax4.grid(True, alpha=0.3)
    
    # Plot 5: "White" genomic ancestry breakdown
    ax5 = axes[1, 1]
    
    # ADNI white breakdown
    adni_white_genomic = adni_white['final_consensus_ancestry'].value_counts()
    gtex_white_genomic = gtex_white['final_consensus_ancestry'].value_counts()
    
    # Prepare data
    all_white_ancestries = set(adni_white_genomic.index) | set(gtex_white_genomic.index)
    
    adni_white_counts = []
    gtex_white_counts = []
    white_labels = []
    
    for ancestry in sorted(all_white_ancestries):
        adni_count = adni_white_genomic.get(ancestry, 0)
        gtex_count = gtex_white_genomic.get(ancestry, 0)
        
        adni_white_counts.append(adni_count)
        gtex_white_counts.append(gtex_count)
        white_labels.append(ancestry)
    
    x = np.arange(len(white_labels))
    width = 0.35
    
    bars1 = ax5.bar(x - width/2, adni_white_counts, width, label='ADNI "white"', color='red', alpha=0.7)
    bars2 = ax5.bar(x + width/2, gtex_white_counts, width, label='GTEx "white"', color='blue', alpha=0.7)
    
    ax5.set_xlabel('Genomic Ancestry')
    ax5.set_ylabel('Number of Samples')
    ax5.set_title('"White" Population → Genomic Ancestry')
    ax5.set_xticks(x)
    ax5.set_xticklabels(white_labels, rotation=45)
    ax5.legend()
    ax5.grid(True, alpha=0.3)
    
    # Add value labels
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax5.annotate(f'{int(height)}',
                           xy=(bar.get_x() + bar.get_width() / 2, height),
                           xytext=(0, 3),
                           textcoords="offset points",
                           ha='center', va='bottom', fontsize=8)
    
    # Plot 6: Confidence comparison
    ax6 = axes[1, 2]
    
    adni_conf = adni_assigned['final_consensus_confidence'].dropna()
    gtex_conf = gtex_assigned['final_consensus_confidence'].dropna()
    
    ax6.hist(adni_conf, bins=20, alpha=0.6, color='red', 
            label=f'ADNI (n={len(adni_conf)}, mean={adni_conf.mean():.3f})', density=True)
    ax6.hist(gtex_conf, bins=20, alpha=0.6, color='blue',
            label=f'GTEx (n={len(gtex_conf)}, mean={gtex_conf.mean():.3f})', density=True)
    
    ax6.axvline(x=0.7, color='black', linestyle='--', label='High Confidence')
    ax6.set_xlabel('Assignment Confidence')
    ax6.set_ylabel('Density')
    ax6.set_title('Assignment Confidence Distribution')
    ax6.legend()
    ax6.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save plot
    plot_file = results_dir / "adni_gtex_comparison.png"
    plt.savefig(plot_file, dpi=300, bbox_inches='tight')
    print(f"Saved comparison plot to {plot_file}")
    
    plt.show()

def main():
    """Main comparison function"""
    print("Starting ADNI vs GTEx Comparison Analysis...")
    
    # Run comparison
    comparison_data = compare_adni_gtex()
    
    # Create visualization
    create_comparison_visualization(comparison_data)
    
    print(f"\nComparison analysis complete!")
    print(f"Results saved to: {results_dir}")

if __name__ == "__main__":
    main()