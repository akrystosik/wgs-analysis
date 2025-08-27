#!/usr/bin/env python
import os
import subprocess
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from tqdm import tqdm
import logging
from datetime import datetime

# Set up logging
log_dir = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/wgs/data/qc/variants/logs"
os.makedirs(log_dir, exist_ok=True)
log_file = os.path.join(log_dir, f"targeted_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")
logging.basicConfig(
    filename=log_file,
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# Define directories
cell_line_dir = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/wgs/data/variants/cell_line_analysis_v5"
output_dir = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/wgs/data/qc/variants/targeted_analysis"
plots_dir = os.path.join(output_dir, "plots")
summary_file = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/wgs/data/qc/variants/bcftools_summary_report.tsv"

# Create output directories
os.makedirs(output_dir, exist_ok=True)
os.makedirs(plots_dir, exist_ok=True)

def load_cell_line_data():
    """Load cell line data from the summary file"""
    if not os.path.exists(summary_file):
        print(f"Summary file not found: {summary_file}")
        logging.error(f"Summary file not found: {summary_file}")
        return None
        
    df = pd.read_csv(summary_file, sep='\t')
    print(f"Loaded data for {len(df)} samples")
    logging.info(f"Loaded data for {len(df)} samples")
    return df

def analyze_gtex_vs_cell_lines(df):
    """Analyze GTEx vs cell lines"""
    print("Analyzing GTEx vs cell lines...")
    logging.info("Analyzing GTEx vs cell lines...")
    
    # Extract GTEx row and cell line rows
    gtex_row = df[df['cell_line'] == 'GTEx']
    cell_line_rows = df[df['cell_line'] != 'GTEx']
    
    if len(gtex_row) == 0:
        print("GTEx data not found in summary file")
        logging.error("GTEx data not found in summary file")
        return
        
    # Calculate average per-sample counts for GTEx
    gtex_samples = gtex_row['num_samples'].values[0]
    gtex_snps_per_sample = gtex_row['num_SNPs'].values[0] / gtex_samples
    gtex_indels_per_sample = gtex_row['num_indels'].values[0] / gtex_samples
    
    print(f"GTEx average per individual: {gtex_snps_per_sample:.2f} SNPs, {gtex_indels_per_sample:.2f} indels")
    logging.info(f"GTEx average per individual: {gtex_snps_per_sample:.2f} SNPs, {gtex_indels_per_sample:.2f} indels")
    
    # Create comparison DataFrame
    comparison_data = []
    
    # Add individual GTEx entry
    comparison_data.append({
        'sample': 'GTEx (per individual)',
        'type': 'GTEx',
        'num_SNPs': gtex_snps_per_sample,
        'num_indels': gtex_indels_per_sample
    })
    
    # Add cell lines
    for _, row in cell_line_rows.iterrows():
        comparison_data.append({
            'sample': row['cell_line'],
            'type': 'Cell Line',
            'num_SNPs': row['num_SNPs'],
            'num_indels': row['num_indels']
        })
    
    comparison_df = pd.DataFrame(comparison_data)
    
    # Save comparison data
    comparison_file = os.path.join(output_dir, "gtex_vs_cell_lines.tsv")
    comparison_df.to_csv(comparison_file, sep='\t', index=False)
    print(f"Comparison data saved to {comparison_file}")
    logging.info(f"Comparison data saved to {comparison_file}")
    
    return comparison_df

def visualize_comparisons(comparison_df):
    """Create visualizations for the comparisons"""
    print("Creating visualizations...")
    logging.info("Creating visualizations...")
    
    # Set seaborn style
    sns.set(style="whitegrid")
    
    # 1. SNP count comparison
    plt.figure(figsize=(12, 8))
    sns.barplot(x='sample', y='num_SNPs', hue='type', data=comparison_df, palette=['blue', 'red'])
    plt.title('SNP Counts: GTEx Average vs Cell Lines')
    plt.xlabel('Sample')
    plt.ylabel('Number of SNPs')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, 'snp_comparison.png'), dpi=300)
    plt.close()
    
    # 2. Indel count comparison
    plt.figure(figsize=(12, 8))
    sns.barplot(x='sample', y='num_indels', hue='type', data=comparison_df, palette=['blue', 'red'])
    plt.title('Indel Counts: GTEx Average vs Cell Lines')
    plt.xlabel('Sample')
    plt.ylabel('Number of Indels')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, 'indel_comparison.png'), dpi=300)
    plt.close()
    
    # 3. SNP/Indel ratio
    comparison_df['snp_indel_ratio'] = comparison_df['num_SNPs'] / comparison_df['num_indels']
    
    plt.figure(figsize=(12, 8))
    sns.barplot(x='sample', y='snp_indel_ratio', hue='type', data=comparison_df, palette=['blue', 'red'])
    plt.title('SNP/Indel Ratio: GTEx Average vs Cell Lines')
    plt.xlabel('Sample')
    plt.ylabel('SNP/Indel Ratio')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, 'ratio_comparison.png'), dpi=300)
    plt.close()
    
    # 4. Cell line deviation from GTEx average
    plt.figure(figsize=(12, 8))
    
    # Get GTEx average
    gtex_snps = comparison_df[comparison_df['type'] == 'GTEx']['num_SNPs'].values[0]
    
    # Calculate fold change for cell lines
    cell_lines = comparison_df[comparison_df['type'] == 'Cell Line'].copy()
    cell_lines['fold_change'] = cell_lines['num_SNPs'] / gtex_snps
    
    # Sort by fold change
    cell_lines = cell_lines.sort_values('fold_change', ascending=False)
    
    # Plot fold change
    sns.barplot(x='sample', y='fold_change', data=cell_lines, palette='YlOrRd')
    plt.axhline(y=1, color='blue', linestyle='--', label='GTEx Average')
    plt.title('Cell Line SNP Count Relative to GTEx Average (Fold Change)')
    plt.xlabel('Cell Line')
    plt.ylabel('Fold Change')
    plt.xticks(rotation=45, ha='right')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, 'fold_change.png'), dpi=300)
    plt.close()
    
    print(f"Visualizations saved to {plots_dir}")
    logging.info(f"Visualizations saved to {plots_dir}")

def analyze_cell_line_patterns(df):
    """Analyze patterns across cell lines"""
    print("Analyzing cell line patterns...")
    logging.info("Analyzing cell line patterns...")
    
    # Filter to only cell lines
    cell_lines = df[df['cell_line'] != 'GTEx'].copy()
    
    if len(cell_lines) == 0:
        print("No cell line data found")
        logging.error("No cell line data found")
        return
    
    # Calculate additional metrics
    cell_lines['snp_indel_ratio'] = cell_lines['num_SNPs'] / cell_lines['num_indels']
    
    # Save cell line analysis
    cell_lines_file = os.path.join(output_dir, "cell_line_analysis.tsv")
    cell_lines.to_csv(cell_lines_file, sep='\t', index=False)
    print(f"Cell line analysis saved to {cell_lines_file}")
    logging.info(f"Cell line analysis saved to {cell_lines_file}")
    
    # Visualize cell line comparison
    plt.figure(figsize=(14, 8))
    sns.scatterplot(x='num_SNPs', y='num_indels', data=cell_lines, s=100)
    
    # Add cell line labels
    for _, row in cell_lines.iterrows():
        plt.text(row['num_SNPs']*1.01, row['num_indels']*1.01, row['cell_line'])
    
    plt.title('Cell Line SNPs vs Indels')
    plt.xlabel('Number of SNPs')
    plt.ylabel('Number of Indels')
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, 'cell_line_comparison.png'), dpi=300)
    plt.close()
    
    # Create a correlation heatmap
    plt.figure(figsize=(10, 8))
    correlation_columns = ['num_SNPs', 'num_indels', 'snp_indel_ratio']
    correlation = cell_lines[correlation_columns].corr()
    sns.heatmap(correlation, annot=True, cmap='coolwarm', vmin=-1, vmax=1)
    plt.title('Correlation Between Variant Metrics')
    plt.tight_layout()
    plt.savefig(os.path.join(plots_dir, 'correlation.png'), dpi=300)
    plt.close()
    
    # Print summary statistics
    print("\nCell Line Summary Statistics:")
    print(f"Average SNPs: {cell_lines['num_SNPs'].mean():,.2f}")
    print(f"Average Indels: {cell_lines['num_indels'].mean():,.2f}")
    print(f"Average SNP/Indel Ratio: {cell_lines['snp_indel_ratio'].mean():.4f}")
    print(f"Cell line with most SNPs: {cell_lines.loc[cell_lines['num_SNPs'].idxmax(), 'cell_line']}")
    print(f"Cell line with most indels: {cell_lines.loc[cell_lines['num_indels'].idxmax(), 'cell_line']}")
    
    logging.info("\nCell Line Summary Statistics:")
    logging.info(f"Average SNPs: {cell_lines['num_SNPs'].mean():,.2f}")
    logging.info(f"Average Indels: {cell_lines['num_indels'].mean():,.2f}")
    logging.info(f"Average SNP/Indel Ratio: {cell_lines['snp_indel_ratio'].mean():.4f}")
    logging.info(f"Cell line with most SNPs: {cell_lines.loc[cell_lines['num_SNPs'].idxmax(), 'cell_line']}")
    logging.info(f"Cell line with most indels: {cell_lines.loc[cell_lines['num_indels'].idxmax(), 'cell_line']}")

def main():
    # Load data from summary file
    df = load_cell_line_data()
    
    if df is not None:
        # Analyze GTEx vs cell lines
        comparison_df = analyze_gtex_vs_cell_lines(df)
        
        if comparison_df is not None:
            # Visualize comparisons
            visualize_comparisons(comparison_df)
        
        # Analyze cell line patterns
        analyze_cell_line_patterns(df)
    
    print(f"Analysis complete. Results saved to {output_dir}")
    logging.info(f"Analysis complete. Results saved to {output_dir}")

if __name__ == "__main__":
    main()