#!/usr/bin/env python3
"""
Real VCF-based PCA analysis using actual MAGE data
Efficiently processes the large VCF file by extracting subsets
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import warnings
import subprocess
import tempfile
import os
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import argparse
import time

warnings.filterwarnings('ignore')

def load_mage_ancestry_data():
    """Load MAGE ancestry assignment data"""
    print("Loading MAGE ancestry data...")
    
    # Load the 1000 Genomes population data
    pop_file = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/mage/20130606_g1k_3202_samples_ped_population.txt"
    pop_df = pd.read_csv(pop_file, sep=' ')
    
    # Load MAGE metadata
    metadata_file = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/rnaseq/mage/sample.metadata.MAGE.v1.0.txt"
    metadata_df = pd.read_csv(metadata_file, sep='\t')
    
    # Create sample to ancestry mapping
    sample_to_ancestry = {}
    sample_to_population = {}
    
    # From population file
    for _, row in pop_df.iterrows():
        sample_id = row['SampleID']
        sample_to_ancestry[sample_id] = row['Superpopulation']
        sample_to_population[sample_id] = row['Population']
    
    # From metadata file
    for _, row in metadata_df.iterrows():
        sample_id = row['sample_kgpID']
        sample_to_ancestry[sample_id] = row['continentalGroup']
        sample_to_population[sample_id] = row['population']
    
    print(f"Loaded ancestry data for {len(sample_to_ancestry)} samples")
    return sample_to_ancestry, sample_to_population

def extract_vcf_samples(vcf_path):
    """Extract all sample names from VCF header"""
    print("Extracting sample names from VCF...")
    
    cmd = f"zcat {vcf_path} | grep '^#CHROM' | head -1"
    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=120)
        
        if result.returncode == 0:
            header_line = result.stdout.strip()
            samples = header_line.split('\t')[9:]  # Skip first 9 VCF columns
            print(f"Found {len(samples)} samples in VCF")
            return samples
        else:
            print(f"Error extracting samples: {result.stderr}")
            return []
            
    except subprocess.TimeoutExpired:
        print("Timeout extracting samples")
        return []

def create_vcf_subset(vcf_path, output_path, chromosome="chr22", max_lines=5000):
    """Create a subset VCF for a specific chromosome"""
    print(f"Creating VCF subset for {chromosome} with max {max_lines} variants...")
    
    cmd = f"""
    zcat {vcf_path} | head -200 | grep '^##' > {output_path} && \
    zcat {vcf_path} | grep '^#CHROM' | head -1 >> {output_path} && \
    zcat {vcf_path} | grep '^{chromosome}\\s' | head -{max_lines} >> {output_path}
    """
    
    try:
        result = subprocess.run(cmd, shell=True, timeout=300)
        
        if result.returncode == 0 and os.path.exists(output_path):
            # Check if we got data
            line_count = sum(1 for line in open(output_path) if not line.startswith('##'))
            print(f"✓ Created subset VCF with {line_count-1} variants")
            return line_count > 1
        else:
            print("Failed to create VCF subset")
            return False
            
    except subprocess.TimeoutExpired:
        print("Timeout creating VCF subset")
        return False

def parse_vcf_to_matrix(vcf_path):
    """Parse VCF file and convert to genotype matrix"""
    print("Parsing VCF to genotype matrix...")
    
    samples = []
    variants = []
    genotype_data = []
    
    with open(vcf_path, 'r') as f:
        for line_num, line in enumerate(f):
            if line_num % 1000 == 0:
                print(f"  Processing line {line_num}...", end='\r')
                
            line = line.strip()
            if not line:
                continue
                
            if line.startswith('##'):
                continue
            elif line.startswith('#CHROM'):
                # Header line with sample names
                fields = line.split('\t')
                samples = fields[9:]  # Skip first 9 VCF columns
                print(f"\nFound {len(samples)} samples")
                continue
            else:
                # Variant line
                fields = line.split('\t')
                if len(fields) < 10:
                    continue
                
                chrom = fields[0]
                pos = fields[1]
                ref = fields[3]
                alt = fields[4]
                
                # Extract genotypes
                genotypes = []
                for i in range(9, len(fields)):
                    gt_field = fields[i].split(':')[0]  # Get GT field
                    
                    # Convert genotype to allele count
                    if gt_field in ['0/0', '0|0']:
                        genotypes.append(0)
                    elif gt_field in ['0/1', '1/0', '0|1', '1|0']:
                        genotypes.append(1)
                    elif gt_field in ['1/1', '1|1']:
                        genotypes.append(2)
                    else:
                        genotypes.append(-1)  # Missing
                
                variants.append(f"{chrom}:{pos}_{ref}>{alt}")
                genotype_data.append(genotypes)
    
    print(f"\nParsed {len(variants)} variants across {len(samples)} samples")
    
    if not genotype_data:
        raise Exception("No genotype data found")
    
    # Convert to numpy array
    gt_matrix = np.array(genotype_data, dtype=float)
    
    # Handle missing data (replace -1 with NaN, then impute with mean)
    gt_matrix[gt_matrix == -1] = np.nan
    
    # Impute missing values with variant mean
    for i in range(gt_matrix.shape[0]):
        variant_data = gt_matrix[i, :]
        if np.any(np.isnan(variant_data)):
            mean_val = np.nanmean(variant_data)
            if not np.isnan(mean_val):
                variant_data[np.isnan(variant_data)] = mean_val
                gt_matrix[i, :] = variant_data
            else:
                # If all missing, set to 0
                variant_data[np.isnan(variant_data)] = 0
                gt_matrix[i, :] = variant_data
    
    # Transpose so samples are rows, variants are columns
    gt_matrix = gt_matrix.T
    
    print(f"Final genotype matrix: {gt_matrix.shape} (samples x variants)")
    print(f"Missing values remaining: {np.sum(np.isnan(gt_matrix))}")
    
    return gt_matrix, samples, variants

def perform_pca_analysis(gt_matrix, n_components=10):
    """Perform PCA on genotype matrix"""
    print(f"Performing PCA with up to {n_components} components...")
    
    # Remove variants with no variation (MAF = 0)
    variant_vars = np.var(gt_matrix, axis=0)
    var_mask = variant_vars > 0
    gt_filtered = gt_matrix[:, var_mask]
    
    print(f"Filtered from {gt_matrix.shape[1]} to {gt_filtered.shape[1]} variable variants")
    
    if gt_filtered.shape[1] == 0:
        raise Exception("No variable variants found")
    
    # Standardize the genotype data
    scaler = StandardScaler()
    gt_scaled = scaler.fit_transform(gt_filtered)
    
    # Determine actual number of components
    max_components = min(n_components, gt_scaled.shape[0]-1, gt_scaled.shape[1])
    
    # Perform PCA
    pca = PCA(n_components=max_components, random_state=42)
    pca_result = pca.fit_transform(gt_scaled)
    
    print(f"Generated {pca_result.shape[1]} principal components")
    print("\nExplained variance ratio:")
    for i, var_ratio in enumerate(pca.explained_variance_ratio_):
        print(f"  PC{i+1}: {var_ratio:.4f} ({var_ratio*100:.2f}%)")
    
    total_var = np.sum(pca.explained_variance_ratio_)
    print(f"Total explained variance: {total_var:.4f} ({total_var*100:.2f}%)")
    
    return pca_result, pca.explained_variance_ratio_

def create_comprehensive_visualizations(pca_result, samples, sample_to_ancestry, sample_to_population, explained_variance, output_dir):
    """Create comprehensive PCA visualizations with real data"""
    print("Creating comprehensive PCA visualizations...")
    
    # Create DataFrame for analysis
    pca_df = pd.DataFrame({
        'sample_id': samples,
        'PC1': pca_result[:, 0],
        'PC2': pca_result[:, 1]
    })
    
    # Add additional PCs if available
    for i in range(2, min(pca_result.shape[1], 6)):
        pca_df[f'PC{i+1}'] = pca_result[:, i]
    
    # Map ancestry and population
    pca_df['ancestry'] = pca_df['sample_id'].map(sample_to_ancestry).fillna('Unknown')
    pca_df['population'] = pca_df['sample_id'].map(sample_to_population).fillna('Unknown')
    
    # Count mappings
    mapped_count = (pca_df['ancestry'] != 'Unknown').sum()
    print(f"Successfully mapped ancestry for {mapped_count}/{len(pca_df)} samples ({mapped_count/len(pca_df)*100:.1f}%)")
    
    print("Ancestry distribution in VCF samples:")
    ancestry_counts = pca_df['ancestry'].value_counts()
    for ancestry, count in ancestry_counts.items():
        print(f"  {ancestry}: {count}")
    
    # Create output directory
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Define colors
    ancestry_colors = {
        'EUR': '#1f77b4', 'AFR': '#ff7f0e', 'EAS': '#2ca02c', 
        'AMR': '#d62728', 'SAS': '#9467bd', 'Unknown': '#7f7f7f'
    }
    
    # 1. Main PCA plot
    plt.figure(figsize=(14, 10))
    
    for ancestry in sorted(pca_df['ancestry'].unique()):
        mask = pca_df['ancestry'] == ancestry
        plt.scatter(pca_df.loc[mask, 'PC1'], pca_df.loc[mask, 'PC2'],
                   c=ancestry_colors.get(ancestry, '#cccccc'),
                   label=f'{ancestry} (n={mask.sum()})',
                   alpha=0.8, s=50, edgecolors='black', linewidth=0.5)
    
    plt.xlabel(f'PC1 ({explained_variance[0]*100:.1f}% variance)', fontsize=12)
    plt.ylabel(f'PC2 ({explained_variance[1]*100:.1f}% variance)', fontsize=12)
    plt.title('PCA of Real VCF Genotype Data\nColored by MAGE Continental Ancestry', 
              fontsize=14, fontweight='bold')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_dir / 'real_vcf_pca_ancestry.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Population-level detailed plot
    plt.figure(figsize=(16, 12))
    
    # Get unique populations with sufficient samples
    pop_counts = pca_df['population'].value_counts()
    major_pops = pop_counts[pop_counts >= 5].index  # Only populations with 5+ samples
    
    colors = plt.cm.tab20(np.linspace(0, 1, len(major_pops)))
    
    for i, population in enumerate(major_pops):
        if population == 'Unknown':
            continue
        mask = pca_df['population'] == population
        plt.scatter(pca_df.loc[mask, 'PC1'], pca_df.loc[mask, 'PC2'],
                   c=[colors[i]], label=f'{population} (n={mask.sum()})',
                   alpha=0.7, s=40, edgecolors='black', linewidth=0.3)
    
    plt.xlabel(f'PC1 ({explained_variance[0]*100:.1f}% variance)', fontsize=12)
    plt.ylabel(f'PC2 ({explained_variance[1]*100:.1f}% variance)', fontsize=12)
    plt.title('PCA of Real VCF Data by Specific Populations\n(Populations with ≥5 samples)', 
              fontsize=14, fontweight='bold')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(output_dir / 'real_vcf_pca_populations.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 3. Multi-component comparison
    if pca_result.shape[1] >= 4:
        fig, axes = plt.subplots(2, 2, figsize=(16, 14))
        axes = axes.flatten()
        
        pc_pairs = [(0, 1), (0, 2), (1, 2), (2, 3)] if pca_result.shape[1] > 3 else [(0, 1), (0, 2), (1, 2), (0, 1)]
        
        for idx, (pc1, pc2) in enumerate(pc_pairs):
            if pc2 >= pca_result.shape[1]:
                continue
                
            ax = axes[idx]
            for ancestry in sorted(pca_df['ancestry'].unique()):
                mask = pca_df['ancestry'] == ancestry
                ax.scatter(pca_df.loc[mask, f'PC{pc1+1}'], pca_df.loc[mask, f'PC{pc2+1}'],
                          c=ancestry_colors.get(ancestry, '#cccccc'),
                          label=ancestry if idx == 0 else "", alpha=0.7, s=25)
            
            ax.set_xlabel(f'PC{pc1+1} ({explained_variance[pc1]*100:.1f}%)')
            ax.set_ylabel(f'PC{pc2+1} ({explained_variance[pc2]*100:.1f}%)')
            ax.set_title(f'PC{pc1+1} vs PC{pc2+1}')
            ax.grid(True, alpha=0.3)
            
            if idx == 0:
                ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        
        plt.suptitle('Multi-Component PCA Analysis of Real VCF Data', fontsize=16, fontweight='bold')
        plt.tight_layout()
        plt.savefig(output_dir / 'real_vcf_pca_multicomponent.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    # 4. Analysis summary dashboard
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    
    # Ancestry pie chart
    ancestry_counts.plot(kind='pie', ax=ax1, autopct='%1.1f%%',
                        colors=[ancestry_colors.get(x, '#cccccc') for x in ancestry_counts.index])
    ax1.set_title('Sample Distribution by Ancestry\n(Real VCF Data)', fontweight='bold')
    ax1.set_ylabel('')
    
    # PC1 histogram by ancestry
    for ancestry in sorted(pca_df['ancestry'].unique()):
        mask = pca_df['ancestry'] == ancestry
        ax2.hist(pca_df.loc[mask, 'PC1'], alpha=0.6, label=ancestry, bins=15,
                color=ancestry_colors.get(ancestry, '#cccccc'))
    ax2.set_xlabel('PC1')
    ax2.set_ylabel('Frequency')
    ax2.set_title('PC1 Distribution by Ancestry', fontweight='bold')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # PC2 histogram by ancestry
    for ancestry in sorted(pca_df['ancestry'].unique()):
        mask = pca_df['ancestry'] == ancestry
        ax3.hist(pca_df.loc[mask, 'PC2'], alpha=0.6, label=ancestry, bins=15,
                color=ancestry_colors.get(ancestry, '#cccccc'))
    ax3.set_xlabel('PC2')
    ax3.set_ylabel('Frequency')
    ax3.set_title('PC2 Distribution by Ancestry', fontweight='bold')
    ax3.legend()
    ax3.grid(True, alpha=0.3)
    
    # Explained variance
    n_show = min(10, len(explained_variance))
    pc_names = [f'PC{i+1}' for i in range(n_show)]
    ax4.bar(pc_names, explained_variance[:n_show] * 100, 
            color='skyblue', edgecolor='black', alpha=0.8)
    ax4.set_ylabel('Explained Variance (%)')
    ax4.set_title('PCA Explained Variance\n(Real VCF Data)', fontweight='bold')
    ax4.tick_params(axis='x', rotation=45)
    ax4.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'real_vcf_pca_dashboard.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✅ Saved all visualizations to: {output_dir}")
    
    return pca_df

def main():
    parser = argparse.ArgumentParser(description='Real VCF-based PCA analysis with MAGE ancestry')
    parser.add_argument('--vcf-path', 
                       default='/mnt/czi-sci-ai/intrinsic-variation-gene-ex/wgs/data/variants/all_vcfs/merged.all.biallelic.maf0.01.vcf.gz',
                       help='Path to VCF file')
    parser.add_argument('--output-dir', default='real_vcf_pca_results',
                       help='Output directory')
    parser.add_argument('--chromosome', default='chr22',
                       help='Chromosome to analyze (default: chr22)')
    parser.add_argument('--max-variants', type=int, default=5000,
                       help='Maximum variants to extract')
    parser.add_argument('--n-components', type=int, default=10,
                       help='Number of PCA components')
    
    args = parser.parse_args()
    
    print("="*70)
    print("REAL VCF-BASED PCA ANALYSIS WITH MAGE ANCESTRY DATA")
    print("="*70)
    print(f"📁 VCF file: {args.vcf_path}")
    print(f"🧬 Analyzing chromosome: {args.chromosome}")
    print(f"📊 Max variants: {args.max_variants}")
    print("="*70)
    
    start_time = time.time()
    
    try:
        # Load ancestry mapping
        sample_to_ancestry, sample_to_population = load_mage_ancestry_data()
        
        # Create temporary subset VCF
        with tempfile.NamedTemporaryFile(suffix='.vcf', delete=False) as tmp_vcf:
            tmp_vcf_path = tmp_vcf.name
        
        try:
            # Extract VCF subset
            if not create_vcf_subset(args.vcf_path, tmp_vcf_path, args.chromosome, args.max_variants):
                raise Exception("Failed to create VCF subset")
            
            # Parse VCF to matrix
            gt_matrix, samples, variants = parse_vcf_to_matrix(tmp_vcf_path)
            
            # Perform PCA
            pca_result, explained_variance = perform_pca_analysis(gt_matrix, args.n_components)
            
            # Create visualizations
            pca_df = create_comprehensive_visualizations(
                pca_result, samples, sample_to_ancestry, sample_to_population, 
                explained_variance, args.output_dir
            )
            
            # Save detailed results
            output_dir = Path(args.output_dir)
            pca_df.to_csv(output_dir / 'real_vcf_pca_results.csv', index=False)
            
            # Create analysis report
            with open(output_dir / 'analysis_report.txt', 'w') as f:
                f.write("REAL VCF PCA ANALYSIS REPORT\n")
                f.write("="*50 + "\n\n")
                f.write(f"VCF file: {args.vcf_path}\n")
                f.write(f"Chromosome analyzed: {args.chromosome}\n")
                f.write(f"Variants used: {len(variants)}\n")
                f.write(f"Samples analyzed: {len(samples)}\n")
                f.write(f"Principal components: {len(explained_variance)}\n\n")
                
                f.write("Ancestry mapping success:\n")
                mapped_count = (pca_df['ancestry'] != 'Unknown').sum()
                f.write(f"  Mapped: {mapped_count}/{len(pca_df)} ({mapped_count/len(pca_df)*100:.1f}%)\n\n")
                
                f.write("Ancestry distribution:\n")
                ancestry_counts = pca_df['ancestry'].value_counts()
                for ancestry, count in ancestry_counts.items():
                    f.write(f"  {ancestry}: {count}\n")
                
                f.write(f"\nExplained variance:\n")
                for i, var_ratio in enumerate(explained_variance):
                    f.write(f"  PC{i+1}: {var_ratio:.4f} ({var_ratio*100:.2f}%)\n")
                
                f.write(f"\nTotal explained variance: {np.sum(explained_variance):.4f} ({np.sum(explained_variance)*100:.2f}%)\n")
            
            elapsed_time = time.time() - start_time
            
            print("\n" + "="*70)
            print("✅ REAL VCF PCA ANALYSIS COMPLETE!")
            print("="*70)
            print(f"⏱️  Total runtime: {elapsed_time:.1f} seconds")
            print(f"🧬 Analyzed {len(variants)} variants from {args.chromosome}")
            print(f"👥 Processed {len(samples)} samples")
            print(f"📈 Generated {len(explained_variance)} principal components")
            print(f"🎯 Ancestry mapping: {mapped_count}/{len(samples)} samples ({mapped_count/len(samples)*100:.1f}%)")
            print(f"📁 Results saved to: {output_dir}")
            print("="*70)
            
        finally:
            # Clean up temporary file
            if os.path.exists(tmp_vcf_path):
                os.unlink(tmp_vcf_path)
        
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()