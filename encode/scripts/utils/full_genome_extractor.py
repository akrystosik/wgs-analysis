#!/usr/bin/env python3
"""
Full genome-wide variant extraction for large-scale PCA analysis
Handles the 40GB VCF efficiently with chromosome-wise processing
"""

import subprocess
import time
import os
import argparse
from datetime import datetime

class GenomeWideExtractor:
    def __init__(self, vcf_path, output_dir="full_genome_extracts"):
        self.vcf_path = vcf_path
        self.output_dir = output_dir
        self.main_chromosomes = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
    def extract_chromosome_variants(self, chromosome, max_variants=5000, skip_existing=True):
        """Extract variants from a specific chromosome"""
        
        output_file = os.path.join(self.output_dir, f"{chromosome}_variants.vcf")
        
        # Skip if file exists and skip_existing is True
        if skip_existing and os.path.exists(output_file):
            print(f"  {chromosome}: Skipping (file exists)")
            return output_file
        
        print(f"  {chromosome}: Extracting up to {max_variants:,} variants...")
        start_time = time.time()
        
        try:
            # Get headers first
            cmd_headers = f"timeout 60 zcat {self.vcf_path} | grep '^#'"
            result = subprocess.run(cmd_headers, shell=True, capture_output=True, text=True)
            
            if result.returncode != 0:
                print(f"    {chromosome}: Failed to get headers")
                return None
            
            headers = result.stdout
            
            # Extract variants for this chromosome
            cmd_variants = f"timeout 300 zcat {self.vcf_path} | grep '^{chromosome}\\t' | head -{max_variants}"
            result = subprocess.run(cmd_variants, shell=True, capture_output=True, text=True)
            
            if result.returncode != 0:
                print(f"    {chromosome}: Failed to extract variants")
                return None
            
            variants = result.stdout
            variant_count = len(variants.splitlines()) if variants else 0
            
            if variant_count == 0:
                print(f"    {chromosome}: No variants found")
                return None
            
            # Write combined VCF
            with open(output_file, 'w') as f:
                f.write(headers)
                f.write(variants)
            
            elapsed = time.time() - start_time
            print(f"    {chromosome}: {variant_count:,} variants extracted in {elapsed:.1f}s")
            
            return output_file
            
        except Exception as e:
            print(f"    {chromosome}: Error - {e}")
            return None
    
    def extract_multiple_chromosomes(self, chromosomes, variants_per_chromosome=5000):
        """Extract variants from multiple chromosomes"""
        
        print(f"\\n{'='*60}")
        print(f"GENOME-WIDE VARIANT EXTRACTION")
        print(f"{'='*60}")
        print(f"Target chromosomes: {len(chromosomes)}")
        print(f"Variants per chromosome: {variants_per_chromosome:,}")
        print(f"Expected total variants: {len(chromosomes) * variants_per_chromosome:,}")
        
        successful_extractions = []
        total_variants = 0
        
        for i, chrom in enumerate(chromosomes, 1):
            print(f"\\nProcessing {i}/{len(chromosomes)}: {chrom}")
            
            output_file = self.extract_chromosome_variants(chrom, variants_per_chromosome)
            if output_file:
                successful_extractions.append((chrom, output_file))
                
                # Count variants in the file
                with open(output_file, 'r') as f:
                    variant_count = sum(1 for line in f if not line.startswith('#'))
                    total_variants += variant_count
        
        print(f"\\n{'='*60}")
        print(f"EXTRACTION SUMMARY")
        print(f"{'='*60}")
        print(f"Successful extractions: {len(successful_extractions)}/{len(chromosomes)}")
        print(f"Total variants extracted: {total_variants:,}")
        
        if successful_extractions:
            print(f"\\nSuccessful chromosomes:")
            for chrom, file_path in successful_extractions:
                print(f"  {chrom}: {file_path}")
        
        return successful_extractions, total_variants
    
    def merge_chromosome_files(self, chromosome_files, output_file="merged_genome_wide.vcf"):
        """Merge multiple chromosome VCF files into one"""
        
        output_path = os.path.join(self.output_dir, output_file)
        
        print(f"\\nMerging {len(chromosome_files)} chromosome files...")
        print(f"Output: {output_path}")
        
        if not chromosome_files:
            print("No files to merge!")
            return None
        
        try:
            # Use the first file as base (includes headers)
            first_file = chromosome_files[0][1]
            print(f"  Using {chromosome_files[0][0]} as header source")
            
            with open(output_path, 'w') as out_f:
                # Copy headers from first file
                with open(first_file, 'r') as f:
                    for line in f:
                        if line.startswith('#'):
                            out_f.write(line)
                        else:
                            out_f.write(line)  # Also copy variants from first file
                
                # Append variants from other files (skip headers)
                for chrom, file_path in chromosome_files[1:]:
                    print(f"  Appending variants from {chrom}")
                    with open(file_path, 'r') as f:
                        for line in f:
                            if not line.startswith('#'):
                                out_f.write(line)
            
            # Count final variants
            with open(output_path, 'r') as f:
                total_variants = sum(1 for line in f if not line.startswith('#'))
            
            print(f"  Merged file created: {total_variants:,} total variants")
            return output_path
            
        except Exception as e:
            print(f"  Merge failed: {e}")
            return None
    
    def extract_genome_wide_variants(self, target_chromosomes=None, variants_per_chromosome=5000, 
                                   merge_output=True):
        """Complete genome-wide extraction workflow"""
        
        if target_chromosomes is None:
            # Use main chromosomes (1-22, X, Y)
            target_chromosomes = self.main_chromosomes
        
        # Extract chromosome files
        chromosome_files, total_variants = self.extract_multiple_chromosomes(
            target_chromosomes, variants_per_chromosome
        )
        
        if not chromosome_files:
            print("\\n❌ No successful extractions!")
            return None
        
        merged_file = None
        if merge_output:
            merged_file = self.merge_chromosome_files(chromosome_files)
        
        print(f"\\n✅ Genome-wide extraction completed!")
        print(f"   Chromosomes processed: {len(chromosome_files)}")
        print(f"   Total variants: {total_variants:,}")
        if merged_file:
            print(f"   Merged file: {merged_file}")
        
        return {
            'chromosome_files': chromosome_files,
            'total_variants': total_variants,
            'merged_file': merged_file
        }

def main():
    parser = argparse.ArgumentParser(description='Extract genome-wide variants from large VCF')
    parser.add_argument('--vcf', default='/mnt/czi-sci-ai/intrinsic-variation-gene-ex/wgs/data/variants/all_vcfs/merged.all.biallelic.maf0.01.vcf.gz',
                       help='Path to input VCF file')
    parser.add_argument('--output-dir', default='full_genome_extracts',
                       help='Output directory for extracted files')
    parser.add_argument('--variants-per-chr', type=int, default=5000,
                       help='Maximum variants per chromosome')
    parser.add_argument('--chromosomes', type=str,
                       help='Comma-separated list of chromosomes (default: chr1-chr22,chrX,chrY)')
    parser.add_argument('--no-merge', action='store_true',
                       help='Skip merging chromosome files')
    
    args = parser.parse_args()
    
    # Parse chromosomes
    if args.chromosomes:
        target_chromosomes = [c.strip() for c in args.chromosomes.split(',')]
    else:
        target_chromosomes = None  # Use default
    
    # Create extractor
    extractor = GenomeWideExtractor(args.vcf, args.output_dir)
    
    # Run extraction
    results = extractor.extract_genome_wide_variants(
        target_chromosomes=target_chromosomes,
        variants_per_chromosome=args.variants_per_chr,
        merge_output=not args.no_merge
    )
    
    return results

if __name__ == "__main__":
    results = main()