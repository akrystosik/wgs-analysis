#!/usr/bin/env python3
"""
Extract a minimal sample from the real VCF for PCA analysis
"""

import subprocess
import tempfile
import os
import time

def extract_minimal_vcf_sample(vcf_path, output_path, max_variants=200, chromosomes=None):
    """Extract minimal VCF sample efficiently using bcftools"""
    print(f"Extracting {max_variants} variants from real VCF...")
    if chromosomes:
        print(f"Target chromosomes: {', '.join(chromosomes)}")
    
    start_time = time.time()
    
    # Use streaming approach for efficient extraction
    if chromosomes:
        # Multi-chromosome extraction using streaming
        variants_per_chr = max_variants // len(chromosomes)
        print(f"  Extracting ~{variants_per_chr} variants per chromosome")
        
        collected_variants = []
        header_lines = []
        
        # Use more efficient grep-based approach
        print("  Using targeted chromosome extraction...")
        
        # Get headers first (including #CHROM line)
        cmd_headers = f"timeout 30 zcat {vcf_path} | grep '^#'"
        result = subprocess.run(cmd_headers, shell=True, capture_output=True, text=True)
        if result.returncode == 0:
            header_lines = result.stdout.splitlines(keepends=True)
        
        chr_counts = {chrom: 0 for chrom in chromosomes}
        
        # Extract variants for each chromosome separately
        for chrom in chromosomes:
            print(f"    Processing chromosome {chrom}...")
            # Use grep to filter lines starting with the chromosome
            cmd_vars = f"timeout 60 zcat {vcf_path} | grep '^{chrom}\t' | head -{variants_per_chr}"
            result = subprocess.run(cmd_vars, shell=True, capture_output=True, text=True)
            
            if result.returncode == 0 and result.stdout:
                lines = result.stdout.splitlines(keepends=True)
                collected_variants.extend(lines)
                chr_counts[chrom] = len(lines)
                print(f"      Found {len(lines)} variants")
            else:
                print(f"      No variants found for {chrom}")
            
            # Break if we have enough total variants
            if sum(chr_counts.values()) >= max_variants:
                break
        
        # Write collected data
        print("  Writing collected variants...")
        with open(output_path, 'w') as f:
            f.writelines(header_lines)
            f.writelines(collected_variants)
        
        print(f"  Collected variants by chromosome:")
        for chrom, count in chr_counts.items():
            print(f"    {chrom}: {count} variants")
    else:
        # Single-chromosome extraction using original streaming method
        print("  Using streaming extraction for all chromosomes...")
        with open(output_path, 'w') as out_file:
            # Get all headers (including #CHROM)
            cmd1 = f"timeout 30 zcat {vcf_path} | grep '^#'"
            result1 = subprocess.run(cmd1, shell=True, capture_output=True, text=True)
            if result1.returncode == 0:
                out_file.write(result1.stdout)
            
            # Get variants
            cmd3 = f"timeout 60 zcat {vcf_path} | grep -v '^#' | head -{max_variants}"
            result3 = subprocess.run(cmd3, shell=True, capture_output=True, text=True)
            if result3.returncode == 0:
                out_file.write(result3.stdout)
    
    elapsed = time.time() - start_time
    
    # Check if we got data
    if os.path.exists(output_path):
        line_count = sum(1 for line in open(output_path))
        print(f"✓ Extracted VCF with {line_count} lines in {elapsed:.1f} seconds")
        return True
    else:
        print("❌ Failed to extract VCF sample")
        return False

def main():
    import argparse
    
    parser = argparse.ArgumentParser(description='Extract VCF sample for PCA analysis')
    parser.add_argument('--max-variants', type=int, default=5000, 
                       help='Maximum number of variants to extract (default: 5000)')
    parser.add_argument('--chromosomes', type=str, 
                       help='Comma-separated list of chromosomes (e.g., chr1,chr2,chr22)')
    parser.add_argument('--output', type=str, default='real_vcf_sample.vcf',
                       help='Output file path (default: real_vcf_sample.vcf)')
    
    args = parser.parse_args()
    
    vcf_path = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/wgs/data/variants/all_vcfs/merged.all.biallelic.maf0.01.vcf.gz"
    
    chromosomes = None
    if args.chromosomes:
        chromosomes = [c.strip() for c in args.chromosomes.split(',')]
    
    print("="*60)
    print("EXTRACTING REAL VCF SAMPLE FOR PCA ANALYSIS")
    print("="*60)
    print(f"Max variants: {args.max_variants}")
    if chromosomes:
        print(f"Chromosomes: {', '.join(chromosomes)}")
    
    if extract_minimal_vcf_sample(vcf_path, args.output, max_variants=args.max_variants, chromosomes=chromosomes):
        print(f"\n✅ Successfully created: {args.output}")
        print("Now ready to run PCA analysis on real data!")
        
        # Show what we extracted
        print("\nSample preview:")
        with open(args.output, 'r') as f:
            lines = f.readlines()
            header_lines = [l for l in lines if l.startswith('#')]
            data_lines = [l for l in lines if not l.startswith('#')]
            
            print(f"  Headers: {len(header_lines)} lines")
            print(f"  Variants: {len(data_lines)} lines")
            
            if len(data_lines) > 0:
                # Parse column header to count samples
                chrom_line = [l for l in header_lines if l.startswith('#CHROM')]
                if chrom_line:
                    columns = chrom_line[0].strip().split('\t')
                    sample_count = len(columns) - 9  # First 9 are VCF standard columns
                    print(f"  Samples: {sample_count}")
                    
                    # Show first few sample names
                    sample_names = columns[9:19]  # First 10 samples
                    print(f"  First samples: {', '.join(sample_names)}")
                    
                # Show chromosome distribution
                chr_counts = {}
                for line in data_lines:
                    chrom = line.split('\t')[0]
                    chr_counts[chrom] = chr_counts.get(chrom, 0) + 1
                
                print(f"  Chromosome distribution:")
                for chrom, count in sorted(chr_counts.items()):
                    print(f"    {chrom}: {count} variants")
    else:
        print("❌ Failed to extract VCF sample")

if __name__ == "__main__":
    main()