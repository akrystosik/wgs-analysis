#!/usr/bin/env python3
"""
Optimized genome-wide variant extraction using streaming approach
"""

import subprocess
import time
import os
import argparse
from collections import defaultdict

def extract_multi_chromosome_variants(vcf_path, chromosomes, variants_per_chr=3000, output_file="multi_chr_sample.vcf"):
    """
    Extract variants from multiple chromosomes in a single pass through the VCF
    Much more efficient than separate chromosome extractions
    """
    
    print("="*60)
    print("OPTIMIZED MULTI-CHROMOSOME EXTRACTION")
    print("="*60)
    print(f"Target chromosomes: {', '.join(chromosomes)}")
    print(f"Variants per chromosome: {variants_per_chr:,}")
    print(f"Expected total: {len(chromosomes) * variants_per_chr:,} variants")
    
    # Initialize counters
    chr_counts = {chrom: 0 for chrom in chromosomes}
    header_lines = []
    collected_variants = []
    
    start_time = time.time()
    
    # Stream through VCF file once
    print("\\nStreaming through 40GB VCF file...")
    cmd = f"zcat {vcf_path}"
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, 
                             stderr=subprocess.PIPE, text=True, bufsize=1)
    
    lines_processed = 0
    variants_found = 0
    
    try:
        for line in process.stdout:
            lines_processed += 1
            
            # Progress reporting
            if lines_processed % 1000000 == 0:
                elapsed = time.time() - start_time
                print(f"  Processed {lines_processed/1000000:.1f}M lines in {elapsed:.1f}s")
                print(f"    Variants collected: {variants_found:,}")
                for chrom in chromosomes:
                    print(f"      {chrom}: {chr_counts[chrom]:,}")
            
            if line.startswith('##'):
                # Header lines
                if not header_lines or len(header_lines) < 200:  # Limit header collection
                    header_lines.append(line)
            
            elif line.startswith('#CHROM'):
                # Column header
                header_lines.append(line)
            
            else:
                # Variant line
                parts = line.split('\\t', 2)  # Only split first two columns for efficiency
                if len(parts) >= 1:
                    chrom = parts[0]
                    
                    # Check if this chromosome is in our target list
                    if chrom in chromosomes and chr_counts[chrom] < variants_per_chr:
                        collected_variants.append(line)
                        chr_counts[chrom] += 1
                        variants_found += 1
                        
                        # Check if we have enough variants total
                        if variants_found >= len(chromosomes) * variants_per_chr:
                            print(f"\\n  Target variant count reached!")
                            break
                        
                        # Check if we have enough for all chromosomes
                        if all(chr_counts[c] >= variants_per_chr for c in chromosomes):
                            print(f"\\n  All chromosomes have sufficient variants!")
                            break
    
    except Exception as e:
        print(f"\\nError during processing: {e}")
    
    finally:
        process.terminate()
        process.wait()
    
    # Write output file
    print(f"\\nWriting output file: {output_file}")
    with open(output_file, 'w') as f:
        # Write headers
        f.writelines(header_lines)
        # Write variants
        f.writelines(collected_variants)
    
    elapsed = time.time() - start_time
    
    # Report results
    print(f"\\n{'='*60}")
    print(f"EXTRACTION RESULTS")
    print(f"{'='*60}")
    print(f"Total time: {elapsed:.1f} seconds ({elapsed/60:.1f} minutes)")
    print(f"Lines processed: {lines_processed:,}")
    print(f"Total variants collected: {variants_found:,}")
    print(f"\\nVariants by chromosome:")
    
    total_collected = 0
    for chrom in chromosomes:
        count = chr_counts[chrom]
        total_collected += count
        status = "✅ Complete" if count >= variants_per_chr else f"⚠️  Partial ({count}/{variants_per_chr})"
        print(f"  {chrom}: {count:,} variants {status}")
    
    print(f"\\nOutput file: {output_file}")
    print(f"File size: {os.path.getsize(output_file)/1024/1024:.1f} MB")
    
    return {
        'output_file': output_file,
        'total_variants': total_collected,
        'chromosome_counts': dict(chr_counts),
        'processing_time': elapsed
    }

def main():
    parser = argparse.ArgumentParser(description='Optimized multi-chromosome variant extraction')
    parser.add_argument('--vcf', 
                       default='/mnt/czi-sci-ai/intrinsic-variation-gene-ex/wgs/data/variants/all_vcfs/merged.all.biallelic.maf0.01.vcf.gz',
                       help='Path to input VCF file')
    parser.add_argument('--chromosomes', type=str, default='chr1,chr2,chr22',
                       help='Comma-separated list of chromosomes')
    parser.add_argument('--variants-per-chr', type=int, default=3000,
                       help='Variants per chromosome')
    parser.add_argument('--output', type=str, default='optimized_multi_chr_sample.vcf',
                       help='Output VCF file')
    
    args = parser.parse_args()
    
    chromosomes = [c.strip() for c in args.chromosomes.split(',')]
    
    results = extract_multi_chromosome_variants(
        args.vcf, 
        chromosomes, 
        args.variants_per_chr, 
        args.output
    )
    
    print(f"\\n✅ Optimized extraction completed!")
    return results

if __name__ == "__main__":
    main()