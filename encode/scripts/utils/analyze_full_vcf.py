#!/usr/bin/env python3
"""
Analyze the full 40GB VCF structure efficiently without full indexing
"""

import subprocess
import time
import os
from collections import defaultdict

def analyze_vcf_structure(vcf_path, sample_lines=100000):
    """Analyze VCF structure by sampling lines efficiently"""
    
    print("="*60)
    print("ANALYZING 40GB VCF STRUCTURE")
    print("="*60)
    
    # Get file info
    file_size = os.path.getsize(vcf_path) / (1024**3)  # GB
    print(f"File size: {file_size:.1f} GB")
    
    # Stream through file and analyze structure
    print(f"Sampling first {sample_lines:,} variant lines...")
    
    cmd = f"timeout 300 zcat {vcf_path}"
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, 
                             stderr=subprocess.PIPE, text=True, bufsize=1)
    
    header_lines = 0
    variant_lines = 0
    chromosome_counts = defaultdict(int)
    sample_count = 0
    first_pos = None
    last_pos = None
    current_chr = None
    
    try:
        for line_num, line in enumerate(process.stdout):
            if line.startswith('##'):
                header_lines += 1
                # Extract contig info
                if line.startswith('##contig='):
                    chrom = line.split('ID=')[1].split(',')[0] if 'ID=' in line else 'unknown'
                    print(f"  Found chromosome: {chrom}")
            
            elif line.startswith('#CHROM'):
                # Count samples
                columns = line.strip().split('\t')
                sample_count = len(columns) - 9  # First 9 are VCF standard columns
                print(f"  Samples: {sample_count}")
                print(f"  Headers complete: {header_lines} lines")
            
            else:
                # Variant line
                variant_lines += 1
                if variant_lines <= sample_lines:
                    parts = line.split('\t')
                    if len(parts) >= 2:
                        chrom = parts[0]
                        pos = int(parts[1])
                        
                        # Track chromosome distribution
                        chromosome_counts[chrom] += 1
                        
                        # Track positions
                        if first_pos is None:
                            first_pos = pos
                            current_chr = chrom
                        last_pos = pos
                
                # Progress reporting
                if variant_lines % 50000 == 0:
                    print(f"    Processed {variant_lines:,} variants...")
                
                # Stop after sampling enough
                if variant_lines >= sample_lines:
                    break
                    
    except Exception as e:
        print(f"Processing interrupted: {e}")
    finally:
        process.terminate()
        process.wait()
    
    print(f"\nAnalysis Results:")
    print(f"  Header lines: {header_lines:,}")
    print(f"  Variant lines sampled: {variant_lines:,}")
    print(f"  Samples: {sample_count:,}")
    
    print(f"\nChromosome distribution (first {sample_lines:,} variants):")
    for chrom, count in sorted(chromosome_counts.items()):
        pct = count / variant_lines * 100 if variant_lines > 0 else 0
        print(f"  {chrom}: {count:,} variants ({pct:.1f}%)")
    
    if first_pos and last_pos:
        print(f"\nPosition range in {current_chr}: {first_pos:,} - {last_pos:,}")
    
    # Estimate total variants
    if variant_lines > 0:
        # This is a rough estimate based on sampling
        print(f"\nEstimations:")
        print(f"  Likely contains multiple chromosomes: {len(chromosome_counts) > 1}")
        if len(chromosome_counts) == 1:
            print(f"  All sampled variants from {list(chromosome_counts.keys())[0]}")
            print(f"  May need deeper sampling to find other chromosomes")
    
    return {
        'file_size_gb': file_size,
        'header_lines': header_lines,
        'variant_lines_sampled': variant_lines,
        'sample_count': sample_count,
        'chromosome_counts': dict(chromosome_counts),
        'first_pos': first_pos,
        'last_pos': last_pos
    }

def estimate_total_variants(vcf_path):
    """Estimate total variants using multiple sampling points"""
    
    print("\nEstimating total variant count...")
    
    # Try to sample from different points in the file
    file_size = os.path.getsize(vcf_path)
    
    # Sample from beginning (we already know this is chr1)
    print("  Sampling from file beginning...")
    
    # Get a rough line count estimate
    print("  Estimating total lines...")
    cmd = f"timeout 120 zcat {vcf_path} | head -1000000 | wc -l"
    try:
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=150)
        if result.returncode == 0:
            first_million_lines = int(result.stdout.strip())
            print(f"    First 1M lines processed successfully")
            
            # Very rough estimate: assume uniform compression
            estimated_total_lines = (first_million_lines * file_size) / (file_size * 0.025)  # Very rough
            print(f"    Rough estimate: {estimated_total_lines/1000000:.1f}M total lines")
    
    except Exception as e:
        print(f"    Line counting failed: {e}")
    
    return None

if __name__ == "__main__":
    vcf_path = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/wgs/data/variants/all_vcfs/merged.all.biallelic.maf0.01.vcf.gz"
    
    # Analyze structure
    results = analyze_vcf_structure(vcf_path, sample_lines=100000)
    
    # Try to estimate total size
    estimate_total_variants(vcf_path)
    
    print(f"\n✅ VCF structure analysis completed!")
    print(f"Ready to design efficient extraction strategy.")