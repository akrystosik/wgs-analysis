#!/usr/bin/env python3
"""
Quick multi-chromosome extractor with early termination for faster processing
"""

import subprocess
import time
import os

def quick_multi_chromosome_extract(vcf_path, chromosomes, variants_per_chr=2000, max_time_minutes=10):
    """Extract variants from multiple chromosomes with time limit"""
    
    print("="*60)
    print("QUICK MULTI-CHROMOSOME EXTRACTION")
    print("="*60)
    print(f"Target chromosomes: {', '.join(chromosomes)}")
    print(f"Variants per chromosome: {variants_per_chr:,}")
    print(f"Time limit: {max_time_minutes} minutes")
    
    chr_counts = {chrom: 0 for chrom in chromosomes}
    header_lines = []
    collected_variants = []
    
    start_time = time.time()
    max_time_seconds = max_time_minutes * 60
    
    # Stream through VCF with time limit
    cmd = f"timeout {max_time_seconds + 30} zcat {vcf_path}"
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, 
                             stderr=subprocess.PIPE, text=True, bufsize=8192)
    
    lines_processed = 0
    variants_found = 0
    
    try:
        for line in process.stdout:
            current_time = time.time()
            elapsed = current_time - start_time
            
            # Time limit check
            if elapsed > max_time_seconds:
                print(f"\\n⏰ Time limit reached ({max_time_minutes} minutes)")
                break
            
            lines_processed += 1
            
            # Progress reporting every 500K lines
            if lines_processed % 500000 == 0:
                print(f"  {elapsed:.1f}s: {lines_processed/1000000:.1f}M lines, {variants_found:,} variants")
                for chrom in chromosomes:
                    if chr_counts[chrom] > 0:
                        print(f"    {chrom}: {chr_counts[chrom]:,}")
            
            if line.startswith('##'):
                if len(header_lines) < 100:  # Limit header collection
                    header_lines.append(line)
            
            elif line.startswith('#CHROM'):
                header_lines.append(line)
            
            else:
                # Variant line - quick chromosome check
                if line[:5] in [chrom[:5] for chrom in chromosomes]:  # Quick prefix check
                    chrom = line.split('\\t', 1)[0]
                    
                    if chrom in chromosomes and chr_counts[chrom] < variants_per_chr:
                        collected_variants.append(line)
                        chr_counts[chrom] += 1
                        variants_found += 1
                        
                        # Early termination if we have enough for all chromosomes
                        if all(chr_counts[c] >= variants_per_chr for c in chromosomes):
                            print(f"\\n🎯 All target chromosomes reached {variants_per_chr:,} variants!")
                            break
    
    except Exception as e:
        print(f"\\nProcessing error: {e}")
    
    finally:
        process.terminate()
        process.wait()
    
    elapsed = time.time() - start_time
    
    # Write output
    output_file = f"quick_multi_chr_{len(chromosomes)}chr_{variants_per_chr}var.vcf"
    print(f"\\nWriting results to {output_file}...")
    
    with open(output_file, 'w') as f:
        f.writelines(header_lines)
        f.writelines(collected_variants)
    
    # Report results
    print(f"\\n{'='*40}")
    print(f"QUICK EXTRACTION RESULTS")
    print(f"{'='*40}")
    print(f"Processing time: {elapsed:.1f}s ({elapsed/60:.1f} min)")
    print(f"Lines processed: {lines_processed:,}")
    print(f"Total variants: {variants_found:,}")
    
    print(f"\\nChromosome results:")
    successful_chromosomes = 0
    for chrom in chromosomes:
        count = chr_counts[chrom]
        if count > 0:
            successful_chromosomes += 1
            status = "✅" if count >= variants_per_chr else "⚠️"
            print(f"  {chrom}: {count:,} variants {status}")
        else:
            print(f"  {chrom}: 0 variants ❌")
    
    file_size_mb = os.path.getsize(output_file) / 1024 / 1024 if os.path.exists(output_file) else 0
    print(f"\\nOutput: {output_file} ({file_size_mb:.1f} MB)")
    print(f"Successful chromosomes: {successful_chromosomes}/{len(chromosomes)}")
    
    return {
        'output_file': output_file,
        'total_variants': variants_found,
        'chromosome_counts': dict(chr_counts),
        'successful_chromosomes': successful_chromosomes,
        'processing_time': elapsed
    }

if __name__ == "__main__":
    vcf_path = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/wgs/data/variants/all_vcfs/merged.all.biallelic.maf0.01.vcf.gz"
    
    # Try extracting from 5 main chromosomes
    chromosomes = ['chr1', 'chr2', 'chr3', 'chr22', 'chrX']
    
    results = quick_multi_chromosome_extract(
        vcf_path, 
        chromosomes, 
        variants_per_chr=3000,
        max_time_minutes=8
    )
    
    print(f"\\n✅ Quick extraction completed!")
    if results['successful_chromosomes'] >= 2:
        print(f"Ready for multi-chromosome PCA analysis!")
    else:
        print(f"Consider using larger time limit or fewer chromosomes.")