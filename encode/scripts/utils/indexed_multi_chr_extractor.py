#!/usr/bin/env python3
"""
Multi-chromosome variant extraction using tabix index for efficient access
"""

import subprocess
import time
import os
import tempfile

def extract_indexed_multi_chromosome(vcf_path, chromosomes, variants_per_chr=5000, output_file="indexed_multi_chr.vcf"):
    """
    Extract variants from multiple chromosomes using tabix index
    Much faster than streaming through the entire file
    """
    
    print("="*60)
    print("INDEXED MULTI-CHROMOSOME EXTRACTION")
    print("="*60)
    print(f"VCF file: {vcf_path}")
    print(f"Target chromosomes: {', '.join(chromosomes)}")
    print(f"Variants per chromosome: {variants_per_chr:,}")
    
    # Check if index exists
    index_file = f"{vcf_path}.tbi"
    if not os.path.exists(index_file):
        print(f"❌ Index file not found: {index_file}")
        return None
    
    start_time = time.time()
    temp_files = []
    chr_results = {}
    
    try:
        # Extract variants from each chromosome
        for i, chrom in enumerate(chromosomes, 1):
            print(f"\\n{i}/{len(chromosomes)}: Extracting from {chrom}")
            
            temp_file = tempfile.NamedTemporaryFile(mode='w', suffix=f'_{chrom}.vcf', delete=False)
            temp_files.append(temp_file.name)
            
            # Use bcftools to extract specific chromosome
            cmd = [
                'bcftools', 'view', vcf_path,
                '--regions', chrom,
                '--max-alleles', '2',
                '--min-alleles', '2'
            ]
            
            print(f"  Running: {' '.join(cmd)}")
            
            # Run bcftools and capture output
            process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            
            variant_count = 0
            header_written = False
            
            with open(temp_file.name, 'w') as f:
                for line in process.stdout:
                    if line.startswith('#'):
                        # Write headers for first file or if we haven't written headers yet
                        if i == 1 or not header_written:
                            f.write(line)
                            header_written = True
                    else:
                        # Variant line
                        if variant_count < variants_per_chr:
                            f.write(line)
                            variant_count += 1
                        else:
                            break
            
            process.wait()
            
            if process.returncode == 0:
                chr_results[chrom] = variant_count
                print(f"  ✅ {chrom}: {variant_count:,} variants extracted")
            else:
                print(f"  ❌ {chrom}: Failed to extract")
                chr_results[chrom] = 0
        
        # Merge all chromosome files
        print(f"\\nMerging {len(temp_files)} chromosome files...")
        
        with open(output_file, 'w') as out_f:
            headers_written = False
            total_variants = 0
            
            for i, (temp_file, chrom) in enumerate(zip(temp_files, chromosomes)):
                if chr_results.get(chrom, 0) > 0:
                    print(f"  Adding {chrom} variants...")
                    
                    with open(temp_file, 'r') as f:
                        for line in f:
                            if line.startswith('#'):
                                # Only write headers from first file
                                if not headers_written:
                                    out_f.write(line)
                            else:
                                # Always write variant lines
                                out_f.write(line)
                                total_variants += 1
                    
                    headers_written = True
        
        elapsed = time.time() - start_time
        
        # Report results
        print(f"\\n{'='*50}")
        print(f"INDEXED EXTRACTION RESULTS")
        print(f"{'='*50}")
        print(f"Processing time: {elapsed:.1f} seconds")
        print(f"Total variants: {total_variants:,}")
        
        successful_chrs = 0
        print(f"\\nChromosome results:")
        for chrom in chromosomes:
            count = chr_results.get(chrom, 0)
            if count > 0:
                successful_chrs += 1
                status = "✅" if count >= variants_per_chr * 0.8 else "⚠️"  # Allow 80% target
                print(f"  {chrom}: {count:,} variants {status}")
            else:
                print(f"  {chrom}: 0 variants ❌")
        
        file_size = os.path.getsize(output_file) / 1024 / 1024
        print(f"\\nOutput file: {output_file} ({file_size:.1f} MB)")
        print(f"Successful chromosomes: {successful_chrs}/{len(chromosomes)}")
        
        return {
            'output_file': output_file,
            'total_variants': total_variants,
            'chromosome_results': chr_results,
            'successful_chromosomes': successful_chrs,
            'processing_time': elapsed
        }
    
    finally:
        # Clean up temp files
        for temp_file in temp_files:
            if os.path.exists(temp_file):
                os.remove(temp_file)

def main():
    vcf_path = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/wgs/data/variants/all_vcfs/merged.all.biallelic.maf0.01.vcf.gz"
    
    # Start with key chromosomes for genome-wide representation
    chromosomes = ['chr1', 'chr2', 'chr22']  # Reduced for faster processing
    
    results = extract_indexed_multi_chromosome(
        vcf_path,
        chromosomes,
        variants_per_chr=6000,  # Higher density for genome-wide coverage
        output_file="indexed_genome_wide_3chr.vcf"
    )
    
    if results and results['successful_chromosomes'] >= 2:
        print(f"\\n🎉 Success! Ready for multi-chromosome PCA analysis")
        print(f"Use: python3 scripts/analysis/enhanced_pca_analysis.py {results['output_file']} results/genome_wide")
    else:
        print(f"\\n⚠️ Limited success. Consider using fewer chromosomes or checking data availability.")
    
    return results

if __name__ == "__main__":
    main()