#!/usr/bin/env python3

import pandas as pd
import subprocess
import os
from pathlib import Path

def check_chromosome_format(file_path, delimiter='\t'):
    """Check chromosome naming format in a file (chr1 vs 1)"""
    try:
        # Read first few lines to detect format
        df = pd.read_csv(file_path, delimiter=delimiter, nrows=1000)
        
        # Get chromosome column (try common names)
        chrom_col = None
        for col in ['Chromosome', 'chr', '#CHROM', 'chrom']:
            if col in df.columns:
                chrom_col = col
                break
        
        if chrom_col is None:
            # If no standard name found, assume first column
            chrom_col = df.columns[0]
        
        # Get unique chromosome names
        chroms = df[chrom_col].unique()
        
        # Check format
        has_chr_prefix = any(str(c).startswith('chr') for c in chroms)
        example_chroms = sorted(list(chroms))[:5]
        
        return {
            'has_chr_prefix': has_chr_prefix,
            'example_chroms': example_chroms,
            'format': 'chr1' if has_chr_prefix else '1'
        }
    except Exception as e:
        return {'error': str(e)}

def download_and_prepare_reference(build='GRCh38', output_dir='data/references'):
    """Download and prepare reference data"""
    
    # Create output directory
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    
    print(f"Preparing reference data in {output_dir}...")
    
    # Download GENCODE GTF
    gtf_url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.basic.annotation.gtf.gz"
    gtf_path = f"{output_dir}/gencode.v44.gtf.gz"
    
    if not os.path.exists(gtf_path):
        print(f"Downloading GENCODE GTF...")
        subprocess.run(['wget', gtf_url, '-O', gtf_path], check=True)
    
    # Extract coding regions - create both chr1 and 1 format versions
    coding_bed = f"{output_dir}/coding_regions.bed"
    coding_bed_nochr = f"{output_dir}/coding_regions.nochr.bed"
    
    print("Extracting coding regions...")
    cmd = f"""zcat {gtf_path} | awk '$3=="CDS"' | cut -f1,4,5 | sort -k1,1 -k2,2n > {coding_bed}"""
    subprocess.run(cmd, shell=True, check=True)
    
    # Create version without 'chr' prefix
    cmd = f"cat {coding_bed} | sed 's/^chr//' > {coding_bed_nochr}"
    subprocess.run(cmd, shell=True, check=True)
    
    return coding_bed, coding_bed_nochr

def verify_coordinates(variants_file, reference_file):
    """Compare coordinate systems between variant and reference files"""
    
    var_format = check_chromosome_format(variants_file)
    ref_format = check_chromosome_format(reference_file)
    
    print("\nChromosome format comparison:")
    print(f"Variants file: {var_format['format']}")
    print(f"Reference file: {ref_format['format']}")
    print("\nExample chromosomes in variants:", var_format['example_chroms'])
    print("Example chromosomes in reference:", ref_format['example_chroms'])
    
    return var_format['format'] == ref_format['format']

def main():
    # Define paths
    variants_path = "data/variants/tumor_only_analysis/potential_somatic.bed"
    
    # Download and prepare reference data
    coding_bed, coding_bed_nochr = download_and_prepare_reference()
    
    # Check variant file format
    print("\nAnalyzing variant file format...")
    var_format = check_chromosome_format(variants_path)
    print(f"Variant file uses {var_format['format']} format")
    
    # Choose appropriate reference file
    reference_to_use = coding_bed if var_format['has_chr_prefix'] else coding_bed_nochr
    print(f"\nUsing reference file: {reference_to_use}")
    
    # Verify coordinate systems match
    print("\nVerifying coordinate systems...")
    formats_match = verify_coordinates(variants_path, reference_to_use)
    
    if formats_match:
        print("\nFormats match! You can proceed with:")
        print(f"bedtools intersect -a {variants_path} -b {reference_to_use} > potential_somatic.coding.bed")
    else:
        print("\nWarning: Format mismatch. You need to convert chromosome formats to match.")
        print("You can convert using:")
        if var_format['has_chr_prefix']:
            print("sed 's/^chr//' input.bed > output.bed  # Remove chr prefix")
        else:
            print("sed 's/^/chr/' input.bed > output.bed  # Add chr prefix")

if __name__ == "__main__":
    main()