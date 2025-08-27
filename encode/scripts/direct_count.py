#!/usr/bin/env python
import os
import subprocess
import pandas as pd
from tqdm import tqdm
import multiprocessing
import argparse
import logging
from datetime import datetime
import random
import time
import re

# Set up logging
log_dir = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/wgs/data/qc/variants/logs"
os.makedirs(log_dir, exist_ok=True)
log_file = os.path.join(log_dir, f"direct_count_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")
logging.basicConfig(
    filename=log_file,
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# Define paths
gtex_vcf = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex-3/project_gene_regulation/data/GTEx/WGS/phg001796.v1.GTEx_v9_WGS_953.genotype-calls-vcf.c1/GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv.vcf"
output_dir = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/wgs/data/qc/variants/direct_count"
summary_file = os.path.join(output_dir, "gtex_sample_stats.tsv")

# Create output directory
os.makedirs(output_dir, exist_ok=True)

def get_gtex_samples():
    """Get list of all samples in the GTEx VCF file"""
    print("Fetching sample IDs from GTEx VCF...")
    logging.info("Fetching sample IDs from GTEx VCF...")
    
    # Use grep instead of bcftools for uncompressed VCF
    cmd = f"grep -m 1 '#CHROM' {gtex_vcf} | cut -f10-"
    try:
        result = subprocess.check_output(cmd, shell=True, text=True)
        samples = result.strip().split('\t')
        print(f"Found {len(samples)} GTEx samples")
        logging.info(f"Found {len(samples)} GTEx samples")
        return samples
    except subprocess.CalledProcessError as e:
        print(f"Error fetching GTEx samples: {e}")
        logging.error(f"Error fetching GTEx samples: {e}")
        return []

def count_variants_direct(sample_id, max_time=3600):
    """Count variants for a sample using grep and awk"""
    output_file = os.path.join(output_dir, f"{sample_id}.counts.txt")
    
    # Check if output file already exists
    if os.path.exists(output_file):
        print(f"Counts for {sample_id} already exist, skipping")
        logging.info(f"Counts for {sample_id} already exist, skipping")
        return True, "Already processed"
    
    print(f"Counting variants for {sample_id} using direct counting...")
    logging.info(f"Counting variants for {sample_id} using direct counting...")
    
    # Get column index for this sample (1-based for awk)
    sample_idx_cmd = f"grep -m 1 '#CHROM' {gtex_vcf} | tr '\\t' '\\n' | grep -n '^{sample_id}$' | cut -d: -f1"
    try:
        sample_idx = int(subprocess.check_output(sample_idx_cmd, shell=True, text=True).strip())
    except:
        error_msg = f"Could not find sample {sample_id} in VCF header"
        print(error_msg)
        logging.error(error_msg)
        return False, error_msg
    
    # Count SNPs (single nucleotide variants where sample has non-reference genotype)
    snp_cmd = f"""grep -v "^#" {gtex_vcf} | awk 'length($4)==1 && length($5)==1 && $10 !~ /0\\/0/ {{ count++ }} END {{ print "SNPs", count }}'"""
    
    # Count indels (length of ref or alt > 1 and sample has non-reference genotype)
    indel_cmd = f"""grep -v "^#" {gtex_vcf} | awk '(length($4)>1 || length($5)>1) && $10 !~ /0\\/0/ {{ count++ }} END {{ print "Indels", count }}'"""
    
    # Combine commands
    cmd = f"echo 'Type\tCount' > {output_file} && {snp_cmd} >> {output_file} && {indel_cmd} >> {output_file}"
    
    try:
        # This might be very slow on a large VCF, consider splitting by chromosome
        process = subprocess.Popen(cmd, shell=True)
        
        # Monitor with timeout
        start_time = time.time()
        while process.poll() is None:
            if time.time() - start_time > max_time:
                process.kill()
                error_msg = f"Timeout processing {sample_id} after {max_time} seconds"
                print(error_msg)
                logging.error(error_msg)
                return False, error_msg
            time.sleep(10)  # Check every 10 seconds
        
        if process.returncode != 0:
            error_msg = f"Error processing {sample_id}, return code: {process.returncode}"
            print(error_msg)
            logging.error(error_msg)
            return False, error_msg
        
        print(f"Successfully processed sample {sample_id}")
        logging.info(f"Successfully processed sample {sample_id}")
        return True, "Success"
    
    except Exception as e:
        error_msg = f"Error processing sample {sample_id}: {e}"
        print(error_msg)
        logging.error(error_msg)
        return False, error_msg

def process_by_chromosome(sample_id, chroms=None, max_time=1800):
    """Process a sample chromosome by chromosome using grep and awk"""
    if chroms is None:
        # Use numbers instead of "chr" prefix
        chroms = [str(i) for i in range(1, 23)] + ["X", "Y"]
    
    output_dir_by_chrom = os.path.join(output_dir, "by_chromosome", sample_id)
    os.makedirs(output_dir_by_chrom, exist_ok=True)
    
    # Final output file
    output_file = os.path.join(output_dir, f"{sample_id}.counts.txt")
    
    # Check if final output already exists
    if os.path.exists(output_file):
        print(f"Counts for {sample_id} already exist, skipping")
        logging.info(f"Counts for {sample_id} already exist, skipping")
        return {
            'sample_id': sample_id,
            'num_snps': 0,
            'num_indels': 0,
            'total_variants': 0
        }
    
    # Get column index for this sample (1-based for awk)
    sample_idx_cmd = f"grep -m 1 '#CHROM' {gtex_vcf} | tr '\\t' '\\n' | grep -n '^{sample_id}$' | cut -d: -f1"
    try:
        sample_idx = int(subprocess.check_output(sample_idx_cmd, shell=True, text=True).strip())
        print(f"Sample {sample_id} is at column {sample_idx}")
        logging.info(f"Sample {sample_id} is at column {sample_idx}")
    except:
        error_msg = f"Could not find sample {sample_id} in VCF header"
        print(error_msg)
        logging.error(error_msg)
        return {
            'sample_id': sample_id,
            'num_snps': 0,
            'num_indels': 0,
            'total_variants': 0
        }
    
    total_snps = 0
    total_indels = 0
    
    for chrom in tqdm(chroms, desc=f"Processing {sample_id} by chromosome"):
        chrom_file = os.path.join(output_dir_by_chrom, f"{chrom}.counts.txt")
        
        # Skip if already processed
        if os.path.exists(chrom_file):
            # Read previous results
            with open(chrom_file, 'r') as f:
                lines = f.readlines()
                for line in lines[1:]:  # Skip header
                    parts = line.strip().split('\t')
                    if parts[0] == "SNPs":
                        total_snps += int(parts[1])
                    elif parts[0] == "Indels":
                        total_indels += int(parts[1])
            print(f"Counts for {sample_id} {chrom} already exist, loaded previous results")
            logging.info(f"Counts for {sample_id} {chrom} already exist, loaded previous results")
            continue
        
        # Count SNPs for this chromosome and sample
        # The awk script filters for: 
        # 1. Lines with this chromosome 
        # 2. Where ref and alt are single nucleotides (SNPs)
        # 3. Where the genotype for this sample is not 0/0 (reference)
        snp_cmd = f"""grep -v "^#" {gtex_vcf} | grep "^{chrom}\\s" | awk 'length($4)==1 && length($5)==1 && ${sample_idx} !~ /0\\/0/ && ${sample_idx} !~ /\\.\\/\\./ {{ count++ }} END {{ print "SNPs\\t" count }}'"""
        
        # Count indels for this chromosome and sample
        indel_cmd = f"""grep -v "^#" {gtex_vcf} | grep "^{chrom}\\s" | awk '(length($4)>1 || length($5)>1) && ${sample_idx} !~ /0\\/0/ && ${sample_idx} !~ /\\.\\/\\./ {{ count++ }} END {{ print "Indels\\t" count }}'"""
        
        # Combined command
        cmd = f"echo -e 'Type\\tCount' > {chrom_file} && {snp_cmd} >> {chrom_file} && {indel_cmd} >> {chrom_file}"
        
        try:
            process = subprocess.Popen(cmd, shell=True)
            
            # Monitor with timeout
            start_time = time.time()
            while process.poll() is None:
                if time.time() - start_time > max_time:
                    process.kill()
                    error_msg = f"Timeout processing {sample_id} {chrom} after {max_time} seconds"
                    print(error_msg)
                    logging.error(error_msg)
                    break
                time.sleep(5)
            
            if process.returncode == 0:
                # Read results
                with open(chrom_file, 'r') as f:
                    lines = f.readlines()
                    for line in lines[1:]:  # Skip header
                        parts = line.strip().split('\t')
                        if parts[0] == "SNPs":
                            snp_count = int(parts[1])
                            total_snps += snp_count
                            print(f"{sample_id} {chrom}: {snp_count} SNPs")
                            logging.info(f"{sample_id} {chrom}: {snp_count} SNPs")
                        elif parts[0] == "Indels":
                            indel_count = int(parts[1])
                            total_indels += indel_count
                            print(f"{sample_id} {chrom}: {indel_count} Indels")
                            logging.info(f"{sample_id} {chrom}: {indel_count} Indels")
            else:
                error_msg = f"Error processing {sample_id} {chrom}, return code: {process.returncode}"
                print(error_msg)
                logging.error(error_msg)
        
        except Exception as e:
            error_msg = f"Error processing {sample_id} {chrom}: {e}"
            print(error_msg)
            logging.error(error_msg)
    
    # Write final results
    with open(output_file, 'w') as f:
        f.write("Type\tCount\n")
        f.write(f"SNPs\t{total_snps}\n")
        f.write(f"Indels\t{total_indels}\n")
    
    print(f"Final counts for {sample_id}: {total_snps} SNPs, {total_indels} Indels")
    logging.info(f"Final counts for {sample_id}: {total_snps} SNPs, {total_indels} Indels")
    
    return {
        'sample_id': sample_id,
        'num_snps': total_snps,
        'num_indels': total_indels,
        'total_variants': total_snps + total_indels
    }

def main():
    parser = argparse.ArgumentParser(description="Count variants in GTEx samples using direct grep/awk")
    parser.add_argument("--samples", type=int, default=5, help="Number of samples to process (default: 5)")
    parser.add_argument("--random", action="store_true", help="Select random samples instead of first N")
    parser.add_argument("--by-chromosome", action="store_true", help="Process by chromosome (recommended)")
    parser.add_argument("--max-time", type=int, default=1800, help="Maximum time per chromosome in seconds (default: 1800)")
    args = parser.parse_args()
    
    # Get list of samples
    all_samples = get_gtex_samples()
    
    if not all_samples:
        print("No samples found, exiting")
        logging.error("No samples found, exiting")
        return
    
    # Select samples to process
    if args.random:
        if args.samples < len(all_samples):
            samples_to_process = random.sample(all_samples, args.samples)
            print(f"Processing {args.samples} random samples out of {len(all_samples)}")
            logging.info(f"Processing {args.samples} random samples out of {len(all_samples)}")
        else:
            samples_to_process = all_samples
            print(f"Processing all {len(all_samples)} samples")
            logging.info(f"Processing all {len(all_samples)} samples")
    else:
        if args.samples < len(all_samples):
            samples_to_process = all_samples[:args.samples]
            print(f"Processing first {args.samples} samples out of {len(all_samples)}")
            logging.info(f"Processing first {args.samples} samples out of {len(all_samples)}")
        else:
            samples_to_process = all_samples
            print(f"Processing all {len(all_samples)} samples")
            logging.info(f"Processing all {len(all_samples)} samples")
    
    # Process samples
    results = []
    
    if args.by_chromosome:
        for sample_id in tqdm(samples_to_process, desc="Processing samples by chromosome"):
            result = process_by_chromosome(sample_id, max_time=args.max_time)
            results.append(result)
    else:
        for sample_id in tqdm(samples_to_process, desc="Processing samples"):
            success, _ = count_variants_direct(sample_id, max_time=args.max_time*24)  # longer timeout for full file
            if success:
                # Parse results
                output_file = os.path.join(output_dir, f"{sample_id}.counts.txt")
                if os.path.exists(output_file):
                    with open(output_file, 'r') as f:
                        lines = f.readlines()
                        snp_count = 0
                        indel_count = 0
                        for line in lines[1:]:  # Skip header
                            parts = line.strip().split('\t')
                            if parts[0] == "SNPs":
                                snp_count = int(parts[1])
                            elif parts[0] == "Indels":
                                indel_count = int(parts[1])
                    
                    results.append({
                        'sample_id': sample_id,
                        'num_snps': snp_count,
                        'num_indels': indel_count,
                        'total_variants': snp_count + indel_count
                    })
    
    # Create summary DataFrame
    if results:
        df = pd.DataFrame(results)
        df.to_csv(summary_file, sep='\t', index=False)
        
        print(f"\nProcessed {len(df)} GTEx samples")
        print(f"Average SNPs per sample: {df['num_snps'].mean():,.2f}")
        print(f"Average indels per sample: {df['num_indels'].mean():,.2f}")
        print(f"Average total variants per sample: {df['total_variants'].mean():,.2f}")
        
        logging.info(f"Processed {len(df)} GTEx samples")
        logging.info(f"Average SNPs per sample: {df['num_snps'].mean():,.2f}")
        logging.info(f"Average indels per sample: {df['num_indels'].mean():,.2f}")
        logging.info(f"Average total variants per sample: {df['total_variants'].mean():,.2f}")
        
        print(f"Summary saved to {summary_file}")
        logging.info(f"Summary saved to {summary_file}")
    else:
        print("No valid results generated")
        logging.error("No valid results generated")

if __name__ == "__main__":
    main()