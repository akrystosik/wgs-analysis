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
log_file = os.path.join(log_dir, f"bcftools_direct_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")
logging.basicConfig(
    filename=log_file,
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# Define paths
#gtex_vcf = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex-3/project_gene_regulation/data/GTEx/WGS/phg001796.v1.GTEx_v9_WGS_953.genotype-calls-vcf.c1/GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv.vcf"
gtex_vcf = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex-3/project_gene_regulation/data/GTEx/WGS/phg001796.v1.GTEx_v9_WGS_phased.genotype-calls-vcf.c1/GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv.SHAPEIT2_phased.vcf.gz"

output_dir = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/wgs/data/qc/variants/bcftools_direct"
summary_file = os.path.join(output_dir, "gtex_sample_stats.tsv")

# Create output directory
os.makedirs(output_dir, exist_ok=True)

def get_gtex_samples():
    """Get list of all samples in the GTEx VCF file"""
    print("Fetching sample IDs from GTEx VCF...")
    logging.info("Fetching sample IDs from GTEx VCF...")
    
    cmd = f"bcftools query -l {gtex_vcf}"
    try:
        result = subprocess.check_output(cmd, shell=True, text=True)
        samples = result.strip().split('\n')
        print(f"Found {len(samples)} GTEx samples")
        logging.info(f"Found {len(samples)} GTEx samples")
        return samples
    except subprocess.CalledProcessError as e:
        print(f"Error fetching GTEx samples: {e}")
        logging.error(f"Error fetching GTEx samples: {e}")
        return []

def count_variants_bcftools(sample_id, max_time=3600):
    """Count variants for a sample using bcftools"""
    output_file = os.path.join(output_dir, f"{sample_id}.stats.txt")
    
    # Check if output file already exists
    if os.path.exists(output_file):
        print(f"Stats for {sample_id} already exist, skipping")
        logging.info(f"Stats for {sample_id} already exist, skipping")
        return True, "Already processed"
    
    print(f"Counting variants for {sample_id} using bcftools...")
    logging.info(f"Counting variants for {sample_id} using bcftools...")
    
    # Use bcftools to directly count variants with GT filter
    cmd = f"bcftools view -s {sample_id} -i 'GT=\"alt\"' {gtex_vcf} | bcftools stats > {output_file}"
    
    try:
        # Start the process
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
            time.sleep(5)  # Check every 5 seconds
        
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

def parse_stats_file(stats_file):
    """Parse bcftools stats output file for a sample"""
    if not os.path.exists(stats_file):
        return None
    
    sample_id = os.path.basename(stats_file).replace('.stats.txt', '')
    
    results = {
        'sample_id': sample_id,
        'num_snps': 0,
        'num_indels': 0,
        'ts_tv_ratio': 0,
        'total_variants': 0
    }
    
    try:
        # Parse the stats file
        with open(stats_file, 'r') as f:
            content = f.read()
            
            # Extract SNP count
            snp_match = re.search(r'number of SNPs:\s+(\d+)', content)
            if snp_match:
                results['num_snps'] = int(snp_match.group(1))
            
            # Extract indel count
            indel_match = re.search(r'number of indels:\s+(\d+)', content)
            if indel_match:
                results['num_indels'] = int(indel_match.group(1))
            
            # Extract Ts/Tv ratio
            tstv_match = re.search(r'ts/tv:\s+([\d\.]+)', content)
            if tstv_match:
                results['ts_tv_ratio'] = float(tstv_match.group(1))
            
            # Calculate total variants
            results['total_variants'] = results['num_snps'] + results['num_indels']
        
        return results
    
    except Exception as e:
        print(f"Error parsing stats file {stats_file}: {e}")
        logging.error(f"Error parsing stats file {stats_file}: {e}")
        return None

def process_sample_wrapper(sample_id):
    """Wrapper function for multiprocessing"""
    return sample_id, count_variants_bcftools(sample_id)

def generate_summary():
    """Generate summary statistics from all processed samples"""
    print("Generating summary statistics...")
    logging.info("Generating summary statistics...")
    
    # Find all stats files
    stats_files = [f for f in os.listdir(output_dir) if f.endswith('.stats.txt')]
    
    if not stats_files:
        print("No stats files found")
        logging.error("No stats files found")
        return
    
    results = []
    
    for stats_file in tqdm(stats_files, desc="Parsing stats files"):
        stats = parse_stats_file(os.path.join(output_dir, stats_file))
        if stats:
            results.append(stats)
    
    # Create summary DataFrame
    df = pd.DataFrame(results)
    
    # Save summary to file
    df.to_csv(summary_file, sep='\t', index=False)
    
    # Print basic statistics
    print(f"\nProcessed {len(df)} GTEx samples")
    print(f"Average SNPs per sample: {df['num_snps'].mean():,.2f}")
    print(f"Average indels per sample: {df['num_indels'].mean():,.2f}")
    print(f"Average total variants per sample: {df['total_variants'].mean():,.2f}")
    print(f"Average Ts/Tv ratio: {df['ts_tv_ratio'].mean():.4f}")
    
    logging.info(f"Processed {len(df)} GTEx samples")
    logging.info(f"Average SNPs per sample: {df['num_snps'].mean():,.2f}")
    logging.info(f"Average indels per sample: {df['num_indels'].mean():,.2f}")
    logging.info(f"Average total variants per sample: {df['total_variants'].mean():,.2f}")
    logging.info(f"Average Ts/Tv ratio: {df['ts_tv_ratio'].mean():.4f}")
    
    print(f"Summary saved to {summary_file}")
    logging.info(f"Summary saved to {summary_file}")

def process_by_chromosome(sample_id, chroms=None):
    """Process a sample chromosome by chromosome"""
    if chroms is None:
        chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    
    output_dir_by_chrom = os.path.join(output_dir, "by_chromosome", sample_id)
    os.makedirs(output_dir_by_chrom, exist_ok=True)
    
    results = {
        'sample_id': sample_id,
        'num_snps': 0,
        'num_indels': 0,
        'total_variants': 0
    }
    
    for chrom in tqdm(chroms, desc=f"Processing {sample_id} by chromosome"):
        output_file = os.path.join(output_dir_by_chrom, f"{chrom}.stats.txt")
        
        # Skip if already processed
        if os.path.exists(output_file):
            print(f"Stats for {sample_id} {chrom} already exist, skipping")
            continue
        
        # Process this chromosome
        cmd = f"bcftools view -s {sample_id} -i 'GT=\"alt\"' -r {chrom} {gtex_vcf} | bcftools stats > {output_file}"
        
        try:
            subprocess.run(cmd, shell=True, check=True)
            
            # Parse chromosome stats
            with open(output_file, 'r') as f:
                content = f.read()
                
                # Extract SNP count
                snp_match = re.search(r'number of SNPs:\s+(\d+)', content)
                if snp_match:
                    snp_count = int(snp_match.group(1))
                    results['num_snps'] += snp_count
                
                # Extract indel count
                indel_match = re.search(r'number of indels:\s+(\d+)', content)
                if indel_match:
                    indel_count = int(indel_match.group(1))
                    results['num_indels'] += indel_count
        
        except Exception as e:
            print(f"Error processing {sample_id} {chrom}: {e}")
            logging.error(f"Error processing {sample_id} {chrom}: {e}")
    
    # Calculate total variants
    results['total_variants'] = results['num_snps'] + results['num_indels']
    
    # Save combined results
    df = pd.DataFrame([results])
    df.to_csv(os.path.join(output_dir, f"{sample_id}.combined.tsv"), sep='\t', index=False)
    
    return results

def main():
    parser = argparse.ArgumentParser(description="Count variants in GTEx samples using bcftools directly")
    parser.add_argument("--samples", type=int, default=10, help="Number of samples to process (default: 10)")
    parser.add_argument("--random", action="store_true", help="Select random samples instead of first N")
    parser.add_argument("--threads", type=int, default=4, help="Number of parallel threads (default: 4)")
    parser.add_argument("--max-time", type=int, default=3600, help="Maximum time per sample in seconds (default: 3600)")
    parser.add_argument("--summary-only", action="store_true", help="Only generate summary from existing results")
    parser.add_argument("--by-chromosome", action="store_true", help="Process samples chromosome by chromosome")
    args = parser.parse_args()
    
    if args.summary_only:
        generate_summary()
        return
    
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
    if args.by_chromosome:
        # Process chromosome by chromosome
        results = []
        for sample_id in tqdm(samples_to_process, desc="Processing samples by chromosome"):
            result = process_by_chromosome(sample_id)
            results.append(result)
        
        # Create summary DataFrame
        df = pd.DataFrame(results)
        df.to_csv(summary_file, sep='\t', index=False)
        
        print(f"\nProcessed {len(df)} GTEx samples by chromosome")
        print(f"Average SNPs per sample: {df['num_snps'].mean():,.2f}")
        print(f"Average indels per sample: {df['num_indels'].mean():,.2f}")
        print(f"Average total variants per sample: {df['total_variants'].mean():,.2f}")
    else:
        # Process samples in parallel
        max_workers = min(args.threads, multiprocessing.cpu_count())
        results = {}
        
        if max_workers > 1:
            print(f"Processing samples using {max_workers} workers")
            logging.info(f"Processing samples using {max_workers} workers")
            
            with multiprocessing.Pool(max_workers) as pool:
                for sample_id, result in pool.imap_unordered(process_sample_wrapper, samples_to_process):
                    results[sample_id] = result
        else:
            print("Processing samples sequentially")
            logging.info("Processing samples sequentially")
            
            for sample_id in tqdm(samples_to_process, desc="Processing samples"):
                results[sample_id] = count_variants_bcftools(sample_id, args.max_time)
        
        # Count successes and failures
        successes = sum(1 for result, _ in results.values() if result)
        failures = sum(1 for result, _ in results.values() if not result)
        
        print(f"Processing complete: {successes} successful, {failures} failed")
        logging.info(f"Processing complete: {successes} successful, {failures} failed")
        
        # Generate summary
        generate_summary()

if __name__ == "__main__":
    main()