#!/usr/bin/env python
import os
import subprocess
import glob
import pandas as pd
from tqdm import tqdm
import time
import logging
from datetime import datetime
import argparse
import multiprocessing

# Set up logging
log_dir = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/wgs/data/qc/variants/logs"
os.makedirs(log_dir, exist_ok=True)
log_file = os.path.join(log_dir, f"gtex_extraction_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")
logging.basicConfig(
    filename=log_file,
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# Define paths
gtex_vcf = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex-3/project_gene_regulation/data/GTEx/WGS/phg001796.v1.GTEx_v9_WGS_953.genotype-calls-vcf.c1/GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv.vcf"
output_dir = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/wgs/data/variants/gtex_samples"
stats_dir = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/wgs/data/qc/variants/gtex_samples"
summary_file = os.path.join(stats_dir, "gtex_sample_stats.tsv")

# Create output directories
os.makedirs(output_dir, exist_ok=True)
os.makedirs(stats_dir, exist_ok=True)

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

def extract_sample(sample_id, max_time=3600):
    """Extract a single sample from the GTEx VCF file"""
    output_file = os.path.join(output_dir, f"{sample_id}.vcf.gz")
    stats_file = os.path.join(stats_dir, f"{sample_id}.stats.txt")
    
    # Skip if both files already exist
    if os.path.exists(output_file) and os.path.exists(stats_file):
        print(f"Sample {sample_id} already processed, skipping")
        logging.info(f"Sample {sample_id} already processed, skipping")
        return True, "Already processed"
    
    print(f"Extracting sample {sample_id}...")
    logging.info(f"Extracting sample {sample_id}...")
    
    # Extract the sample
    extract_cmd = f"bcftools view -s {sample_id} -Oz -o {output_file} {gtex_vcf}"
    
    try:
        # Start the extraction process
        process = subprocess.Popen(extract_cmd, shell=True)
        
        # Monitor the process with timeout
        start_time = time.time()
        while process.poll() is None:
            if time.time() - start_time > max_time:
                process.kill()
                error_msg = f"Timeout extracting sample {sample_id} after {max_time} seconds"
                print(error_msg)
                logging.error(error_msg)
                return False, error_msg
            time.sleep(1)
        
        # Check extraction result
        if process.returncode != 0:
            error_msg = f"Error extracting sample {sample_id}, return code: {process.returncode}"
            print(error_msg)
            logging.error(error_msg)
            return False, error_msg
        
        # Index the extracted VCF
        index_cmd = f"bcftools index {output_file}"
        subprocess.run(index_cmd, shell=True, check=True)
        
        # Run bcftools stats on the extracted VCF
        stats_cmd = f"bcftools stats -s {sample_id} {output_file} > {stats_file}"
        subprocess.run(stats_cmd, shell=True, check=True)
        
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
        'ts_tv_ratio': 0
    }
    
    try:
        # Extract stats using grep
        grep_cmd = f"grep '^SN' {stats_file}"
        grep_output = subprocess.check_output(grep_cmd, shell=True, text=True)
        
        for line in grep_output.splitlines():
            if 'number of SNPs:' in line:
                results['num_snps'] = int(line.split('\t')[-1])
            elif 'number of indels:' in line:
                results['num_indels'] = int(line.split('\t')[-1])
            elif 'ts/tv:' in line:
                try:
                    results['ts_tv_ratio'] = float(line.split('\t')[-1])
                except ValueError:
                    results['ts_tv_ratio'] = 0
        
        return results
    
    except Exception as e:
        print(f"Error parsing stats file {stats_file}: {e}")
        logging.error(f"Error parsing stats file {stats_file}: {e}")
        return None

def process_sample_wrapper(sample_id):
    """Wrapper function for multiprocessing"""
    return sample_id, extract_sample(sample_id)

def generate_summary():
    """Generate summary statistics from all processed samples"""
    print("Generating summary statistics...")
    logging.info("Generating summary statistics...")
    
    # Find all stats files
    stats_files = glob.glob(os.path.join(stats_dir, "*.stats.txt"))
    
    if not stats_files:
        print("No stats files found")
        logging.error("No stats files found")
        return
    
    results = []
    
    for stats_file in tqdm(stats_files, desc="Parsing stats files"):
        stats = parse_stats_file(stats_file)
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
    print(f"Average Ts/Tv ratio: {df['ts_tv_ratio'].mean():.4f}")
    
    logging.info(f"Processed {len(df)} GTEx samples")
    logging.info(f"Average SNPs per sample: {df['num_snps'].mean():,.2f}")
    logging.info(f"Average indels per sample: {df['num_indels'].mean():,.2f}")
    logging.info(f"Average Ts/Tv ratio: {df['ts_tv_ratio'].mean():.4f}")
    
    print(f"Summary saved to {summary_file}")
    logging.info(f"Summary saved to {summary_file}")

def main():
    parser = argparse.ArgumentParser(description="Extract and analyze individual GTEx samples")
    parser.add_argument("--samples", type=int, default=100, help="Number of samples to process (default: 100)")
    parser.add_argument("--threads", type=int, default=4, help="Number of parallel threads (default: 4)")
    parser.add_argument("--max-time", type=int, default=3600, help="Maximum time per sample in seconds (default: 3600)")
    parser.add_argument("--summary-only", action="store_true", help="Only generate summary from existing stats files")
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
    
    # Limit to specified number if needed
    if args.samples < len(all_samples):
        samples_to_process = all_samples[:args.samples]
        print(f"Processing first {args.samples} samples out of {len(all_samples)}")
        logging.info(f"Processing first {args.samples} samples out of {len(all_samples)}")
    else:
        samples_to_process = all_samples
        print(f"Processing all {len(all_samples)} samples")
        logging.info(f"Processing all {len(all_samples)} samples")
    
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
            results[sample_id] = extract_sample(sample_id, args.max_time)
    
    # Count successes and failures
    successes = sum(1 for result, _ in results.values() if result)
    failures = sum(1 for result, _ in results.values() if not result)
    
    print(f"Processing complete: {successes} successful, {failures} failed")
    logging.info(f"Processing complete: {successes} successful, {failures} failed")
    
    # Generate summary
    generate_summary()

if __name__ == "__main__":
    main()