import os
import subprocess
import glob
import re
import pandas as pd
from concurrent.futures import ProcessPoolExecutor
from tqdm import tqdm
import time
import logging
from datetime import datetime

# Set up logging
log_dir = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/wgs/data/qc/variants/logs"
os.makedirs(log_dir, exist_ok=True)
log_file = os.path.join(log_dir, f"bcftools_stats_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log")
logging.basicConfig(
    filename=log_file,
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

# Define input and output directories
cell_line_dir = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/wgs/data/variants/cell_line_analysis_v5"
gtex_dir = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex-3/project_gene_regulation/data/GTEx/WGS/phg001796.v1.GTEx_v9_WGS_953.genotype-calls-vcf.c1"
output_dir = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/wgs/data/qc/variants"
summary_file = os.path.join(output_dir, "bcftools_summary_report.tsv")

# Settings
max_workers = 32  # Use up to 32 cores
timeout_seconds = 1800  # 30 minute timeout per file
skip_existing = True    # Skip files that already have stats output
force_regenerate = False  # Set to True to regenerate any empty or problematic stats files

# Create output directory if it doesn't exist
os.makedirs(output_dir, exist_ok=True)

# Function to find all VCF files in the specific directory structure
def find_vcf_files():
    all_vcfs = []
    
    # Check cell line directory
    if os.path.exists(cell_line_dir):
        print(f"Searching for cell line VCF files in: {cell_line_dir}")
        logging.info(f"Searching for cell line VCF files in: {cell_line_dir}")
        
        # List all subdirectories (cell lines)
        cell_lines = [d for d in os.listdir(cell_line_dir) if os.path.isdir(os.path.join(cell_line_dir, d))]
        
        for cell_line in cell_lines:
            final_dir = os.path.join(cell_line_dir, cell_line, "final")
            if os.path.exists(final_dir):
                # Only look for combined_variants.vcf.gz files, not .tbi files
                vcf_pattern = os.path.join(final_dir, "combined_variants.vcf.gz")
                if os.path.exists(vcf_pattern):
                    all_vcfs.append(vcf_pattern)
                    print(f"  Found VCF file for {cell_line}: {vcf_pattern}")
                    logging.info(f"  Found VCF file for {cell_line}: {vcf_pattern}")
                else:
                    print(f"  No VCF file found for {cell_line}")
                    logging.info(f"  No VCF file found for {cell_line}")
    else:
        print(f"Warning: Cell line directory does not exist: {cell_line_dir}")
        logging.warning(f"Cell line directory does not exist: {cell_line_dir}")
    
    # Check GTEx directory
    if os.path.exists(gtex_dir):
        print(f"Searching for GTEx VCF files in: {gtex_dir}")
        logging.info(f"Searching for GTEx VCF files in: {gtex_dir}")
        
        # If gtex_dir is a file itself
        if os.path.isfile(gtex_dir) and (gtex_dir.endswith('.vcf') or gtex_dir.endswith('.vcf.gz')):
            all_vcfs.append(gtex_dir)
            print(f"  Found GTEx VCF file: {gtex_dir}")
            logging.info(f"  Found GTEx VCF file: {gtex_dir}")
        else:
            # Search for VCF files in GTEx directory (only .vcf and .vcf.gz, not .tbi)
            gtex_vcfs = []
            for file in os.listdir(gtex_dir):
                if file.endswith('.vcf') or file.endswith('.vcf.gz'):
                    if not file.endswith('.tbi'):  # Skip index files
                        gtex_vcfs.append(os.path.join(gtex_dir, file))
            
            if gtex_vcfs:
                all_vcfs.extend(gtex_vcfs)
                print(f"  Found {len(gtex_vcfs)} GTEx VCF file(s)")
                logging.info(f"  Found {len(gtex_vcfs)} GTEx VCF file(s)")
                for vcf in gtex_vcfs:
                    print(f"    {vcf}")
                    logging.info(f"    {vcf}")
            else:
                print("  No VCF files found directly in GTEx directory, attempting recursive search")
                logging.info("  No VCF files found directly in GTEx directory, attempting recursive search")
                
                # Try a recursive search, excluding .tbi files
                for root, dirs, files in os.walk(gtex_dir):
                    for file in files:
                        if (file.endswith('.vcf') or file.endswith('.vcf.gz')) and not file.endswith('.tbi'):
                            all_vcfs.append(os.path.join(root, file))
                            print(f"    Found {os.path.join(root, file)}")
                            logging.info(f"    Found {os.path.join(root, file)}")
    else:
        print(f"Warning: GTEx directory does not exist: {gtex_dir}")
        logging.warning(f"GTEx directory does not exist: {gtex_dir}")
    
    print(f"Found a total of {len(all_vcfs)} VCF files")
    logging.info(f"Found a total of {len(all_vcfs)} VCF files")
    
    return all_vcfs

# Function to check if a stats file is valid
def is_valid_stats_file(stats_file):
    if not os.path.exists(stats_file):
        return False
        
    # Check if file is empty
    if os.path.getsize(stats_file) == 0:
        return False
        
    # Check if file contains SN lines (summary numbers)
    try:
        grep_cmd = f"grep -c '^SN' {stats_file}"
        count = int(subprocess.check_output(grep_cmd, shell=True, text=True).strip())
        return count > 0
    except:
        return False

# Function to run bcftools stats on a VCF file
def run_bcftools_stats(vcf_file):
    filename = os.path.basename(vcf_file)
    cell_line = os.path.basename(os.path.dirname(os.path.dirname(vcf_file))) if "cell_line_analysis_v5" in vcf_file else "GTEx"
    
    # Create a more informative output filename that includes cell line info
    if "cell_line_analysis_v5" in vcf_file:
        output_file = os.path.join(output_dir, f"{cell_line}_{filename}.stats.txt")
    else:
        output_file = os.path.join(output_dir, f"{filename}.stats.txt")
    
    # Skip if file already exists and is valid and skip_existing is True
    if skip_existing and os.path.exists(output_file) and is_valid_stats_file(output_file) and not force_regenerate:
        print(f"Skipping {filename} for {cell_line} - valid stats file already exists")
        logging.info(f"Skipping {filename} for {cell_line} - valid stats file already exists")
        return (vcf_file, output_file, True)
    
    # File exists but is invalid or force_regenerate is True
    if os.path.exists(output_file) and (not is_valid_stats_file(output_file) or force_regenerate):
        print(f"Regenerating stats for {filename} for {cell_line} - previous file was invalid or regeneration forced")
        logging.info(f"Regenerating stats for {filename} for {cell_line} - previous file was invalid or regeneration forced")
    
    # Fix the bcftools command - DO NOT use -1 parameter as it's causing issues
    # Instead, use standard bcftools stats command
    cmd = f"bcftools stats {vcf_file} > {output_file}"
    
    try:
        logging.info(f"Starting bcftools stats on {filename} for {cell_line}")
        print(f"Processing {filename} for {cell_line}...")
        
        # Run with timeout to prevent hanging
        process = subprocess.Popen(cmd, shell=True, executable='/bin/bash')
        
        # Track start time
        start_time = time.time()
        
        # Wait for process with timeout
        while process.poll() is None:
            if time.time() - start_time > timeout_seconds:
                process.kill()
                error_msg = f"Timeout processing {filename} for {cell_line} after {timeout_seconds} seconds"
                print(error_msg)
                logging.error(error_msg)
                return (vcf_file, output_file, False)
            time.sleep(1)
        
        if process.returncode == 0:
            # Verify the output file is valid
            if is_valid_stats_file(output_file):
                logging.info(f"Successfully processed {filename} for {cell_line}")
                return (vcf_file, output_file, True)
            else:
                error_msg = f"Generated stats file for {filename} for {cell_line} is invalid"
                print(error_msg)
                logging.error(error_msg)
                return (vcf_file, output_file, False)
        else:
            error_msg = f"Error processing {filename} for {cell_line}, return code: {process.returncode}"
            print(error_msg)
            logging.error(error_msg)
            return (vcf_file, output_file, False)
            
    except Exception as e:
        error_msg = f"Error processing {vcf_file}: {e}"
        print(error_msg)
        logging.error(error_msg)
        return (vcf_file, output_file, False)

# Function to parse bcftools stats output file efficiently
def parse_stats_file(stats_file):
    if not os.path.exists(stats_file) or not is_valid_stats_file(stats_file):
        logging.warning(f"Skipping invalid stats file: {stats_file}")
        return None
    
    # Extract cell line from filename
    filename = os.path.basename(stats_file).replace('.stats.txt', '')
    cell_line = filename.split('_')[0] if '_' in filename else "Unknown"
    
    # Try to find the original VCF file
    vcf_filename = filename
    if '_combined_variants.vcf' in filename:
        vcf_filename = filename.split('_', 1)[1]  # Remove cell line prefix
    
    # Try different possible locations for the original VCF
    possible_paths = [
        stats_file.replace('.stats.txt', ''),
        os.path.join(cell_line_dir, cell_line, "final", vcf_filename),
        os.path.join(gtex_dir, vcf_filename)
    ]
    
    vcf_file = None
    for path in possible_paths:
        if os.path.exists(path):
            vcf_file = path
            break
    
    file_size_mb = 0
    if vcf_file and os.path.exists(vcf_file):
        file_size_mb = os.path.getsize(vcf_file) / (1024 * 1024)
    
    results = {
        'filename': filename,
        'cell_line': cell_line,
        'num_samples': 0,
        'num_records': 0,
        'num_SNPs': 0,
        'num_indels': 0,
        'num_multiallelic': 0,
        'ts_tv_ratio': 0,
        'file_size_MB': file_size_mb
    }
    
    try:
        # Use grep to extract only the summary numbers
        grep_cmd = f"grep '^SN' {stats_file}"
        grep_output = subprocess.check_output(grep_cmd, shell=True, text=True)
        
        # Extract values
        for line in grep_output.splitlines():
            if 'number of samples:' in line:
                results['num_samples'] = int(line.split('\t')[-1])
            elif 'number of records:' in line:
                results['num_records'] = int(line.split('\t')[-1])
            elif 'number of SNPs:' in line:
                results['num_SNPs'] = int(line.split('\t')[-1])
            elif 'number of indels:' in line:
                results['num_indels'] = int(line.split('\t')[-1])
            elif 'number of multiallelic sites:' in line:
                results['num_multiallelic'] = int(line.split('\t')[-1])
            elif 'ts/tv:' in line:
                results['ts_tv_ratio'] = float(line.split('\t')[-1])
                
        return results
    except Exception as e:
        error_msg = f"Error parsing stats file {stats_file}: {e}"
        print(error_msg)
        logging.error(error_msg)
        return None

# Find all VCF files based on the specific directory structure
all_vcfs = find_vcf_files()

if not all_vcfs:
    print("No VCF files found. Please check the directory paths.")
    logging.error("No VCF files found. Please check the directory paths.")
    exit(1)

print(f"Starting to process {len(all_vcfs)} VCF files...")
logging.info(f"Starting to process {len(all_vcfs)} VCF files...")

# Process VCF files in parallel with progress bar
results = []

# Use parallel processing with the optimized number of workers
with ProcessPoolExecutor(max_workers=max_workers) as executor:
    futures = [executor.submit(run_bcftools_stats, vcf) for vcf in all_vcfs]
    
    # Track progress with tqdm
    for future in tqdm(futures, total=len(futures), desc="Processing VCF files"):
        results.append(future.result())

# Count successful and failed processes
successful = sum(1 for _, _, success in results if success)
failed = sum(1 for _, _, success in results if not success)

print(f"Processing complete: {successful} successful, {failed} failed.")
logging.info(f"Processing complete: {successful} successful, {failed} failed.")

# Parse all stats files and create summary report
stats_files = [output_file for _, output_file, success in results if success]
parsed_stats = []

print("Creating summary report...")
logging.info("Creating summary report...")

for stats_file in tqdm(stats_files, desc="Parsing stats files"):
    stats = parse_stats_file(stats_file)
    if stats:
        parsed_stats.append(stats)

# Create summary DataFrame and save to TSV
if parsed_stats:
    df = pd.DataFrame(parsed_stats)
    df.to_csv(summary_file, sep='\t', index=False)
    print(f"Summary report created at: {summary_file}")
    logging.info(f"Summary report created at: {summary_file}")
    
    # Print basic summary statistics
    print("\nSummary Statistics:")
    print(f"Total VCF files processed: {len(parsed_stats)}")
    print(f"Total cell lines: {df['cell_line'].nunique()}")
    print(f"Total samples across all files: {df['num_samples'].sum():,}")
    print(f"Total SNPs across all files: {df['num_SNPs'].sum():,}")
    print(f"Total Indels across all files: {df['num_indels'].sum():,}")
    
    # Average stats
    print(f"Average samples per file: {df['num_samples'].mean():.2f}")
    print(f"Average SNPs per file: {df['num_SNPs'].mean():,.0f}")
    print(f"Average Indels per file: {df['num_indels'].mean():,.0f}")
    print(f"Average Ts/Tv ratio: {df['ts_tv_ratio'].mean():.4f}")
    
    # Stats by cell line
    print("\nStats by Cell Line:")
    cell_line_stats = df.groupby('cell_line').agg({
        'num_samples': 'mean',
        'num_SNPs': 'sum',
        'num_indels': 'sum',
        'ts_tv_ratio': 'mean'
    }).reset_index()
    
    for _, row in cell_line_stats.iterrows():
        print(f"  {row['cell_line']}: {row['num_SNPs']:,} SNPs, {row['num_indels']:,} Indels, Ts/Tv: {row['ts_tv_ratio']:.4f}")
    
    # Log the same information
    logging.info("\nSummary Statistics:")
    logging.info(f"Total VCF files processed: {len(parsed_stats)}")
    logging.info(f"Total cell lines: {df['cell_line'].nunique()}")
    logging.info(f"Total samples across all files: {df['num_samples'].sum():,}")
    logging.info(f"Total SNPs across all files: {df['num_SNPs'].sum():,}")
    logging.info(f"Total Indels across all files: {df['num_indels'].sum():,}")
    logging.info(f"Average samples per file: {df['num_samples'].mean():.2f}")
    logging.info(f"Average SNPs per file: {df['num_SNPs'].mean():,.0f}")
    logging.info(f"Average Indels per file: {df['num_indels'].mean():,.0f}")
    logging.info(f"Average Ts/Tv ratio: {df['ts_tv_ratio'].mean():.4f}")
    
    # Print top files by size
    print("\nLargest VCF Files:")
    largest_files = df.sort_values('file_size_MB', ascending=False).head(5)
    logging.info("\nLargest VCF Files:")
    
    for _, row in largest_files.iterrows():
        print(f"{row['cell_line']} - {row['filename']}: {row['file_size_MB']:.2f} MB")
        logging.info(f"{row['cell_line']} - {row['filename']}: {row['file_size_MB']:.2f} MB")
        
else:
    print("No valid stats files were generated.")
    logging.error("No valid stats files were generated.")

print(f"All bcftools stats processing completed. Results saved to {output_dir}")
logging.info(f"All bcftools stats processing completed. Results saved to {output_dir}")
print(f"Log file created: {log_file}")