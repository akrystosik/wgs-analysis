#!/usr/bin/env python3

import subprocess
import logging
import os
from pathlib import Path

# Setup logging
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def subsample_fastq_pair(r1_path, r2_path, output_dir, num_reads=100000):
    """
    Subsample a pair of FASTQ files using seqtk.
    Uses the same seed for both files to maintain pair consistency.
    """
    try:
        # Check for seqtk
        if subprocess.run(['which', 'seqtk'], capture_output=True).returncode != 0:
            logger.info("Installing seqtk...")
            subprocess.run(['apt-get', 'install', '-y', 'seqtk'], check=True)
        
        # Create output directory
        os.makedirs(output_dir, exist_ok=True)
        
        # Generate output paths
        r1_output = os.path.join(output_dir, f"test_subset_{os.path.basename(r1_path)}")
        r2_output = os.path.join(output_dir, f"test_subset_{os.path.basename(r2_path)}")
        
        # Subsample R1
        logger.info(f"Subsampling R1 to {num_reads} reads...")
        cmd_r1 = f"seqtk sample -s100 {r1_path} {num_reads} | gzip > {r1_output}"
        subprocess.run(cmd_r1, shell=True, check=True)
        
        # Subsample R2 with same seed
        logger.info(f"Subsampling R2 to {num_reads} reads...")
        cmd_r2 = f"seqtk sample -s100 {r2_path} {num_reads} | gzip > {r2_output}"
        subprocess.run(cmd_r2, shell=True, check=True)
        
        # Verify outputs
        for output_file in [r1_output, r2_output]:
            if not os.path.exists(output_file) or os.path.getsize(output_file) == 0:
                raise Exception(f"Failed to create subsampled FASTQ: {output_file}")
            logger.info(f"Successfully created: {output_file}")
            # Log file size
            size_mb = os.path.getsize(output_file) / (1024 * 1024)
            logger.info(f"File size: {size_mb:.2f} MB")
        
        return r1_output, r2_output
        
    except Exception as e:
        logger.error(f"Subsampling failed: {str(e)}")
        raise

def main():
    # Input files
    r1_path = "encode_analysis/data/fastq/K562/K562_WGS_pair1_R1.fastq.gz"
    r2_path = "encode_analysis/data/fastq/K562/K562_WGS_pair1_R2.fastq.gz"
    
    # Output directory for test data
    output_dir = "encode_analysis/data/test_data/K562"
    
    # Create subsampled files
    subsampled_r1, subsampled_r2 = subsample_fastq_pair(r1_path, r2_path, output_dir)
    
    logger.info("Subsampling complete. You can now run the pipeline test with:")
    logger.info(f"python pipeline_test.py --pair-id K562_test --input-fastq {subsampled_r1} {subsampled_r2} --force --force-qc")

if __name__ == "__main__":
    main()