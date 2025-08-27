# encode_analysis/pipeline/pipeline.py

import subprocess
import logging
import os
import sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import shutil             
import time               
from datetime import datetime
import traceback 
from typing import Dict, List, Optional, Any
import yaml
import argparse

# Add the pipeline directory to Python path
current_dir = os.path.dirname(os.path.abspath(__file__))
pipeline_dir = current_dir  # pipeline.py is now at the pipeline level
sys.path.insert(0, os.path.dirname(current_dir))  # Add encode_analysis to path

# Import local modules
from core.monitoring import WGSMonitor
from core.exceptions import QCGateException
from qc.bam import BAMQCRunner
from qc.fastq import FastQMultiQC
from qc.vcf import run_vcf_qc
from pipeline_io.fastq import process_fastqs
from pipeline_io.downloader import FastqDownloader
from pipeline_io.verify import verify_fastq_pairs
from utils.naming import validate_pipeline_names, get_standardized_names

class WGSPipeline:
    def __init__(self, args):
        # First validate the pair_id format before any processing
        validate_pipeline_names(args.pair_id)
        
        self.args = args
        self.pair_id = args.pair_id
        # Get standardized names that will be used throughout the pipeline
        self.names = get_standardized_names(self.pair_id)
        self.force = args.force
        self.force_qc = args.force_qc
        self.input_fastq = args.input_fastq
        self.input_bam = args.input_bam
        self.environment = args.environment
        
        # Load config file
        self.config = self._load_config(args.config)
        
        # Setup paths
        self.base_dir = Path(self.config['base_dir'])
        self.reference_dir = Path(self.base_dir / self.config['data_dirs']['reference'])
        self.reference_path = Path(self.reference_dir / self.config['reference']['fasta'])

        # Create sample-specific directories
        self.sample_dir = self.base_dir / 'data/bam' / self.names['sample_dir']
        self.sample_dir.mkdir(parents=True, exist_ok=True)
                
        # Setup logging
        self._setup_logging()
        self.logger = logging.getLogger(__name__)
        
        # Verify reference path
        if not self.reference_path.exists():
            raise FileNotFoundError(f"Reference genome not found: {self.reference_path}")
        
        # Initialize components
        self.library_type = self._determine_library_type()
        self.monitor = WGSMonitor(self.config)
        self._setup_directories()
        self.downloader = FastqDownloader(self.config)
        self.fastqc_runner = FastQMultiQC(self.config)
        self.bam_qc_runner = BAMQCRunner(self.config)

    def _load_config(self, config_path: str) -> Dict:
        """Load pipeline configuration file."""
        if not os.path.exists(config_path):
            raise FileNotFoundError(f"Config file not found: {config_path}")
        
        with open(config_path) as f:
            return yaml.safe_load(f)

    def _setup_logging(self):
        """Configure logging to both console and file with timestamps."""
        log_config = self.config.get('logging', {})
        
        # Create a timestamp-based log filename
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        log_file = f"{self.base_dir}/logs/pipeline_{self.pair_id}_{timestamp}.log"
        
        # Create logs directory if it doesn't exist
        os.makedirs(os.path.dirname(log_file), exist_ok=True)
        
        # Set up logging format
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        
        # Set up file handler with DEBUG level
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(logging.DEBUG)
        file_handler.setFormatter(formatter)
        
        # Set up console handler with configured level
        console_handler = logging.StreamHandler()
        console_handler.setLevel(getattr(logging, log_config.get('level', 'INFO')))
        console_handler.setFormatter(formatter)
        
        # Configure root logger
        root_logger = logging.getLogger()
        root_logger.setLevel(logging.DEBUG)  # Allow DEBUG messages to be handled
        root_logger.addHandler(file_handler)
        root_logger.addHandler(console_handler)
        
        # Log the start of logging
        logging.info(f"Pipeline log file created at: {log_file}")    

    def _setup_directories(self):
        """Create necessary directory structure."""
        for dir_name, dir_path in self.config['data_dirs'].items():
            full_path = self.base_dir / dir_path
            full_path.mkdir(parents=True, exist_ok=True)

    def _determine_library_type(self) -> str:
        """Determine library type from args or metadata."""
        if hasattr(self.args, 'library_type') and self.args.library_type:
            return self.args.library_type
        return "pcr"  # Default to PCR library type

    def _get_qc_thresholds(self) -> Dict:
        """Get QC thresholds for current library type."""
        return self.config['library_types'][self.library_type]

    def _run_variant_calling_script(self, script_path: str, args: str) -> bool:
        """Execute variant calling script with improved process monitoring."""
        try:
            cmd = [script_path] + args.split()
            self.logger.info(f"Executing: {' '.join(cmd)}")
            self.logger.info(f"Working directory: {os.getcwd()}")
            
            # Check script existence and permissions
            if not os.path.exists(script_path):
                self.logger.error(f"Script not found: {script_path}")
                return False
            if not os.access(script_path, os.X_OK):
                self.logger.error(f"Script is not executable: {script_path}")
                return False
                
            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
                bufsize=1  # Line buffered
            )
            
            self.logger.info(f"Process started with PID: {process.pid}")
            
            # Monitor output without blocking
            while True:
                stdout_line = process.stdout.readline()
                stderr_line = process.stderr.readline()
                
                if stdout_line:
                    self.logger.info(f"SUBPROCESS STDOUT: {stdout_line.strip()}")
                if stderr_line:
                    # Only log as warning if it's actually an error
                    if 'error' in stderr_line.lower():
                        self.logger.warning(f"SUBPROCESS STDERR: {stderr_line.strip()}")
                    else:
                        self.logger.info(f"SUBPROCESS STDERR: {stderr_line.strip()}")
                
                # Check if process has finished
                if not stdout_line and not stderr_line and process.poll() is not None:
                    break
                    
                # Add periodic status check with timeout
                if process.poll() is not None:
                    self.logger.info(f"Process exited with code {process.returncode}")
                    if process.returncode != 0:
                        self.logger.error(f"Process exited with non-zero code {process.returncode}")
                        return False
                    break
            
            # Log successful completion
            self.logger.info(f"Process completed successfully")
            return True
                
        except Exception as e:
            self.logger.error(f"Script execution failed: {str(e)}")
            self.logger.error(f"Traceback: {traceback.format_exc()}")
            return False

    def run_variant_calling(self, marked_bam: str) -> str:
        """Run appropriate variant caller based on environment."""
        try:
            scripts_dir = self.base_dir / 'pipeline/scripts'

            
            if self.environment == 'deepvariant':
                # Run DeepVariant with full pair_id
                deepvariant_script = scripts_dir / 'run_deepvariant_analysis.sh'
                self.logger.info(f"Running DeepVariant for {self.pair_id}")
                success = self._run_variant_calling_script(str(deepvariant_script), self.pair_id)
                if not success and not self.args.force_qc:
                    raise RuntimeError("DeepVariant failed")
                elif not success:
                    self.logger.warning("DeepVariant failed but continuing due to --force-qc flag")
                
                # Use standardized naming for output path
                vcf_output = self.base_dir / 'data/variants/deepvariant' / f"{self.pair_id}.deepvariant.vcf.gz"
                    
            elif self.environment == 'deepsomatic':
                # Run DeepSomatic with full pair_id
                deepsomatic_script = scripts_dir / 'run_deepsomatic_analysis.sh'
                self.logger.info(f"Running DeepSomatic for {self.pair_id}")
                success = self._run_variant_calling_script(str(deepsomatic_script), self.pair_id)
                if not success and not self.args.force_qc:
                    raise RuntimeError("DeepSomatic failed")
                elif not success:
                    self.logger.warning("DeepSomatic failed but continuing due to --force-qc flag")
                
                # Use standardized naming for output path
                vcf_output = self.base_dir / 'data/variants/deepsomatic' / f"{self.pair_id}.deepsomatic.vcf.gz"
                
            else:
                raise ValueError(f"Unknown environment: {self.environment}. Must be 'deepvariant' or 'deepsomatic'")
            
            if vcf_output.exists():
                self.logger.info(f"Variant calling completed: {vcf_output}")
                return str(vcf_output)
            else:
                raise RuntimeError("Variant calling did not produce expected output")
            
        except Exception as e:
            self.logger.error(f"Variant calling failed: {str(e)}")
            if self.args.force_qc:
                self.logger.warning("Continuing despite variant calling failure due to --force-qc flag")
                return None
            raise

    def _verify_file(self, filepath, expected_min_size=None):
        """Verify file exists and optionally check minimum size."""
        path = Path(filepath)
        if not path.exists():
            self.logger.error(f"File does not exist: {filepath}")
            return False
            
        size_mb = path.stat().st_size / (1024**2)
        self.logger.info(f"File {filepath} exists with size: {size_mb:.2f} MB")
        
        if expected_min_size and size_mb < expected_min_size:
            self.logger.warning(f"File {filepath} is smaller than expected: {size_mb:.2f} MB < {expected_min_size} MB")
            return False
            
        return True

    def _log_system_resources(self):
        """Log system resource usage information."""
        try:
            # Log disk space
            disk_usage = shutil.disk_usage(self.base_dir)
            disk_free_gb = disk_usage.free / (1024**3)
            disk_total_gb = disk_usage.total / (1024**3)
            disk_percent = 100 * (disk_usage.total - disk_usage.free) / disk_usage.total
            
            self.logger.info(f"Disk space: {disk_free_gb:.2f} GB free of {disk_total_gb:.2f} GB total ({disk_percent:.1f}% used)")
            
            # Memory info if on Linux
            if os.path.exists('/proc/meminfo'):
                with open('/proc/meminfo', 'r') as f:
                    meminfo = f.read()
                
                # Extract memory information
                mem_total = int(meminfo.split('MemTotal:')[1].split('kB')[0].strip()) / 1024**2
                mem_free = int(meminfo.split('MemFree:')[1].split('kB')[0].strip()) / 1024**2
                mem_available = int(meminfo.split('MemAvailable:')[1].split('kB')[0].strip()) / 1024**2
                
                self.logger.info(f"Memory: {mem_available:.2f} GB available of {mem_total:.2f} GB total")
        except Exception as e:
            self.logger.warning(f"Failed to log system resources: {str(e)}")

    def run(self):
        """Main pipeline execution with enhanced validation and consistent naming."""
        try:
            self.monitor.start_pipeline(self.pair_id)
            self._log_system_resources()        
            # Process input data
            if self.input_fastq:
                self.logger.info("Starting FASTQ processing workflow")
                
                # Download FASTQs if URLs provided
                if all(url.startswith('http') for url in self.input_fastq):
                    self.logger.info(f"Downloading FASTQs from URLs: {self.input_fastq}")
                    fastq_results = self.downloader.process_fastqs(self.input_fastq)
                    self.logger.info(f"Download results: {fastq_results}")
                    self.input_fastq = [str(path) for path in fastq_results.values()]
                    self.logger.info(f"Using downloaded FASTQs: {self.input_fastq}")
                
                # Process FASTQs to BAM using standardized naming
                self.logger.info(f"Starting FASTQ to BAM conversion with files: {self.input_fastq}")
                # Add file size logging
                for fastq in self.input_fastq:
                    if os.path.exists(fastq):
                        size_gb = os.path.getsize(fastq) / (1024**3)
                        self.logger.info(f"FASTQ file {fastq} exists, size: {size_gb:.2f} GB")
                    else:
                        self.logger.warning(f"FASTQ file {fastq} does not exist")
                        
                marked_bam = process_fastqs(
                    self.config,
                    self.input_fastq, 
                    self.pair_id,
                    force=self.force
                )
                self.logger.info(f"FASTQ to BAM conversion completed. BAM file: {marked_bam}")
                
                
                # Run BAM QC before variant calling
                if not self.bam_qc_runner.run_bam_qc(
                    marked_bam,
                    self.names['sample_name'],
                    self._get_qc_thresholds()
                )['passed'] and not self.force_qc:
                    raise QCGateException(f"BAM QC failed for {marked_bam}")

            # In the run() method, modify the "else" branch where it handles input BAM files:
            else:
                self.logger.info("Starting with provided BAM file")
                marked_bam = self.input_bam
                
                # BAM index check/creation
                bam_index = f"{marked_bam}.bai"
                if not os.path.exists(bam_index):
                    self.logger.info(f"BAM index not found, creating index for {marked_bam}")
                    try:
                        subprocess.run(["samtools", "index", marked_bam], check=True)
                        self.logger.info(f"Successfully created BAM index: {bam_index}")
                    except subprocess.CalledProcessError as e:
                        self.logger.error(f"Failed to create BAM index: {str(e)}")
                        if not self.args.force_qc:
                            raise RuntimeError(f"Failed to create BAM index for {marked_bam}")
                
            
            # Run variant calling
            self._log_system_resources()
            vcf_output = self.run_variant_calling(marked_bam)
            
            # Run VCF QC on output
            if vcf_output and not run_vcf_qc(vcf_output) and not self.force_qc:
                raise QCGateException(f"VCF QC failed for {vcf_output}")
            self._log_system_resources()  # Log at end            
            self.monitor.end_pipeline(self.pair_id, 'completed')
            
        except QCGateException as e:
            self.logger.error(f"Pipeline failed QC: {str(e)}")
            self.monitor.end_pipeline(self.pair_id, 'failed_qc')
            raise
        except Exception as e:
            self.logger.error(f"Pipeline failed: {str(e)}")
            self.logger.error(f"Traceback: {traceback.format_exc()}")
            self._log_system_resources()  # Log on failure
            self.monitor.end_pipeline(self.pair_id, 'failed')
            raise

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='WGS Pipeline for ENCODE analysis')
    parser.add_argument('--pair-id', required=True, help='Sample pair ID')
    parser.add_argument('--input-fastq', nargs='+', help='Input FASTQ files or URLs')
    parser.add_argument('--input-bam', help='Input BAM file')
    parser.add_argument('--config', required=True, help='Pipeline configuration file')
    parser.add_argument('--library-type', default='pcr', choices=['pcr', 'pcr_free'], help='Library type')
    parser.add_argument('--force', action='store_true', help='Force overwrite existing files')
    parser.add_argument('--force-qc', action='store_true', help='Continue despite QC failures')
    parser.add_argument('--environment', required=True, choices=['deepvariant', 'deepsomatic'],
                      help='Specify which variant caller environment is available')
    
    args = parser.parse_args()
    
    pipeline = WGSPipeline(args)
    pipeline.run()