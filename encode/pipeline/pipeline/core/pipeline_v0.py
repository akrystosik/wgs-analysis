#!/usr/bin/env python3
#encode_analysis/pipeline_v1_fasta_bam/core/pipeline.py

import subprocess
import logging
import os
from pathlib import Path
import psutil
from concurrent.futures import ThreadPoolExecutor
import shutil             
import time               
from datetime import datetime
import traceback 
from typing import Dict, List, Optional, Any

# Core functionality
from pipeline.core.monitoring import WGSMonitor
from pipeline.core.exceptions import QCGateException

# Utilities
from pipeline.utils.helpers import (
    setup_directories,
    run_command,
    parse_bam_metrics,
    parse_vcf_metrics,
    generate_html_report,
    load_config, 
    setup_logging
)

# QC modules
from pipeline.qc.bam import BAMQCRunner
from pipeline.qc.fastq import FastQMultiQC
from pipeline.qc.vcf import run_vcf_qc

# IO modules
from pipeline.io.fastq import process_fastqs
from pipeline.io.downloader import FastqDownloader
from pipeline.io.verify import verify_fastq_pairs



class WGSPipeline:
    def __init__(self, args):
        self.args = args  # Assign the entire args object
        self.config = load_config(args.config)
        self.pair_id = args.pair_id
        self.force = args.force
        self.force_qc = args.force_qc
        self.input_fastq = args.input_fastq
        self.input_bam = args.input_bam
        
        self.datetime = datetime  
        setup_logging(self.config)
        self.logger = logging.getLogger(__name__)
        
        self.logger.debug("DEBUG: Starting pipeline.py")
        self.logger.debug(f"DEBUG: datetime imported: {self.datetime}")
        self.logger.debug("DEBUG: Imports completed")
        
        self.base_dir = Path(self.config['base_dir'])
        self.reference_dir = Path(self.base_dir / self.config['data_dirs']['reference'])
        
        self.reference_path = Path(self.reference_dir / self.config['reference']['fasta'])
        
        self.logger.info(f"Using reference directory: {self.reference_dir}")
        self.logger.info(f"Using reference genome: {self.reference_path}")
        self.library_type = self._determine_library_type()
        
        
        # Verify reference path exists
        if not self.reference_path.exists():
            raise FileNotFoundError(
                f"Reference genome not found: {self.reference_path}\n"
                f"Please ensure the reference genome is located at this path."
            )
        
        # Continue with other initialization
        self.force_qc = self.args.force_qc
        self.monitor = WGSMonitor(self.config)
        setup_directories(self.config)
        self.downloader = FastqDownloader(self.config)
        self.fastqc_runner = FastQMultiQC(self.config)

    def _determine_library_type(self) -> str:
        """Determine library type from metadata or arguments."""
        if hasattr(self.args, 'library_type'):
            return self.args.library_type
        # Default to more stringent PCR settings if unknown
        return "pcr"
  

    def _verify_reference(self, reference_path: str = None) -> bool:
        """Verify reference genome and indices exist with improved path handling."""
        try:
            # Use class reference path if none provided
            reference_path = reference_path or self.reference_path
            reference_path = Path(reference_path)
            
            self.logger.info(f"Verifying reference genome: {reference_path}")
            
            # Check main reference file
            if not reference_path.exists():
                self.logger.error(f"Reference genome not found: {reference_path}")
                self.logger.info("Expected directory structure:")
                self.logger.info(f"  Base directory: {self.base_dir}")
                self.logger.info(f"  Reference directory: {self.reference_dir}")
                self.logger.info(f"  Reference file: {self.config['reference']['fasta']}")
                return False
            
            # Check BWA indices with proper path handling
            required_indices = ['.amb', '.ann', '.bwt', '.pac', '.sa']
            missing_indices = []
            
            for idx in required_indices:
                idx_path = reference_path.with_suffix(reference_path.suffix + idx)
                if not idx_path.exists():
                    missing_indices.append(idx)
                    self.logger.error(f"Missing BWA index: {idx_path}")
            
            if missing_indices:
                self.logger.error(f"Missing indices: {', '.join(missing_indices)}")
                
                # Attempt to build missing indices
                if self.args.force or input("Would you like to build missing indices? [y/N] ").lower() == 'y':
                    if self._build_bwa_indices(str(reference_path)):
                        self.logger.info("Successfully built missing BWA indices")
                        return True
                return False
                
            self.logger.info("Reference genome and all indices verified successfully")
            return True
            
        except Exception as e:
            self.logger.error(f"Reference verification failed: {str(e)}")
            return False

    def _build_bwa_indices(self, reference_path: str) -> bool:
        """Build BWA indices with improved error handling."""
        try:
            self.logger.info(f"Building BWA indices for {reference_path}")
            
            # Verify BWA is installed
            if subprocess.run(['which', 'bwa'], capture_output=True).returncode != 0:
                self.logger.error("BWA not found in PATH. Please install BWA first.")
                return False
            
            # Build indices
            cmd = f"bwa index {reference_path}"
            process = subprocess.Popen(
                cmd,
                shell=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True
            )
            
            # Monitor index building progress
            while True:
                output = process.stdout.readline()
                if output == '' and process.poll() is not None:
                    break
                if output:
                    self.logger.info(output.strip())
            
            rc = process.poll()
            if rc != 0:
                _, stderr = process.communicate()
                self.logger.error(f"BWA index build failed: {stderr}")
                return False
                
            self.logger.info("BWA indices built successfully")
            return True
            
        except Exception as e:
            self.logger.error(f"BWA index building failed: {str(e)}")
            return False

    def _monitor_resources(self, process_name):
        """Monitor system resources during processing."""
        while True:
            try:
                cpu_percent = psutil.cpu_percent(interval=1, percpu=True)
                mem = psutil.virtual_memory()
                
                logging.info(
                    f"Resource Usage ({process_name}) - "
                    f"CPU Usage: {sum(cpu_percent)/len(cpu_percent):.1f}% "
                    f"Memory Used: {mem.percent}% "
                    f"({mem.used/1024/1024/1024:.1f}GB / {mem.total/1024/1024/1024:.1f}GB)"
                )
                
                time.sleep(30)  # Update every 30 seconds
            except Exception as e:
                logging.error(f"Resource monitoring error: {str(e)}")
                break
                    
    def validate_inputs(self):
        """Validate all input parameters and paths."""
        if not self.args.pair_id:
            raise ValueError("--pair-id is required")
            
        if not self.args.input_bam and not self.args.input_fastq:
            raise ValueError("Either --input-bam or --input-fastq must be provided")
            
        # Validate input paths exist
        if self.args.input_bam and not os.path.exists(self.args.input_bam):
            raise FileNotFoundError(f"Input BAM file not found: {self.args.input_bam}")
            
        if self.args.input_fastq:
            for fastq in self.args.input_fastq:
                if not (fastq.startswith('http') or os.path.exists(fastq)):
                    raise FileNotFoundError(f"FASTQ file not found or invalid URL: {fastq}")


    def _process_fastq_file(self, fastq_file, qc_dir):
        """Process a single FASTQ file with FastQC."""
        try:
            self.logger.debug("DEBUG: Entering _process_fastq_file")
            self.logger.debug(f"DEBUG: datetime available in method: {self.datetime}")  # Corrected to use self.datetime

            # Calculate threads_per_process based on configuration
            total_threads = self.config['tools']['fastqc']['threads']
            max_processes = self.config['tools']['fastqc']['execution']['max_processes']
            threads_per_process = total_threads // max_processes
            
            start_time = time.time()
            file_size = os.path.getsize(fastq_file) / (1024 * 1024 * 1024)
            base_name = fastq_file.name.split('.')[0]
            
            self.logger.info(f"Starting FastQC for {fastq_file.name}")
            self.logger.info(f"File size: {file_size:.2f}GB")
            self.logger.info(f"Using base name: {base_name}")
            
            # Build command
            cmd = [
                "fastqc",
                str(fastq_file),
                f"--outdir={str(qc_dir)}",
                f"--threads={threads_per_process}",
                "--noextract"  # Removed --quiet to see output
            ]
            
            self.logger.info(f"Running command: {' '.join(cmd)}")
            
            # Run FastQC with output monitoring
            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
                bufsize=1
            )
            
            # Monitor output in real-time
            while True:
                stdout_line = process.stdout.readline()
                stderr_line = process.stderr.readline()
                
                if stdout_line:
                    self.logger.info(f"FastQC stdout: {stdout_line.strip()}")
                if stderr_line:
                    self.logger.info(f"FastQC stderr: {stderr_line.strip()}")
                
                if stdout_line == '' and stderr_line == '' and process.poll() is not None:
                    break

            # Get final return code
            return_code = process.poll()
            if return_code != 0:
                self.logger.error(f"FastQC process returned {return_code}")
                stdout, stderr = process.communicate()
                self.logger.error(f"Final stdout: {stdout}")
                self.logger.error(f"Final stderr: {stderr}")
                return False
            
            # Verify outputs
            fastqc_html = Path(qc_dir / f"{base_name}_fastqc.html")
            fastqc_zip = Path(qc_dir / f"{base_name}_fastqc.zip")
            
            self.logger.info("Checking FastQC outputs:")
            self.logger.info(f"  HTML path: {fastqc_html}")
            self.logger.info(f"  ZIP path: {fastqc_zip}")
            self.logger.info(f"  HTML exists: {fastqc_html.exists()}")
            self.logger.info(f"  ZIP exists: {fastqc_zip.exists()}")
            
            # List directory contents
            self.logger.info("Directory contents:")
            for f in qc_dir.glob('*'):
                self.logger.info(f"  {f}")
            
            if not (fastqc_html.exists() and fastqc_zip.exists()):
                self.logger.error(f"FastQC outputs not found for {fastq_file.name}")
                return False
            
            elapsed_time = (time.time() - start_time) / 60
            self.logger.info(f"Completed {fastq_file.name} in {elapsed_time:.2f} minutes")
            return True
                
        except Exception as e:
            self.logger.error(f"Error processing {fastq_file}: {str(e)}")
            self.logger.error(f"Exception type: {type(e)}")
            self.logger.error(f"Traceback: {traceback.format_exc()}")
            return False


    def run_fastqc(self, fastq_files: List[Path]) -> bool:
        """Run FastQC on a list of FASTQ files with optimized settings."""
        self.logger.info(f"Running FastQC on {len(fastq_files)} files")
        threads_per_process = max(1, self.config['tools']['fastqc']['threads'] // self.config['tools']['fastqc']['execution']['max_processes'])
        java_opts = self._build_java_opts()
        
        def process_single_fastq(fastq_file: Path) -> bool:
            try:
                cmd = [
                    "fastqc",
                    str(fastq_file),
                    f"--outdir={str(self.fastqc_dir)}",
                    f"--threads={threads_per_process}",
                    "--noextract",
                    "--nogroup"     # Removed --quiet to see output
                ]
                
                self.logger.info(f"Processing FastQC for {fastq_file}")
                self.logger.debug(f"Running command: {' '.join(cmd)} with JAVA_TOOL_OPTIONS: {' '.join(java_opts)}")
                
                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    check=True,
                    env={**os.environ, 'JAVA_TOOL_OPTIONS': ' '.join(java_opts)}
                )
                
                self.logger.debug(f"FastQC stdout: {result.stdout}")
                self.logger.debug(f"FastQC stderr: {result.stderr}")
                
                expected_output = self.fastqc_dir / f"{fastq_file.stem}_fastqc.html"
                if not expected_output.exists():
                    self.logger.error(f"FastQC output not found for {fastq_file}")
                    return False
                    
                self.logger.debug(f"FastQC completed for {fastq_file}")
                return True
                    
            except subprocess.CalledProcessError as e:
                self.logger.error(f"FastQC failed for {fastq_file}: {e.stderr}")
                return False
            except Exception as e:
                self.logger.error(f"Unexpected error processing {fastq_file}: {str(e)}")
                return False

        try:
            max_workers = min(len(fastq_files), os.cpu_count() or 1)
            success = True
            
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                futures = {
                    executor.submit(process_single_fastq, fastq_file): fastq_file 
                    for fastq_file in fastq_files
                }
                
                for future in as_completed(futures):
                    fastq_file = futures[future]
                    try:
                        if not future.result():
                            success = False
                            if not getattr(self, 'force_qc', False):
                                self.logger.warning(f"FastQC failed for {fastq_file}. Stopping execution.")
                                executor.shutdown(cancel_futures=True)
                                break
                    except Exception as e:
                        success = False
                        self.logger.error(f"FastQC execution failed for {fastq_file}: {str(e)}")
                        if not getattr(self, 'force_qc', False):
                            executor.shutdown(cancel_futures=True)
                            break
            
            if success and self.config.get('tools', {}).get('fastqc', {}).get('cleanup_temp', True):
                self._cleanup_temp_files()
                
            return success
            
        except Exception as e:
            self.logger.error(f"Failed to run FastQC processing: {str(e)}")
            return False


    def _validate_bam(self, step_name: str, bam_file: str) -> bool:
        """Validate BAM file with QC metrics."""
        try:
            self.logger.info(f"Running BAM QC validation for step: {step_name}")
            
            # Run QC using run_bam_qc
            qc_runner = BAMQCRunner(self.config)
            metrics, evaluation = qc_runner.run_bam_qc(bam_file, self.args.pair_id, step_name)

            # Add debug statements
            self.logger.debug(f"Metrics: {metrics}")
            self.logger.debug(f"Evaluation: {evaluation}")

            if not evaluation['passed'] and not self.args.force_qc:
                failures = '\n'.join(evaluation['failures'])
                raise QCGateException(
                    f"BAM QC failed at {step_name}:\n{failures}"
                )
            elif not evaluation['passed']:
                self.logger.warning(
                    f"BAM QC failed but continuing due to --force-qc flag at {step_name}"
                )
            
            return evaluation['passed']
            
        except Exception as e:
            self.logger.error(f"BAM validation failed: {str(e)}")
            if self.args.force_qc:
                self.logger.warning("Continuing despite BAM validation failure due to --force-qc flag")
                return True
            raise    


    def process_bam_data(self, input_bam: str):
        """Process BAM file using qc.bam.py's BAMQCRunner."""
        self.monitor.start_step('bam_processing', self.args.pair_id)
        
        try:
            runner = BAMQCRunner(self.config)
            
            # Step 1: Process BAM (sort by queryname, fixmate, sort by coordinate, markdup)
            marked_bam = runner.process_bam(original_bam=input_bam, pair_id=self.args.pair_id, force=self.args.force)
            
            # Step 2: Run QC on the marked BAM using run_bam_qc
            metrics, evaluation = runner.run_bam_qc(marked_bam, self.args.pair_id, step='post_markdup')
            
            # Add debug statements
            self.logger.debug(f"Metrics: {metrics}")
            self.logger.debug(f"Evaluation: {evaluation}")

            if not evaluation['passed'] and not self.args.force_qc:
                self.logger.error("BAM processing QC failed.")
                raise RuntimeError("BAM processing QC failed.")
            
            self.monitor.end_step('bam_processing', self.args.pair_id, status='completed')
            return marked_bam
            
        except Exception as e:
            self.monitor.end_step('bam_processing', self.args.pair_id, status='failed')
            self.logger.error(f"BAM processing failed: {str(e)}")
            raise


    def _verify_bam_for_variant_calling(self, bam_file: str) -> bool:
        """Verify BAM file is ready for variant calling."""
        try:
            # Check BAM file exists
            if not os.path.exists(bam_file):
                self.logger.error(f"BAM file not found: {bam_file}")
                return False
                
            # Check index exists
            if not os.path.exists(f"{bam_file}.bai"):
                self.logger.info(f"Creating index for {bam_file}")
                result = subprocess.run(['samtools', 'index', bam_file], 
                                    capture_output=True, text=True)
                if result.returncode != 0:
                    self.logger.error(f"Failed to create BAM index: {result.stderr}")
                    return False
            
            # Verify BAM is sorted
            header = subprocess.run(['samtools', 'view', '-H', bam_file],
                                capture_output=True, text=True)
            if 'SO:coordinate' not in header.stdout:
                self.logger.error(f"BAM file is not coordinate sorted: {bam_file}")
                return False
                
            self.logger.info(f"BAM file {bam_file} is ready for variant calling")
            return True
                
        except Exception as e:
            self.logger.error(f"Error verifying BAM file: {str(e)}")
            return False


    def _determine_library_type(self) -> str:
        """Determine library type from metadata or arguments."""
        if hasattr(self.args, 'library_type'):
            return self.args.library_type
        # Default to more stringent PCR settings if unknown
        return "pcr"
    
    def get_variant_calling_params(self) -> Dict:
        """Get library-specific variant calling parameters."""
        base_params = self.config['tools']['deepvariant']['base_settings']
        library_params = self.config['library_types'][self.library_type]['variant_calling']['deepvariant']
        return {**base_params, **library_params}
    
    def get_qc_thresholds(self) -> Dict:
        """Get library-specific QC thresholds."""
        return self.config['library_types'][self.library_type]['qc']['thresholds']
    
    def run_qc(self, bam_path: str) -> Dict:
        """Run QC with library-specific metrics and thresholds."""
        qc_config = self.config['library_types'][self.library_type]['qc']
        thresholds = self.get_qc_thresholds()
        
        metrics = {}
        for metric in qc_config['metrics']['required']:
            metrics[metric] = self._run_qc_metric(
                metric=metric,
                bam_path=bam_path,
                threshold=thresholds.get(metric, {})
            )
        
        return metrics
    
    def generate_qc_report(self, metrics: Dict) -> str:
        """Generate QC report with library-specific sections and plots."""
        report_config = self.config['qc_reporting'][self.library_type]
        thresholds = self.get_qc_thresholds()
        
        report = QCReport(
            metrics=metrics,
            thresholds=thresholds,
            sections=report_config['report_sections'],
            plots=report_config['plots']
        )
        
        return report.generate()

    def run_variant_calling(self, marked_bam: str) -> str:
        """Run variant calling with appropriate caller based on sample type."""
        try:
            vcf_output = self.variants_dir / f"{self.pair_id}.vcf.gz"
            gvcf_output = self.variants_dir / f"{self.pair_id}.g.vcf.gz"

            if hasattr(self.args, 'is_cancer') and self.args.is_cancer:
                # Use DeepSomatic for cancer samples
                from pipeline.core.run_deepsomatic import DeepSomaticRunner
                runner = DeepSomaticRunner(self.config)
                self.logger.info("Running DeepSomatic for cancer sample")
                success = runner.run_somatic_calling(
                    tumor_bam=marked_bam,
                    pair_id=self.pair_id
                )
            else:
                # Use DeepVariant for non-cancer samples
                from pipeline.core.local_runner import LocalDeepVariantRunner
                runner = LocalDeepVariantRunner(self.config)
                self.logger.info("Running DeepVariant for non-cancer sample")
                success = runner.run_deepvariant(
                    pair_id=self.pair_id,
                    bam_file=marked_bam,
                    reference_path=str(self.reference_path),
                    output_vcf=str(vcf_output),
                    output_gvcf=str(gvcf_output)
                )

            if not success and not self.args.force_qc:
                raise RuntimeError("Variant calling failed")
            elif not success:
                self.logger.warning("Variant calling failed but continuing due to --force-qc flag")

            # Verify and process output
            if vcf_output.exists():
                self.logger.info(f"Variant calling completed: {vcf_output}")
                
                # Index VCF
                index_cmd = ["bcftools", "index", "-t", str(vcf_output)]
                if subprocess.run(index_cmd, check=True).returncode != 0:
                    raise RuntimeError("VCF indexing failed")
                
                # Filter VCF
                filtered_vcf = self.filter_vcf(str(vcf_output))
                if filtered_vcf:
                    self.logger.info(f"Generated filtered VCF: {filtered_vcf}")
                    return filtered_vcf
                
            raise RuntimeError("Variant calling did not produce expected output")

        except Exception as e:
            self.logger.error(f"Variant calling failed: {str(e)}")
            if self.args.force_qc:
                self.logger.warning("Continuing despite variant calling failure due to --force-qc flag")
                return None
            raise

    def filter_vcf(self, vcf_path: str) -> str:
        """Filter VCF file based on criteria defined in config.yaml."""
        self.logger.info(f"Filtering VCF file: {vcf_path}")
        try:
            vcf_filters = self.config.get('vcf_filters', {})
            min_depth = vcf_filters.get('min_depth', 10)
            min_quality = vcf_filters.get('min_quality', 20)
            ts_tv_ratio_range = vcf_filters.get('ts_tv_ratio', [2.0, 2.1])
            het_hom_ratio_range = vcf_filters.get('het_hom_ratio', [1.5, 2.1])
            
            # Convert vcf_path to Path object
            vcf_path = Path(vcf_path)
            
            # Construct bcftools filter expression
            filter_expr = (
                f"QUAL >= {min_quality} && FORMAT/DP >= {min_depth} && "
                f"TS_TV >= {ts_tv_ratio_range[0]} && TS_TV <= {ts_tv_ratio_range[1]} && "
                f"HET_HOM >= {het_hom_ratio_range[0]} && HET_HOM <= {het_hom_ratio_range[1]}"
            )
            filtered_vcf = vcf_path.with_suffix(".filtered.vcf.gz")
            filter_cmd = [
                "bcftools", "filter",
                "-i", filter_expr,
                "-Oz",
                "-o", str(filtered_vcf),
                str(vcf_path)
            ]
            
            result = subprocess.run(filter_cmd, capture_output=True, text=True)
            if result.returncode != 0:
                self.logger.error(f"VCF filtering failed: {result.stderr}")
                raise RuntimeError(f"VCF filtering failed: {result.stderr}")
            
            # Index the filtered VCF
            index_cmd = ["bcftools", "index", str(filtered_vcf)]
            result = subprocess.run(index_cmd, capture_output=True, text=True)
            if result.returncode != 0:
                self.logger.error(f"Failed to index filtered VCF: {result.stderr}")
                raise RuntimeError(f"VCF indexing failed: {result.stderr}")
            
            return str(filtered_vcf)
            
        except Exception as e:
            self.logger.error(f"VCF filtering failed: {str(e)}")
            if not self.args.force_qc:
                raise
            return str(vcf_path)


    def generate_final_report(self):
        """Generate comprehensive QC report."""
        report_dir = Path(self.base_dir) / self.config['data_dirs']['qc'] / 'reports'
        report_dir.mkdir(parents=True, exist_ok=True)
        
        try:
            # Collect all metrics
            metrics = {}
            
            # FASTQ Metrics
            try:
                metrics['fastq'] = self.monitor.get_fastqc_metrics(self.config, self.args.pair_id)
            except Exception as e:
                self.logger.error(f"Failed to retrieve FastQC metrics: {str(e)}")
                metrics['fastq'] = {"error": "FastQC metrics not available."}
            
            # Alignment Metrics
            try:
                metrics['alignment'] = self.monitor.get_bam_qc_metrics(
                    self.config, 
                    self.args.pair_id, 
                    'post_markdup'
                )
            except Exception as e:
                self.logger.error(f"Failed to retrieve BAM QC metrics: {str(e)}")
                metrics['alignment'] = {"error": "BAM QC metrics not available."}
            
            # Variant Metrics
            try:
                vcf_path = Path(self.base_dir) / self.config['data_dirs']['variants'] / f"{self.args.pair_id}.filtered.vcf.gz"
                self.logger.info(f"Checking for VCF at: {vcf_path}")
                self.logger.info(f"Base dir: {self.base_dir}")
                self.logger.info(f"Variants dir: {self.config['data_dirs']['variants']}")                
                
                if vcf_path.exists():
                    self.logger.info(f"Found VCF file at: {vcf_path}")
                    metrics['variants'] = parse_vcf_metrics(str(vcf_path.resolve()))
                else:
                    self.logger.error(f"VCF file not found at: {vcf_path}")
                    metrics['variants'] = {"error": f"VCF file not found at: {vcf_path}"}
            except Exception as e:
                self.logger.error(f"Failed to retrieve VCF metrics: {str(e)}")
                metrics['variants'] = {"error": f"Variant metrics error: {str(e)}"}
            
            # Generate HTML report
            report_path = report_dir / f"{self.args.pair_id}_final_report.html"
            report_content = generate_html_report(metrics, self.args.pair_id)
            
            with open(report_path, 'w') as f:
                f.write(report_content)
                
            self.logger.info(f"Final report generated: {report_path}")
            
        except Exception as e:
            self.logger.error(f"Report generation encountered an error: {str(e)}")
            if not self.args.force_qc:
                raise
            self.logger.warning("Report generated with partial data due to errors.")


    def run(self):
        """Main pipeline execution with comprehensive error handling."""
        try:
            # Validate inputs
            self.validate_inputs()
            self.monitor.start_pipeline(self.pair_id)
            
            # Process input data
            if self.input_fastq:
                try:
                    self.logger.info("Starting FASTQ processing workflow")
                    input_bam = process_fastqs(
                        self.config,
                        self.input_fastq, 
                        self.pair_id,
                        force=self.force
                    )
                    
                    self.logger.info("Starting BAM processing workflow")
                    marked_bam = self.process_bam_data(input_bam)
                    
                except Exception as e:
                    self.logger.error(f"FASTQ/BAM processing failed: {str(e)}")
                    if not self.force_qc:
                        raise
                    self.logger.warning("Continuing despite error due to --force-qc flag")
                    marked_bam = input_bam
                    
            else:
                self.logger.info("Starting with provided BAM file")
                self._validate_bam("input_bam", self.input_bam)
                marked_bam = self.input_bam
                
            try:
                self.logger.info("Starting variant calling workflow")
                vcf_output = self.run_variant_calling(marked_bam)
                
                if vcf_output:
                    self.logger.info(f"Variant calling completed. Filtered VCF: {vcf_output}")
                else:
                    self.logger.warning("Variant calling was skipped or failed.")
                
            except Exception as e:
                self.logger.error(f"Variant calling failed: {str(e)}")
                if not self.force_qc:
                    raise
                self.logger.warning("Skipping variant calling due to error and --force-qc flag")
                vcf_output = None
                
        except QCGateException as e:
            self.logger.error(f"Pipeline failed QC gates: {str(e)}")
            self.monitor.end_pipeline(self.pair_id, status='failed_qc')
        except Exception as e:
            self.logger.error(f"Pipeline failed: {str(e)}")
            self.monitor.end_pipeline(self.pair_id, status='failed')
        finally:
            try:
                if vcf_output:
                    self.logger.info(f"Final VCF output: {vcf_output}")
                else:
                    self.logger.warning("Variant calling was skipped or failed. Final VCF output is unavailable.")
                
                self.logger.info("Generating final analysis report")
                self.generate_final_report()
                
                if hasattr(self, 'pair_id'):
                    self._cleanup_temp_files(self.pair_id)
                    
                if vcf_output:
                    self.logger.info(f"Pipeline completed successfully for {self.pair_id}")
                    self.monitor.end_pipeline(self.pair_id, 'completed')
                else:
                    self.logger.warning(f"Pipeline completed with warnings for {self.pair_id}")
                    self.monitor.end_pipeline(self.pair_id, 'completed_with_warnings')
                    
            except Exception as e:
                self.logger.error(f"Cleanup/finalization failed: {str(e)}")

    def _cleanup_temp_files(self, pair_id: str):
        """Clean up temporary files."""
        try:
            for cleanup_dir in ['variants', 'bam', 'fastq']:
                cleanup_path = self.base_dir / f'data/{cleanup_dir}'
                if cleanup_path.exists():
                    for pattern in [f"{pair_id}*.tmp", f"{pair_id}*.temp", f"tmp.{pair_id}*"]:
                        for temp_file in cleanup_path.glob(pattern):
                            try:
                                temp_file.unlink()
                                self.logger.debug(f"Removed temp file: {temp_file}")
                            except Exception as e:
                                self.logger.warning(f"Failed to remove {temp_file}: {str(e)}")
        except Exception as e:
            self.logger.warning(f"Cleanup failed for {pair_id}: {str(e)}")