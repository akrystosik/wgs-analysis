from typing import Dict, List, Optional, Any    
import logging
from pathlib import Path  # Add this line
import subprocess  # Add this for running commands
import shutil     # Add this for file operations

class LocalDeepVariantRunner:
    def __init__(self, pipeline_config: Dict[str, Any]):
        self.config = pipeline_config
        self.logger = logging.getLogger(__name__)
        
    def run_deepvariant(self,
                       pair_id: str,
                       bam_file: str,
                       reference_path: str,
                       output_vcf: str,
                       output_gvcf: str,
                       num_shards: int = 128) -> bool:
        try:
            # Verify inputs
            for f in [bam_file, reference_path]:
                if not Path(f).exists():
                    raise FileNotFoundError(f"File not found: {f}")

            # Setup directories
            output_dir = Path(output_vcf).parent
            output_dir.mkdir(parents=True, exist_ok=True)
            intermediate_dir = output_dir / f'intermediate_{pair_id}'
            if intermediate_dir.exists():
                shutil.rmtree(intermediate_dir)
            intermediate_dir.mkdir()
            
            # Run DeepVariant
            cmd = [
                "/opt/deepvariant/bin/run_deepvariant",
                f"--model_type=WGS",
                f"--ref={reference_path}",
                f"--reads={bam_file}",
                f"--output_vcf={output_vcf}",
                f"--output_gvcf={output_gvcf}",
                f"--intermediate_results_dir={str(intermediate_dir)}",
                f"--num_shards={num_shards}",
                "--verbosity=2"
            ]
            
            self.logger.info(f"Starting DeepVariant for {pair_id}")
            self.logger.info(f"Command: {' '.join(cmd)}")
            
            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True
            )
            
            # Monitor progress
            while True:
                output = process.stdout.readline()
                if output == '' and process.poll() is not None:
                    break
                if output:
                    self.logger.info(output.strip())
            
            rc = process.poll()
            if rc != 0:
                _, stderr = process.communicate()
                self.logger.error(f"DeepVariant failed: {stderr}")
                return False
                
            # Index VCF
            index_cmd = ["bcftools", "index", "-t", output_vcf]
            result = subprocess.run(index_cmd, capture_output=True, text=True)
            if result.returncode != 0:
                self.logger.error(f"VCF indexing failed: {result.stderr}")
                return False
            
            # Create filtered VCF
            filtered_vcf = str(output_vcf).replace('.vcf.gz', '.filtered.vcf.gz')
            filter_cmd = [
                "bcftools", "filter",
                "-i", f"QUAL>={self.config['vcf_filters']['min_quality']}",
                "-Oz", "-o", filtered_vcf,
                output_vcf
            ]
            result = subprocess.run(filter_cmd, capture_output=True, text=True)
            if result.returncode != 0:
                self.logger.error(f"VCF filtering failed: {result.stderr}")
                return False
                
            # Index filtered VCF
            index_cmd = ["bcftools", "index", "-t", filtered_vcf]
            result = subprocess.run(index_cmd, capture_output=True, text=True)
            if result.returncode != 0:
                self.logger.error(f"Filtered VCF indexing failed: {result.stderr}")
                return False
            
            self.logger.info("DeepVariant completed successfully")
            return True
            
        except Exception as e:
            self.logger.error(f"DeepVariant execution failed: {str(e)}")
            raise
        finally:
            # Cleanup
            if 'intermediate_dir' in locals() and intermediate_dir.exists():
                try:
                    shutil.rmtree(intermediate_dir)
                except Exception as e:
                    self.logger.warning(f"Cleanup failed: {str(e)}")