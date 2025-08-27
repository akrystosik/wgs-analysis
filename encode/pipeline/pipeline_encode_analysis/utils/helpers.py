# In pipeline/utils/helpers.py

from typing import List, Any

import logging
import os
import requests
import subprocess
import yaml  # Added import for YAML
from pathlib import Path
from typing import Dict, List
from datetime import datetime  # Added import for datetime

logger = logging.getLogger(__name__)

def run_command(command: str, description: str = ""):
    """
    Execute a shell command and handle errors.
    
    Args:
        command (str): The shell command to execute.
        description (str): Description of the command for logging purposes.
    
    Raises:
        subprocess.CalledProcessError: If the command fails.
    """
    logger.debug(f"Executing command: {command}")
    try:
        result = subprocess.run(
            command, shell=True, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        logger.info(f"{description} completed successfully.")
        logger.debug(f"Command output: {result.stdout}")
    except subprocess.CalledProcessError as e:
        logger.error(f"{description} failed with error: {e.stderr.strip()}")
        raise

def download_file(url: str, dest_path: str, chunk_size: int = 1024*1024):
    """
    Download a file from a URL to a local destination.
    
    Args:
        url (str): URL of the file to download.
        dest_path (str): Local path where the file will be saved.
        chunk_size (int): Chunk size for streaming download.
    
    Raises:
        requests.HTTPError: If the download fails.
    """
    logger.debug(f"Initiating download from {url} to {dest_path}")
    response = requests.get(url, stream=True)
    try:
        response.raise_for_status()
    except requests.HTTPError as e:
        logger.error(f"Failed to download {url}: {e}")
        raise
    
    Path(os.path.dirname(dest_path)).mkdir(parents=True, exist_ok=True)
    
    with open(dest_path, 'wb') as f:
        for chunk in response.iter_content(chunk_size=chunk_size):
            if chunk:  # filter out keep-alive new chunks
                f.write(chunk)
    
    logger.info(f"Download completed: {dest_path}")


def verify_fastq_pairs(fastq_inputs: List[str], fastq_dir: str) -> Dict[str, str]:
    """Verify that FASTQ files are correctly paired (read1 and read2)."""
    logger.debug(f"Verifying FASTQ pairs: {fastq_inputs}")
    
    if len(fastq_inputs) != 2:
        logger.error("Expected exactly two FASTQ files for pairing (read1 and read2).")
        raise ValueError("Expected exactly two FASTQ files for pairing (read1 and read2).")
    
    read1_input, read2_input = fastq_inputs
    fastq_dir = Path(fastq_dir)
    fastq_dir.mkdir(parents=True, exist_ok=True)
    
    def process_fastq_input(input_path: str, read_num: int) -> str:
        if input_path.startswith('https://www.encodeproject.org/'):
            filename = Path(input_path.split('@@download/')[1]).name
            local_path = fastq_dir / filename
            if not local_path.exists():
                logger.info(f"Downloading Read{read_num} from {input_path}")
                # ENCODE requires specific headers
                headers = {'accept': 'application/json'}
                response = requests.get(input_path, headers=headers, stream=True)
                if response.status_code == 403:
                    # Try without headers if 403
                    response = requests.get(input_path, stream=True)
                response.raise_for_status()
                
                with open(local_path, 'wb') as f:
                    for chunk in response.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)
            return str(local_path)
        return input_path

    try:
        read1_path = process_fastq_input(read1_input, 1)
        read2_path = process_fastq_input(read2_input, 2)
        logger.info(f"FASTQ pairs verified:\nRead1: {read1_path}\nRead2: {read2_path}")
        return {'read1': read1_path, 'read2': read2_path}
        
    except Exception as e:
        logger.error(f"Failed to verify/download FASTQ pairs: {str(e)}")
        raise


def setup_logging(config):
    """Setup logging configuration"""
    # Create timestamped log filename
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    log_dir = config.get('logging', {}).get('log_dir', 'logs')
    log_file = Path(log_dir) / f"pipeline_{timestamp}.log"
    
    # Create log directory if it doesn't exist
    Path(log_dir).mkdir(parents=True, exist_ok=True)
    
    # Configure logging
    log_level = config.get('logging', {}).get('log_level', 'INFO')
    logging.basicConfig(
        level=getattr(logging, log_level),
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(str(log_file)),
            logging.StreamHandler()  # This will also print to console
        ],
        force=True  # This ensures logging is reconfigured
    )
    
    # Log initial setup information
    logging.info(f"Starting new pipeline run. Log file: {log_file}")
    
def load_config(config_path):
    """Load configuration from YAML file"""
    with open(config_path) as f:
        return yaml.safe_load(f)

def setup_directories(config):
    """Create necessary directories if they don't exist"""
    for dir_name, dir_path in config['data_dirs'].items():
        path = Path(config['base_dir']) / dir_path
        path.mkdir(parents=True, exist_ok=True)
        logging.info(f"Created/verified directory: {path}")

def parse_bam_metrics(metrics_file):
    """Parse BAM metrics from file"""
    try:
        with open(metrics_file) as f:
            return yaml.safe_load(f)
    except Exception as e:
        logging.error(f"Failed to parse BAM metrics: {str(e)}")
        return {}

def parse_vcf_metrics(vcf_path):
    """Parse VCF metrics from file."""
    try:
        # Ensure absolute path
        vcf_path = Path(vcf_path)
        if not vcf_path.is_absolute():
            vcf_path = Path("/mnt/czi-sci-ai/intrinsic-variation-gene-ex/encode_analysis/data/variants") / vcf_path.name
            
        self.logger.info(f"Using absolute VCF path: {vcf_path}")
        
        if not vcf_path.exists():
            raise FileNotFoundError(f"VCF file not found: {vcf_path}")
            
        cmd = ["bcftools", "stats", str(vcf_path)]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        return parse_vcf_stats(result.stdout)
            
    except Exception as e:
        self.logger.error(f"Failed to parse VCF metrics: {str(e)}")
        return {}
    
    
def parse_vcf_stats(stats_output):
    """Parse bcftools stats output"""
    metrics = {}
    try:
        for line in stats_output.split('\n'):
            if line.startswith('SN'):
                parts = line.split('\t')
                if len(parts) >= 3:
                    key = parts[2].strip(':')
                    if len(parts) > 3:
                        value = parts[3]
                        try:
                            value = float(value)
                        except ValueError:
                            pass
                        metrics[key] = value
        return metrics
    except Exception as e:
        logging.error(f"Failed to parse VCF stats: {str(e)}")
        return {}

def generate_html_report(metrics, pair_id):
    """Generate HTML report from metrics"""
    if not isinstance(metrics, dict):
        metrics = {'error': 'No metrics available'}
    
    html_template = f"""
    <html>
    <head>
        <title>Pipeline Report - {pair_id}</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 20px; }}
            h1 {{ color: #333; }}
            .section {{ margin: 20px 0; padding: 10px; border: 1px solid #ddd; }}
            .error {{ color: red; }}
            .warning {{ color: orange; }}
            .success {{ color: green; }}
            table {{ border-collapse: collapse; width: 100%; }}
            th, td {{ border: 1px solid #ddd; padding: 8px; }}
            th {{ background-color: #f2f2f2; }}
        </style>
    </head>
    <body>
        <h1>Pipeline Report - {pair_id}</h1>
        
        <div class="section">
            <h2>Alignment Metrics</h2>
            {generate_alignment_section(metrics.get('alignment', {}))}
        </div>
        
        <div class="section">
            <h2>Variant Metrics</h2>
            {generate_variant_section(metrics.get('variants', {}))}
        </div>
    </body>
    </html>
    """
    return html_template

def generate_alignment_section(step_metrics):
    """Generate HTML for alignment metrics"""
    if not isinstance(step_metrics, dict):
        return "<p class='error'>No alignment metrics available</p>"
    
    html = "<table>"
    html += "<tr><th>Metric</th><th>Value</th></tr>"
    for k, v in step_metrics.items():
        html += f"<tr><td>{k}</td><td>{v}</td></tr>"
    html += "</table>"
    return html

def generate_variant_section(variant_metrics):
    """Generate HTML for variant metrics"""
    if not isinstance(variant_metrics, dict):
        return "<p class='error'>No variant metrics available</p>"
    
    html = "<table>"
    html += "<tr><th>Metric</th><th>Value</th></tr>"
    for k, v in variant_metrics.items():
        html += f"<tr><td>{k}</td><td>{v}</td></tr>"
    html += "</table>"
    return html


def convert_paths_to_str(cmd_list: List[Any]) -> List[str]:
    """
    Convert all Path objects in a list to strings.

    Parameters:
        cmd_list (List[Any]): List containing strings or Path objects.

    Returns:
        List[str]: List with all elements as strings.
    """
    return [str(item) if isinstance(item, Path) else item for item in cmd_list]
