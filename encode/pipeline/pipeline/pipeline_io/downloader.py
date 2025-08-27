#pipeline/io/downloader.py
# #!/usr/bin/env python3

import os
import logging
import hashlib
import requests
import shutil
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm
import subprocess
from typing import List, Dict, Optional, Union, Tuple

logger = logging.getLogger(__name__)



class FastqDownloader:
    """Handles downloading and verification of FASTQ and reference files."""


    
# In FastqDownloader.__init__:
    def __init__(self, config: dict):
        """Initialize with configuration."""
        self.config = config
        self.base_dir = Path(config['base_dir'])
        self.fastq_dir = self.base_dir / config['data_dirs']['fastq']
        self.ref_dir = self.base_dir / config['data_dirs']['reference']
        self.temp_dir = self.base_dir / 'temp'
        self.max_workers = config.get('download_workers', 4)
        self.logger = logging.getLogger(__name__)  # Add this line
        
        # Create necessary directories
        self.fastq_dir.mkdir(parents=True, exist_ok=True)
        self.ref_dir.mkdir(parents=True, exist_ok=True)
        self.temp_dir.mkdir(parents=True, exist_ok=True)



    def _is_url(self, path: str) -> bool:
        """Check if a path is a URL."""
        return path.startswith(('http://', 'https://', 'ftp://', 's3://'))

    def _get_file_size(self, path: Union[str, Path]) -> Optional[int]:
        """Get file size for either URLs or local files."""
        if isinstance(path, str) and self._is_url(path):
            try:
                response = requests.head(path, allow_redirects=True)
                return int(response.headers.get('content-length', 0))
            except Exception as e:
                logger.warning(f"Could not get file size for URL {path}: {str(e)}")
                return None
        else:
            try:
                return Path(path).stat().st_size
            except Exception as e:
                logger.warning(f"Could not get file size for file {path}: {str(e)}")
                return None

    def _handle_fastq(self, input_path: Union[str, Path], pair_number: int) -> Tuple[Path, bool]:
        """
        Handle a single FASTQ file, whether local or remote.
        Returns (local_path, needs_download)
        """
        if isinstance(input_path, str) and self._is_url(input_path):
            local_path = self.fastq_dir / Path(input_path).name
            return local_path, True
        else:
            input_path = Path(input_path)
            if not input_path.exists():
                raise FileNotFoundError(f"FASTQ file not found: {input_path}")
            
            # If file is already in fastq_dir, use it directly
            if input_path.parent == self.fastq_dir:
                return input_path, False
            
            # Otherwise, copy or symlink to fastq_dir
            local_path = self.fastq_dir / input_path.name
            if not local_path.exists():
                if self.config.get('symlink_fastq', False):
                    local_path.symlink_to(input_path)
                    logger.info(f"Created symlink for FASTQ {pair_number}: {local_path}")
                else:
                    shutil.copy2(input_path, local_path)
                    logger.info(f"Copied FASTQ {pair_number} to: {local_path}")
            return local_path, False


    def _download_file(self, url: str, output_path: Path, desc: str = None) -> bool:
        """Download a file with progress bar and verification."""
        try:
            headers = {'accept': 'application/json'} if 'encodeproject.org' in url else {}
            remote_size = self._get_file_size(url)
            if output_path.exists():
                local_size = output_path.stat().st_size
                if remote_size and local_size == remote_size:
                    logger.info(f"File {output_path} already exists and is complete")
                    return True
                    
            # Download file with wget for better resume capability
            cmd = [
                'wget', '-c', '--progress=dot:giga',
                '-O', str(output_path),
                url
            ]
            
            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True
            )
            
            # Monitor download progress
            with tqdm(total=remote_size, unit='B', unit_scale=True, desc=desc) as pbar:
                while True:
                    output = process.stderr.readline()
                    if output == '' and process.poll() is not None:
                        break
                    if output:
                        # Update progress bar based on wget output
                        if '%' in output:
                            try:
                                current = int(float(output.split('%')[0]) * remote_size / 100)
                                pbar.n = current
                                pbar.refresh()
                            except:
                                pass
                
            if process.returncode == 0:
                logger.info(f"Successfully downloaded {url}")
                return True
            else:
                logger.error(f"Failed to download {url}")
                return False
                
        except Exception as e:
            logger.error(f"Error downloading {url}: {str(e)}")
            return False


    def process_fastqs(self, fastq_paths: List[Union[str, Path]]) -> Dict[str, Path]:
        """Process FASTQ files and ensure they're available."""
        if len(fastq_paths) != 2:
            raise ValueError(f"Expected exactly 2 FASTQ files, got {len(fastq_paths)}")
        
        self.logger.info(f"Starting FASTQ processing for: {fastq_paths}")
        results = {}
        download_tasks = []
        
        # First, handle all files and identify which need downloading
        for i, path in enumerate(fastq_paths, 1):
            try:
                self.logger.info(f"Processing FASTQ {i}: {path}")
                local_path, needs_download = self._handle_fastq(path, i)
                if needs_download:
                    self.logger.info(f"File {path} needs download, adding to tasks")
                    download_tasks.append((path, local_path))
                results[str(path)] = local_path
                self.logger.info(f"FASTQ {i} processed: local_path={local_path}, needs_download={needs_download}")
            except Exception as e:
                self.logger.error(f"Failed to process FASTQ {i}: {path}")
                self.logger.error(f"Error details: {str(e)}")
                raise
        
        # Download any remote files in parallel
        if download_tasks:
            logger.info(f"Downloading {len(download_tasks)} FASTQ files...")
            with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
                futures = {}
                for url, output_path in download_tasks:
                    desc = f"Downloading {output_path.name}"
                    future = executor.submit(self._download_file, url, output_path, desc)
                    futures[future] = url
                
                # Wait for all downloads and check results
                for future in futures:
                    url = futures[future]
                    if not future.result():
                        raise Exception(f"Failed to download {url}")
        
        return results

    def download_reference(self, ref_url: str) -> Optional[Path]:
        """Download and prepare reference genome."""
        logger.info("Downloading reference genome...")
        
        # Create encode subdirectory
        encode_ref_dir = self.ref_dir / 'encode'
        encode_ref_dir.mkdir(exist_ok=True)
        
        ref_path = encode_ref_dir / Path(ref_url).name
        if not self._download_file(ref_url, ref_path, "Downloading reference"):
            return None

        # Extract if gzipped
        if ref_path.suffix == '.gz':
            logger.info("Extracting reference genome...")
            extracted_path = ref_path.with_suffix('')
            if not extracted_path.exists():
                subprocess.run(['gunzip', '-k', str(ref_path)], check=True)
            ref_path = extracted_path

        # Create indices
        logger.info("Creating reference indices...")
        subprocess.run(['samtools', 'faidx', str(ref_path)], check=True)
        subprocess.run(['bwa', 'index', str(ref_path)], check=True)

        return ref_path

    def cleanup(self):
        """Clean up temporary files."""
        if self.temp_dir.exists():
            for file in self.temp_dir.iterdir():
                file.unlink()
            self.temp_dir.rmdir()