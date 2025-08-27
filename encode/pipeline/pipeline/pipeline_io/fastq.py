# Improved pipeline_io/fastq.py processing functions

import os
import subprocess
import logging
import time
from pathlib import Path
import math
import psutil
import shutil  
import traceback  

logger = logging.getLogger(__name__)

def process_fastqs(config, fastq_files, pair_id, force=False):
    """Process FASTQ files with enhanced resource management and progress monitoring."""
    logger.info(f"Starting enhanced FASTQ processing for {len(fastq_files)} files")
    
    # Calculate optimal resources based on file sizes and system capabilities
    fastq_sizes = []
    total_size_gb = 0
    for fastq in fastq_files:
        if os.path.exists(fastq):
            size_gb = os.path.getsize(fastq) / (1024**3)
            fastq_sizes.append(size_gb)
            total_size_gb += size_gb
            logger.info(f"FASTQ file {fastq} size: {size_gb:.2f} GB")
        else:
            logger.error(f"FASTQ file not found: {fastq}")
            raise FileNotFoundError(f"FASTQ file not found: {fastq}")
    
    # Determine optimal thread count based on file size and available CPUs
    available_cpus = psutil.cpu_count(logical=True)
    available_memory_gb = psutil.virtual_memory().available / (1024**3)
    
    # Calculate optimal thread allocation
    bwa_threads = min(max(16, int(available_cpus * 0.7)), 64)
    sort_threads = min(max(8, int(available_cpus * 0.3)), 32)
    
    # Determine memory per thread for sorting based on file size
    # Formula: larger files need more memory per thread
    estimated_bam_size_gb = total_size_gb * 2  # BAM is typically 2x larger than FASTQ
    optimal_memory_per_thread_gb = min(
        max(2, int(estimated_bam_size_gb / 50)),  # At least 2GB, scales with file size
        int(available_memory_gb / sort_threads * 0.8)  # Don't use more than 80% of available memory
    )
    
    # Log resource allocation plan
    logger.info(f"Resource allocation plan:")
    logger.info(f"  - Total FASTQ size: {total_size_gb:.2f} GB")
    logger.info(f"  - Available resources: {available_cpus} CPUs, {available_memory_gb:.2f} GB memory")
    logger.info(f"  - BWA threads: {bwa_threads}")
    logger.info(f"  - Sort threads: {sort_threads}")
    logger.info(f"  - Memory per sort thread: {optimal_memory_per_thread_gb}G")
    
    # Get sample name from pair_id (simplified here - use your actual logic)
    sample_name = pair_id.split('_')[0]
    
    # Create sample directory
    sample_dir = Path(config['base_dir']) / 'data/bam' / sample_name
    sample_dir.mkdir(parents=True, exist_ok=True)
    
    # Create timestamps for temp directories to avoid conflicts
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    temp_dir = sample_dir / f"temp_{timestamp}"
    temp_dir.mkdir(parents=True, exist_ok=True)
    logger.info(f"Created temporary directory: {temp_dir}")
    
    # Output BAM file path
    sorted_bam = sample_dir / f"{pair_id}.sorted.bam"
    marked_bam = sample_dir / f"{pair_id}.marked_duplicates.bam"
    
    # Skip processing if marked BAM already exists and force is False
    if os.path.exists(marked_bam) and not force:
        logger.info(f"Marked duplicates BAM already exists: {marked_bam}")
        return str(marked_bam)
    
    logger.info(f"Creating of marked duplicates BAM")
    
    # Step 1: Align with BWA and pipe to sorting with progress monitoring
    reference_path = Path(config['base_dir']) / config['data_dirs']['reference'] / config['reference']['fasta']
    
    # Add read group for this sample
    read_group = f"@RG\\tID:{pair_id}\\tSM:{pair_id}\\tLB:lib1\\tPL:ILLUMINA"
    
    logger.info(f"Starting alignment pipeline...")
    
    # Create a checkpoint file
    alignment_started_file = temp_dir / "alignment_started"
    with open(alignment_started_file, 'w') as f:
        f.write(f"Started at {time.ctime()}")
    
    # Create BWA command
    bwa_cmd = [
        "bwa", "mem",
        "-t", str(bwa_threads),
        "-R", read_group,
        str(reference_path),
    ] + fastq_files
    
    # Create samtools sort command with optimized parameters
    sort_cmd = [
        "samtools", "sort",
        "-@", str(sort_threads),
        "-T", str(temp_dir / "sort"),
        "-m", f"{optimal_memory_per_thread_gb}G",
        "-o", str(sorted_bam)
    ]
    
    # Log full commands
    logger.info(f"BWA command: {' '.join(bwa_cmd)}")
    logger.info(f"Samtools sort command: {' '.join(sort_cmd)}")
    
    try:
        # Start alignment and sorting with progress monitoring
        start_time = time.time()
        last_update_time = start_time
        
        bwa_process = subprocess.Popen(
            bwa_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=False,
            bufsize=-1  # Use system default buffering
        )
        
        sort_process = subprocess.Popen(
            sort_cmd,
            stdin=bwa_process.stdout,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=-1
        )
        
        # Close the pipe in the parent process
        bwa_process.stdout.close()
        
        # Create threads to monitor stderr output without blocking
        def monitor_stderr(process, name):
            while True:
                line = process.stderr.readline()
                if not line and process.poll() is not None:
                    break
                if line:
                    logger.info(f"{name} output: {line.strip()}")
        
        import threading
        bwa_monitor = threading.Thread(target=monitor_stderr, args=(bwa_process, "BWA"))
        sort_monitor = threading.Thread(target=monitor_stderr, args=(sort_process, "Sort"))
        
        bwa_monitor.daemon = True
        sort_monitor.daemon = True
        
        bwa_monitor.start()
        sort_monitor.start()
        
        # Monitor progress of temporary files
        while bwa_process.poll() is None or sort_process.poll() is None:
            # Sleep briefly to avoid excessive CPU usage
            time.sleep(30)
            
            current_time = time.time()
            # Log progress every 10 minutes
            if current_time - last_update_time > 600:
                elapsed_minutes = (current_time - start_time) / 60
                
                # Calculate progress based on temp files
                temp_files = list(temp_dir.glob("sort.*.bam"))
                total_temp_size = sum(f.stat().st_size for f in temp_files)
                
                # Calculate disk usage
                disk_usage = shutil.disk_usage(str(sample_dir))
                disk_percent = disk_usage.used / disk_usage.total * 100
                
                logger.info(f"Progress update after {elapsed_minutes:.1f} minutes:")
                logger.info(f"  - Temporary files: {len(temp_files)} files, {total_temp_size / 1024**3:.2f} GB")
                logger.info(f"  - Disk usage: {disk_percent:.1f}% used")
                logger.info(f"  - BWA running: {bwa_process.poll() is None}")
                logger.info(f"  - Sort running: {sort_process.poll() is None}")
                
                # Check if processes are making progress
                if len(temp_files) > 0:
                    with open(temp_dir / "last_progress", 'w') as f:
                        f.write(f"{time.time()},{len(temp_files)},{total_temp_size}")
                
                last_update_time = current_time
        
        # Check for errors
        if bwa_process.returncode != 0:
            bwa_error = bwa_process.stderr.read()
            logger.error(f"BWA failed with return code {bwa_process.returncode}: {bwa_error}")
            raise RuntimeError(f"BWA alignment failed: {bwa_error}")
        
        if sort_process.returncode != 0:
            sort_error = sort_process.stderr.read()
            logger.error(f"Samtools sort failed with return code {sort_process.returncode}: {sort_error}")
            raise RuntimeError(f"Samtools sort failed: {sort_error}")
        
        # Verify the sorted BAM exists and has a reasonable size
        if not os.path.exists(sorted_bam):
            logger.error(f"Sorted BAM file not found after alignment: {sorted_bam}")
            raise FileNotFoundError(f"Sorted BAM file not found: {sorted_bam}")
        
        sorted_bam_size = os.path.getsize(sorted_bam) / (1024**3)
        logger.info(f"Alignment and sorting completed. Sorted BAM size: {sorted_bam_size:.2f} GB")
        
        if sorted_bam_size < 0.1:
            logger.warning(f"Sorted BAM file is suspiciously small: {sorted_bam_size:.2f} GB")
        
        # Step 2: Mark duplicates
        logger.info("Starting duplicate marking...")
        
        # Create MarkDuplicates command with metrics
        metrics_file = sample_dir / f"{pair_id}.duplicate_metrics.txt"
        
        mark_cmd = [
            "picard", "MarkDuplicates",
            f"INPUT={sorted_bam}",
            f"OUTPUT={marked_bam}",
            f"METRICS_FILE={metrics_file}",
            "VALIDATION_STRINGENCY=LENIENT",
            "ASSUME_SORTED=true",
            f"TMP_DIR={temp_dir}",
            f"MAX_RECORDS_IN_RAM={10000000}"  # Adjust based on available memory
        ]
        
        logger.info(f"Mark duplicates command: {' '.join(mark_cmd)}")
        
        mark_process = subprocess.Popen(
            mark_cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True
        )
        
        # Monitor MarkDuplicates output
        for line in mark_process.stdout:
            logger.info(f"MarkDuplicates: {line.strip()}")
        
        mark_process.wait()
        
        if mark_process.returncode != 0:
            logger.error(f"MarkDuplicates failed with return code {mark_process.returncode}")
            raise RuntimeError("MarkDuplicates failed")
        
        # Verify marked BAM file
        if not os.path.exists(marked_bam):
            logger.error(f"Marked BAM file not found after processing: {marked_bam}")
            raise FileNotFoundError(f"Marked BAM file not found: {marked_bam}")
        
        marked_bam_size = os.path.getsize(marked_bam) / (1024**3)
        logger.info(f"Duplicate marking completed. Marked BAM size: {marked_bam_size:.2f} GB")
        
        # Step 3: Index the marked BAM
        logger.info("Indexing marked BAM file...")
        index_cmd = ["samtools", "index", str(marked_bam)]
        subprocess.run(index_cmd, check=True)
        
        # Calculate processing time
        total_time = (time.time() - start_time) / 60
        logger.info(f"FASTQ processing completed in {total_time:.1f} minutes")
        
        # Optionally clean up temporary files if needed
        if config.get('cleanup_temp', True):
            logger.info(f"Cleaning up temporary directory: {temp_dir}")
            try:
                shutil.rmtree(temp_dir)
            except Exception as e:
                logger.warning(f"Failed to clean up temporary directory: {e}")
        
        return str(marked_bam)
        
    except Exception as e:
        logger.error(f"FASTQ processing failed: {str(e)}")
        logger.error(f"Traceback: {traceback.format_exc()}")
        
        # Write failure information
        with open(temp_dir / "processing_failed", 'w') as f:
            f.write(f"Failed at {time.ctime()}: {str(e)}")
        
        raise