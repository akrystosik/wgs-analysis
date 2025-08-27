#!/usr/bin/env python3
import os
import json
import logging
import time
from pathlib import Path
from datetime import datetime
from typing import Dict, Any, Optional
from utils.naming import get_standardized_names  # Add this import

logger = logging.getLogger(__name__)

class WGSMonitor:
    """Monitor and record pipeline progress and resource usage."""
    def __init__(self, config: Dict):
        self.config = config
        self.base_dir = Path(config['base_dir'])
        self.monitor_dir = self.base_dir / 'monitoring'
        self.monitor_dir.mkdir(parents=True, exist_ok=True)
        
    def start_pipeline(self, pair_id: str):
        """Record pipeline start."""
        self.logger = logging.getLogger(__name__)
        self.logger.info(f"Started pipeline monitoring for {pair_id}")
        
        monitor_file = self.monitor_dir / f"{pair_id}_pipeline_status.json"
        status = {
            "pair_id": pair_id,
            "pipeline_status": {
                "start_time": str(datetime.now()),
                "status": "running",
                "steps": {},
                "resources": [],
                "errors": []
            },
            "step_status": {}
        }
        
        self._update_monitor_file(monitor_file, status)
    
    def end_pipeline(self, pair_id: str, status: str):
        """Record pipeline completion."""
        self.logger.info(f"Pipeline {pair_id} completed with status: {status}")
        
        monitor_file = self.monitor_dir / f"{pair_id}_pipeline_status.json"
        if monitor_file.exists():
            status_data = self._read_monitor_file(monitor_file)
            status_data["pipeline_status"]["status"] = status
            status_data["pipeline_status"]["end_time"] = str(datetime.now())
            
            # Calculate duration
            start_time = datetime.fromisoformat(status_data["pipeline_status"]["start_time"])
            end_time = datetime.fromisoformat(status_data["pipeline_status"]["end_time"])
            duration = (end_time - start_time).total_seconds()
            
            self.logger.info(f"Duration: {duration} seconds")
            self._update_monitor_file(monitor_file, status_data)
        
        # Generate summary
        self._generate_summary(pair_id)
    
    def start_step(self, step_name: str, pair_id: str):
        """Record the start of a pipeline step."""
        monitor_file = self.monitor_dir / f"{pair_id}_pipeline_status.json"
        if monitor_file.exists():
            status_data = self._read_monitor_file(monitor_file)
            status_data["step_status"][step_name] = {
                "start_time": str(datetime.now()),
                "status": "running"
            }
            self._update_monitor_file(monitor_file, status_data)
    
    def end_step(self, step_name: str, pair_id: str, status: str):
        """Record the completion of a pipeline step."""
        monitor_file = self.monitor_dir / f"{pair_id}_pipeline_status.json"
        if monitor_file.exists():
            status_data = self._read_monitor_file(monitor_file)
            if step_name in status_data["step_status"]:
                status_data["step_status"][step_name]["status"] = status
                status_data["step_status"][step_name]["end_time"] = str(datetime.now())
                self._update_monitor_file(monitor_file, status_data)
    
    def _read_monitor_file(self, file_path: Path) -> Dict:
        """Read status from monitor file."""
        try:
            if file_path.exists():
                with open(file_path) as f:
                    return json.load(f)
            return {}
        except Exception as e:
            self.logger.error(f"Failed to read monitor file: {str(e)}")
            return {}
    
    def _update_monitor_file(self, file_path: Path, status: Dict):
        """Update status in monitor file."""
        try:
            with open(file_path, 'w') as f:
                json.dump(status, f, indent=2)
        except Exception as e:
            self.logger.error(f"Failed to update monitor file: {str(e)}")
 
    
    def _generate_summary(self, pair_id: str):
        """
        Generate a comprehensive summary of the pipeline run, including:
        - Pipeline status and timing
        - Standardized file locations
        - Step-by-step history
        """
        try:
            monitor_file = self.monitor_dir / f"{pair_id}_pipeline_status.json"
            summary_file = self.monitor_dir / f"{pair_id}_summary.json"
            
            if monitor_file.exists():
                status_data = self._read_monitor_file(monitor_file)
                names = get_standardized_names(pair_id)
                
                # Create enhanced summary combining both approaches
                summary = {
                    "pair_id": pair_id,
                    "pipeline_status": status_data["pipeline_status"]["status"],
                    "start_time": status_data["pipeline_status"]["start_time"],
                    "steps": status_data["step_status"],
                    
                    # Add standardized output file locations
                    "output_files": {
                        "bam": str(self.base_dir / 'data/bam' / names['sample_dir'] / names['marked_bam']),
                        "vcf": str(self.base_dir / 'data/variants/deepvariant' / names['vcf_output'])
                    }
                }
                
                # Add timing information if available
                if "end_time" in status_data["pipeline_status"]:
                    summary["end_time"] = status_data["pipeline_status"]["end_time"]
                    start_time = datetime.fromisoformat(status_data["pipeline_status"]["start_time"])
                    end_time = datetime.fromisoformat(status_data["pipeline_status"]["end_time"])
                    summary["duration_seconds"] = (end_time - start_time).total_seconds()
                
                # Write enhanced summary
                with open(summary_file, 'w') as f:
                    json.dump(summary, f, indent=2)
                    
                self.logger.info(f"Pipeline summary generated: {summary_file}")
                
        except Exception as e:
            self.logger.error(f"Failed to generate summary: {str(e)}")