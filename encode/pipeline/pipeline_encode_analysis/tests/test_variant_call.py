#!/usr/bin/env python3
#encode_analysis/pipeline_v1_fasta_bam/tests/test_variant_call.py

import logging
import os
from pathlib import Path
from pipeline.core.local_runner import LocalDeepVariantRunner

# Setup logging
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def test_variant_call():
    """Test DeepVariant calling with minimal setup."""
    # Test parameters
    bam_file = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/encode_analysis/data/bam/ENCSR456SNK_marked.bam"
    output_dir = "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/encode_analysis/data/variants"
    reference = "encode_analysis/data/reference/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
    
    # Basic config
    config = {
        "base_dir": "/mnt/czi-sci-ai/intrinsic-variation-gene-ex/encode_analysis/data",
        "data_dirs": {
            "variants": "variants",
            "reference": "reference"
        },
        "tools": {
            "deepvariant": {
                "num_shards": 128,
                "model_type": "WGS"
            }
        }
    }
    
    try:
        # Verify files exist
        for path in [bam_file, reference]:
            if not os.path.exists(path):
                raise FileNotFoundError(f"Required file not found: {path}")
        
        logger.info("Starting variant calling test")
        logger.debug(f"Using BAM file: {bam_file}")
        logger.debug(f"Using reference: {reference}")
        
        runner = LocalDeepVariantRunner(config)
        success = runner.run_deepvariant(
            pair_id="ENCSR456SNK",
            bam_file=bam_file,
            reference_path=reference,
            output_vcf=f"{output_dir}/ENCSR456SNK.vcf.gz",
            output_gvcf=f"{output_dir}/ENCSR456SNK.g.vcf.gz"
        )
        
        logger.info(f"Variant calling {'succeeded' if success else 'failed'}")
        
    except Exception as e:
        logger.error(f"Test failed: {str(e)}")
        logger.error(f"Error details: {type(e).__name__}")
        raise

if __name__ == "__main__":
    test_variant_call()
