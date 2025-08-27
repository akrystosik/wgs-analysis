# ENCODE Analysis Pipeline

This repository contains a Whole Genome Sequencing (WGS) pipeline for ENCODE analysis. The pipeline handles processing from FASTQ to variant calls using DeepVariant and DeepSomatic.

## Repository Structure

- `encode_analysis/`: Main project directory
  - `pipeline/`: Core pipeline implementation
  - `scripts/`: Utility scripts
  - `configs/`: Configuration files
  - `logs/`: Log output directory

## Usage

Basic usage:

```bash
python3 -m encode_analysis.pipeline.pipeline --pair-id SAMPLE_ID --input-bam INPUT_BAM --config CONFIG_FILE --environment deepvariant
```

For more details, see the documentation.
