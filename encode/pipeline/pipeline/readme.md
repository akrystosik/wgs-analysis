# ENCODE WGS Analysis Pipeline

A comprehensive pipeline for processing whole genome sequencing (WGS) data from ENCODE FASTQ files through alignment, QC, and variant calling, designed to run in a high-performance computing environment.

## Directory Structure

```
encode_analysis/
├── pipeline/
│   ├── core/                # Core pipeline functionality
│   │   ├── pipeline.py      # Main pipeline class
│   │   ├── monitoring.py    # Pipeline monitoring
│   │   ├── local_runner.py  # Local DeepVariant execution
│   │   └── exceptions.py    # Custom exceptions
│   ├── utils/               # Utility functions
│   │   └── helpers.py
│   ├── io/                  # Input/output operations
│   │   ├── fastq.py        # FASTQ processing
│   │   ├── downloader.py   # Data download utilities
│   │   └── verify.py       # File verification
│   ├── qc/                  # Quality control modules
│   │   ├── bam.py          # BAM QC
│   │   ├── fastq.py        # FASTQ QC
│   │   └── vcf.py          # VCF QC
│   └── cli.py              # Command-line interface
├── data/                    # Data directory
│   ├── fastq/              # FASTQ files
│   ├── bam/                # BAM outputs
│   ├── qc/                 # QC reports
│   ├── variants/           # Variant calls
│   └── reference/          # Reference genome
└── config.yaml             # Pipeline configuration
```

## Prerequisites

### System Requirements

- **CPU:** 128 cores recommended
- **RAM:** 2TB recommended
- **Storage:** 1TB+ free space
- **GPU:** NVIDIA H100 or similar for DeepVariant

### Software Requirements

- **Python:** 3.8+
- **BWA:** v0.7.17+
- **Samtools:** v1.15+
- **FastQC**
- **DeepVariant:** v1.6.0

## Configuration

The pipeline is configured through a YAML file (`config.yaml`) with settings for:

- **Data Directories:** Paths for input and output data
- **Tool Parameters:** Configuration for BWA, Samtools, FastQC, DeepVariant
- **QC Thresholds:** Metrics thresholds for quality control steps
- **Resource Allocation:** CPU and memory specifications for different pipeline steps

## Usage

### Basic Usage

```bash
python pipeline/cli.py \
  --pair-id <sample_id> \
  --input-fastq <read1.fastq.gz> <read2.fastq.gz> \
  --config config.yaml
```

Or with a pre-aligned BAM file:
```bash
python pipeline/cli.py \
  --pair-id <sample_id> \
  --input-bam <input.bam> \
  --config config.yaml
```

### Command Line Arguments

- `--pair-id`: **(Required)** Sample identifier
- `--input-fastq`: Paths to paired FASTQ files
- `--input-bam`: Path to input BAM file
- `--config`: Path to configuration YAML file
- `--force`: Force rerun of all steps
- `--force-qc`: Continue pipeline execution even if QC checks fail

## Pipeline Steps

1. ### FASTQ Processing
   - FastQC quality assessment
   - Read quality validation
   - Output: FastQC reports in `data/qc/fastqc/`

2. ### Alignment
   - BWA-MEM alignment
   - SAMtools sorting and indexing
   - Output: Sorted BAM files in `data/bam/`

3. ### BAM QC
   - Coverage analysis
   - Mapping statistics
   - Insert size metrics
   - Output: QC reports in `data/qc/alignment/`

4. ### Variant Calling
   - Local DeepVariant execution
   - VCF filtering
   - Output: Filtered VCF files in `data/variants/`

## QC Metrics and Thresholds

### Alignment QC
- Minimum mapping rate: 95%
- Properly paired reads: ≥90%
- Insert size range: 200-600bp

### Coverage QC
- Minimum mean coverage: 20x
- Target coverage: 30x
- Coverage uniformity: ≥90%

### Variant QC
- Minimum depth: 10x
- Minimum quality: 20
- Ts/Tv ratio: 2.0-2.1

## Output Files

### QC Reports
- `data/qc/fastqc/<sample>_fastqc.html`
- `data/qc/alignment/<sample>_metrics.json`
- `data/qc/alignment/<sample>_report.html`

### Pipeline Outputs
- `data/bam/<sample>.bam`
- `data/bam/<sample>.bam.bai`
- `data/variants/<sample>.vcf.gz`
- `data/variants/<sample>.filtered.vcf.gz`
- `data/variants/<sample>.filtered.vcf.gz.tbi`

## Error Handling

- QC failures halt execution unless `--force-qc` is used
- Resource monitoring tracks CPU and memory usage
- Comprehensive logging in the `logs/` directory

## Monitoring

Progress monitoring available through:
- Real-time logging
- Resource usage tracking
- Final QC reports