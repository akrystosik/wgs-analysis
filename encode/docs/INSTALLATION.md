# Installation Guide

## Prerequisites

### System Requirements
- Linux/Unix operating system
- Python 3.8+
- 64+ GB RAM (recommended for WGS processing)
- 96+ CPU cores (for optimal performance)
- 1TB+ storage space

### Required Software

#### Core Tools
- **DeepVariant v1.6.0** - Germline variant calling
- **DeepSomatic v1.6.0** - Somatic variant calling  
- **BWA-MEM v0.7.17** - Read alignment
- **FastQC v0.11.9** - Quality control
- **Picard Tools** - BAM processing
- **BCFtools** - VCF manipulation
- **Samtools** - BAM manipulation

#### Python Dependencies
```bash
pip install -r requirements.txt
```

## Installation Steps

### 1. Clone Repository
```bash
git clone https://github.com/akrystosik/encode_analysis.git
cd encode_analysis
```

### 2. Install Dependencies
```bash
# Install Python packages
pip install -r requirements.txt

# Install system tools (example for Ubuntu/Debian)
sudo apt-get update
sudo apt-get install samtools bcftools fastqc bwa
```

### 3. Configure DeepVariant/DeepSomatic
Download and install Google's DeepVariant and DeepSomatic:

```bash
# DeepVariant
curl -O https://github.com/google/deepvariant/releases/download/v1.6.0/deepvariant-1.6.0.tar.gz
tar -xzf deepvariant-1.6.0.tar.gz

# DeepSomatic  
curl -O https://github.com/google/deepvariant/releases/download/v1.6.0/deepsomatic-1.6.0.tar.gz
tar -xzf deepsomatic-1.6.0.tar.gz
```

### 4. Set Up Configuration
```bash
cp config/pipeline_config.yaml.example config/pipeline_config.yaml
# Edit config/pipeline_config.yaml with your paths
```

### 5. Verify Installation
```bash
python -m pipeline.tests.test_variant_call
```

## Docker Installation (Alternative)

```bash
# Build Docker image
docker build -t encode_analysis .

# Run pipeline in container
docker run -v /path/to/data:/data encode_analysis \
    python -m pipeline.pipeline --input-bam /data/sample.bam
```

## Reference Genome Setup

Download GRCh38 reference:
```bash
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
bwa index GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
```

## Troubleshooting

See [TROUBLESHOOTING.md](TROUBLESHOOTING.md) for common issues and solutions.