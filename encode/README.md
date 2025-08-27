# ENCODE Variant Calling Pipeline

## Overview

This directory contains the complete ENCODE variant calling pipeline used to process whole genome sequencing data from ENCODE cell lines and GTEx samples. The pipeline produces high-quality variant calls that serve as input for the ancestry inference analysis in `../pca/`.

## Pipeline Components

### Core Pipeline (`pipeline/`)
- **DeepVariant**: Germline variant calling
- **DeepSomatic**: Somatic variant calling for cell lines
- **Quality Control**: Comprehensive QC and validation
- **Monitoring**: Pipeline execution tracking

### Analysis Scripts (`scripts/`)
- Variant concordance analysis
- Statistical summarization
- Comparison tools between callers
- GTEx batch processing
- ENCODE data integration

### Utilities (`utils/`)
- VCF processing and filtering
- Data extraction tools
- Format conversion utilities
- Validation helpers

### Data Processing (`data/`)
- **BAM files**: Aligned sequencing data
- **VCF files**: Variant call results
- **Reference**: GRCh38 reference genome
- **QC reports**: Quality control outputs

## Key Features

- **Multi-caller approach**: DeepVariant + DeepSomatic
- **Comprehensive QC**: Multiple validation layers
- **Scalable processing**: Handles thousands of samples
- **Standards compliance**: ENCODE/GTEx compatible outputs
- **Integration ready**: Outputs feed into ancestry pipeline

## Usage

The pipeline processes raw sequencing data through variant calling to produce filtered VCF files. These VCFs are then used by the ancestry inference pipeline in `../pca/`.

## Outputs

- High-quality variant calls for 2,331+ samples
- Comprehensive QC metrics and reports
- Standardized metadata for downstream analysis
- Cell line and GTEx population data

## Connection to Ancestry Pipeline

This variant calling pipeline produces the genomic variant data that is processed by the ancestry inference pipeline (`../pca/`) to perform population structure analysis and ancestry assignment.

---

**Status**: Production pipeline  
**Samples**: 2,331+ ENCODE cell lines and GTEx individuals  
**Variants**: 13.9M high-quality variants genome-wide