# WGS Analysis Pipeline

Multi-omics whole genome sequencing variant analysis pipeline with paired RNA-seq datasets for integrated genomics research.

## Overview

This repository contains a comprehensive WGS analysis pipeline focused on datasets with paired transcriptomic data, enabling multi-omics integration and analysis. All datasets include both genomic variants and corresponding RNA-seq expression data.

## Datasets

### Multi-omics Datasets (WGS + RNA-seq)
- **ENCODE Paired Cell Lines** - 7 cell lines with comprehensive multi-omics data
- **ADNI** - Alzheimer's Disease Neuroimaging Initiative cohort  
- **GTEx** - Genotype-Tissue Expression project samples
- **MAGE** - Multi-ancestry Analysis of Gene Expression study

## Pipeline Components

### Variant Calling
- **DeepVariant** - Germline variant calling
- **DeepSomatic** - Somatic variant detection  
- Quality control and filtering workflows

### Population Genetics
- **Ancestry Inference** - PCA-based population structure analysis
- **ADMIXTURE** - Ancestry composition estimation
- Multi-population comparative analysis

### Data Integration
- Cross-platform variant harmonization
- Multi-dataset quality assurance
- Standardized metadata with Croissant JSON-LD

## Metadata

Scientific metadata following **Croissant JSON-LD** specification:
- [ADNI WGS Variants](metadata/adni_wgs_variants.jsonld)
- [ENCODE Paired WGS](metadata/encode_paired_wgs_variants.jsonld) 
- [GTEx WGS Variants](metadata/gtex_wgs_variants.jsonld)
- [MAGE WGS Variants](metadata/mage_wgs_variants.jsonld)

All metadata files validated against CZI Science Data Registry standards.

## Quick Start

```bash
# Clone repository
git clone https://github.com/akrystosik/wgs-analysis.git
cd wgs-analysis

# Install dependencies
pip install -r requirements.txt

# Run variant calling pipeline
bash scripts/run_variant_calling.sh

# Perform ancestry analysis
python pca/scripts/ancestry_analysis.py
```

## Structure

```
wgs-analysis/
├── metadata/           # Croissant JSON-LD metadata files
├── scripts/           # Analysis and utility scripts  
├── pca/              # Ancestry inference pipeline
├── config/           # Pipeline configuration
├── docs/             # Documentation
└── results/          # Analysis outputs
```

## Documentation

- [Installation Guide](docs/INSTALLATION.md)
- [Methods](docs/ENCODE_WGS_METHODS_FINAL.md)
- [Results Summary](docs/RESULTS_SUMMARY.md)
- [Troubleshooting](docs/TROUBLESHOOTING.md)

## Citation

If you use this pipeline or datasets, please cite:

```bibtex
@software{wgs_analysis_pipeline,
  title = {Multi-omics WGS Analysis Pipeline},
  author = {Krystosik, Amy},
  year = {2025},
  url = {https://github.com/akrystosik/wgs-analysis}
}
```
---
**Note**: This pipeline focuses specifically on WGS datasets with paired RNA-seq data for comprehensive multi-omics analysis.
