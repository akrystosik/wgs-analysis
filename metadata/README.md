# WGS Dataset Metadata

This directory contains Croissant JSON-LD metadata files describing whole genome sequencing (WGS) variant datasets used in the multi-omics analysis pipeline.

## Metadata Files

### Open Access Datasets

#### `encode_paired_wgs_variants_croissant.jsonld`
- **Dataset**: ENCODE paired cell lines with RNA-seq
- **Samples**: 7 cell lines (A549, Caki2, GM23248, HepG2, K562, NCI-H460, Panc1)
- **License**: CC0 (fully open access)
- **Variants**: ~10,001 filtered variants
- **Purpose**: Multi-omics WGS-RNA-seq integration analysis

#### `encode_wgs_variants_croissant.jsonld`
- **Dataset**: Complete ENCODE WGS collection
- **License**: CC0 (fully open access)
- **Purpose**: Comprehensive ENCODE genomics analysis

### Protected Access Datasets

#### `adni_wgs_variants_croissant.jsonld`
- **Dataset**: Alzheimer's Disease Neuroimaging Initiative
- **Participants**: ~1,700 individuals
- **Access**: Requires ADNI Data Use Agreement
- **Repository**: https://ida.loni.usc.edu/
- **Purpose**: Alzheimer's disease genetics research

#### `gtex_wgs_variants_croissant.jsonld`
- **Dataset**: Genotype-Tissue Expression Project
- **Participants**: 838 individuals across multiple tissues
- **Access**: Requires dbGaP authorization (phs000424)
- **Repository**: https://www.ncbi.nlm.nih.gov/projects/gap/
- **Purpose**: eQTL and tissue-specific expression analysis

#### `mage_wgs_variants_croissant.jsonld`
- **Dataset**: Multi-ancestry Analysis of Gene Expression
- **Participants**: 731 individuals from 26 global populations
- **Source**: 1000 Genomes Project
- **Access**: Available through 1000 Genomes consortium
- **Purpose**: Population genetics and ancestry inference

## Technical Specifications

### Format Compliance
- **Standard**: Croissant JSON-LD v1.0 specification
- **Validation**: Passed MLCommons Croissant validator
- **Schema**: Compatible with Schema.org genomics extensions

### Genomics Extensions
All metadata files include standardized genomics fields:
- `sc:organism`: Species (Homo sapiens)
- `sc:assayType`: Assay methodology (Whole genome sequencing)  
- `sc:dataType`: Data content (genetic variants)
- `sc:accessLevel`: Access requirements (open/protected)

### Data Access Strategy
- **Open datasets**: Include file metadata and download links
- **Protected datasets**: Reference official repositories only
- **Compliance**: Respects data use agreements and access requirements

## Validation

All metadata files have been validated using the official CZI Science Data Registry validators:
```bash
python validate_vcf_croissant.py
# Result: 4/4 files passed validation ✅
```

## Citation

If using these datasets, please cite the original sources as specified in each metadata file's `citation` field.

## Last Updated

Generated: 2025-08-27
Validator: https://github.com/chanzuckerberg/sci-data-registry/tree/main/validators