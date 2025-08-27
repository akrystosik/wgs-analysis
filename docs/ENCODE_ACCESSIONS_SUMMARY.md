# ENCODE Accessions and Assay Types Summary

## Overview
This document catalogs all ENCODE accessions identified from the project files, organized by data type and assay category.

## RNA-seq Data (ENCFF accessions)

### Cell Line RNA-seq Files
| Cell Line | ENCFF Accession | Assay Type | Experiment Accession | File Type | Assembly | Notes |
|-----------|-----------------|------------|---------------------|-----------|----------|-------|
| A549      | ENCFF244DNJ     | polyA plus RNA-seq | ENCSR000CON | gene quantifications (tsv) | GRCh38 | Lung carcinoma, stranded PE76 |
| K562      | ENCFF171FQU     | total RNA-seq | ENCSR000AEL | gene quantifications (tsv) | GRCh38 | Chronic myelogenous leukemia |
| HepG2     | ENCFF863QWG     | total RNA-seq | ENCSR245ATJ | gene quantifications (tsv) | GRCh38 | Hepatocellular carcinoma |
| GM23248   | ENCFF640FPG     | total RNA-seq | ENCSR797BPP | gene quantifications (tsv) | hg19 | Skin fibroblast, normal |
| Caki2     | ENCFF685WJV     | total RNA-seq | ENCSR584JXD | gene quantifications (tsv) | GRCh38 | Clear cell renal carcinoma |
| NCI-H460  | ENCFF876SRX     | total RNA-seq | ENCSR164OCT | gene quantifications (tsv) | GRCh38 | Large cell lung carcinoma |
| Panc1     | ENCFF710IFD     | total RNA-seq | ENCSR128CYL | gene quantifications (tsv) | GRCh38 | Pancreatic ductal carcinoma |

**Key RNA-seq Assay Details:**
- **EFO:0009922** - RNA sequencing assay (primary ontology ID)
- **OBI:0002571** - polyA plus RNA-seq (A549 specific)
- **OBI:0001271** - RNA-seq (general, used by most cell lines)
- **Processing**: RSEM quantification, GENCODE annotation
- **Selection**: rRNA depleted, >200 nucleotides
- **Platform**: Illumina (various models)
- **Read configuration**: Typically paired-end, strand-specific

## WGS Variant Data (ENCFF accessions)

### Legacy ENCODE VCF Files (hg19→hg38 lifted)
| ENCFF Accession | Variants | Data Type | Notes |
|-----------------|----------|-----------|-------|
| ENCFF752OAX     | 3,775,298 | Comprehensive variant set | Primary large variant dataset |
| ENCFF785JVR     | 3,544 | Structural variants | Specialized variant calls |
| ENCFF574MDJ     | 301 | Small variant set | Limited scope variants |
| ENCFF863MPP     | 296 | Small variant set | Limited scope variants |

**WGS Assay Details:**
- **Assay Type**: Whole genome sequencing (WGS)
- **Assay Term ID**: OBI:0002117
- **Technology**: Illumina sequencing
- **Variant Calling**: Various methods (pre-standardization era)
- **Assembly**: Originally hg19, lifted to hg38/GRCh38

## Experiment-Level Accessions (ENCSR)

### RNA-seq Experiments
| Experiment ID | Cell Line | Assay Description | Lab | Award |
|---------------|-----------|-------------------|-----|-------|
| ENCSR000CON   | A549      | Stranded PE76 polyA+ RNA-seq | ENCODE Consortium | ENCODE4 |
| ENCSR000AEL   | K562      | Total RNA-seq | ENCODE Consortium | ENCODE4 |
| ENCSR245ATJ   | HepG2     | Total RNA-seq | ENCODE Consortium | ENCODE4 |
| ENCSR797BPP   | GM23248   | Total RNA-seq | ENCODE Consortium | ENCODE3 |
| ENCSR584JXD   | Caki2     | Total RNA-seq | ENCODE Consortium | ENCODE4 |
| ENCSR164OCT   | NCI-H460  | Total RNA-seq | ENCODE Consortium | ENCODE4 |
| ENCSR128CYL   | Panc1     | Total RNA-seq | ENCODE Consortium | ENCODE4 |

### WGS/Variant Calling Experiments
| Experiment ID | Sample/Cell Line | Assay Description | Technology | Lab | Notes |
|---------------|------------------|-------------------|------------|-----|-------|
| ENCSR456SNK   | Transverse colon tissue | 10X Chromium WGS | 10X Genomics | Unknown | EN-TEx project |
| ENCSR121TMQ   | [Cell line - referenced in scripts] | WGS | Unknown | Unknown | Referenced in analysis scripts |
| ENCSR053AXS   | K562 | 10X Chromium WGS | 10X Genomics | Unknown | Referenced in project notes |
| ENCSR674PQI   | GM23248 | WGS experiment | Unknown | Unknown | Referenced in project notes |
| ENCSR420NDH   | [Unknown sample] | WGS | Unknown | Unknown | Referenced in project notes |
| **ENCSR521ELB** | **A549** | **WGS** | **Illumina** | **Feng Yue, PSU** | **✅ Confirmed via API lookup** |
| **ENCSR654YMF** | **NCI-H460** | **WGS** | **Illumina** | **Feng Yue, PSU** | **✅ Confirmed via API lookup** |
| **ENCSR696BCR** | **Panc1** | **WGS** | **Illumina** | **Feng Yue, PSU** | **✅ Confirmed via API lookup** |
| **ENCSR319QHO** | **HepG2** | **10X Chromium WGS** | **10X Genomics** | **Alexander Urban, Stanford** | **✅ Confirmed via API lookup** |
| **ENCSR236KYS** | **sknmc** | **WGS** | **Illumina** | **Feng Yue, PSU** | **✅ Confirmed via API lookup** |

## WGS FASTQ File Details (From Pipeline Logs)

### Complete FASTQ-to-Cell Line Mapping
The following ENCFF FASTQ file accessions were extracted from WGS processing pipeline logs:

| Cell Line | FASTQ File 1 | FASTQ File 2 | Pipeline Log Source | ENCSR Status |
|-----------|--------------|--------------|---------------------|--------------|
| **A549** | ENCFF122NPY | ENCFF846WHK | pipeline_A549_WGS_pair1_20250214_195958.log | **ENCSR521ELB** ✅ |
| **GM23248** | ENCFF826GTR | ENCFF851EYG | pipeline_GM23248_WGS_pair1_20250214_223049.log | **ENCSR674PQI** ✅ |
| **NCI-H460** | ENCFF022XPK | ENCFF534EUU | pipeline_NCI-H460_WGS_pair1_20250215_000713.log | **ENCSR654YMF** ✅ |
| **Panc1** | ENCFF477JTA | ENCFF896PZG | pipeline_Panc1_WGS_pair1_20250219_224929.log | **ENCSR696BCR** ✅ |
| **HepG2** | ENCFF320KMG | ENCFF045JFV | pipeline_20241111_172508.log | **ENCSR319QHO** ✅ |
| **sknmc** | ENCFF212SDN | ENCFF675TKC | pipeline_sknmc_WGS_pair1_20250221_142307.log | **ENCSR236KYS** ✅ |
| **K562** | ENCFF066GQD, ENCFF506TKC | ENCFF004THU, ENCFF313MGL | 10X Chromium FASTQ → DeepVariant/DeepSomatic pipeline | **ENCSR053AXS** ✅ |
| **Caki2** | *Not found* | *Not found* | Log files checked, no FASTQ info | *no paired rnaseq* |
| **T47D** | *Not found* | *Not found* | Log files checked, no FASTQ info | *no paired rnaseq* |

### WGS Processing Status Summary
- **Total cell lines processed**: 9 (A549, Caki2, GM23248, HepG2, K562, NCI-H460, Panc1, T47D, sknmc)
- **FASTQ files identified**: 6 cell lines (12 FASTQ files total)
- **Missing FASTQ logs**: 3 cell lines (K562, Caki2, T47D)
- **Confirmed ENCSR mapping**: 7 cell lines ✅ (A549→ENCSR521ELB, GM23248→ENCSR674PQI, HepG2→ENCSR319QHO, K562→ENCSR053AXS, NCI-H460→ENCSR654YMF, Panc1→ENCSR696BCR, sknmc→ENCSR236KYS)
- **ENCSR lookup needed**: 2 cell lines (Caki2, T47D)

### Next Steps for Complete Documentation
To complete the WGS experiment mapping (7/9 cell lines now complete):
1. ~~**Use ENCFF accessions to find ENCSR experiments** via ENCODE portal queries~~ ✅ **COMPLETED** (5/5 remaining lookups successful)
2. **Remaining tasks for Caki2 and T47D**: 
   - Check archived logs for missing FASTQ information 
   - Search alternative log sources or processing scripts
   - Query ENCODE portal directly by cell line biosample if FASTQ files unavailable
3. **Validate experiment-to-cell-line mapping** using ENCODE biosample metadata ✅ **PARTIALLY COMPLETE** (7/9 validated)

## Assay Type Ontology Mapping

### Primary Assay Types Used
| Ontology ID | Term Name | Usage Context |
|-------------|-----------|---------------|
| EFO:0009922 | RNA sequencing assay | Primary RNA-seq ontology for all cell lines |
| OBI:0002117 | whole genome sequencing assay | WGS experiments |
| OBI:0002571 | polyA plus RNA-seq | A549 specific assay |
| OBI:0001271 | RNA-seq | General RNA-seq term |

### Supporting Ontology Terms
| Ontology ID | Term Name | Usage |
|-------------|-----------|-------|
| SO:0000871  | polyadenylated mRNA | Nucleic acid type for polyA+ |
| SO:0000356  | RNA | General RNA nucleic acid type |
| SO:0000352  | DNA | DNA nucleic acid type (WGS) |
| SO:0000252  | rRNA | For depletion specification |

## Data Processing Pipelines

### RNA-seq Processing
- **Pipeline**: ENCODE uniform processing pipeline
- **Software**: RSEM (versions 1.2.19, 1.2.31)
- **Reference**: GENCODE v24/v29, GRCh38/hg38
- **Output**: TPM (Transcripts Per Million) quantification
- **QC**: Standard ENCODE quality metrics

### WGS Processing (Project-specific)
- **Alignment**: BWA-MEM to GRCh38
- **Variant Calling**: DeepVariant (germline), DeepSomatic (somatic)
- **Post-processing**: Picard MarkDuplicates
- **QC**: FastQC, samtools stats, variant statistics

## File Organization Patterns

### RNA-seq Files
```
ENCFF + 6-character ID + .tsv (gene quantifications)
Example: ENCFF244DNJ.tsv
```

### WGS Files
```
ENCFF + 6-character ID + .vcf.gz (variant calls)
Example: ENCFF752OAX.vcf.gz
```

### Experiment IDs
```
ENCSR + 6-character ID (experiment containers)
Example: ENCSR000CON
```

## Summary Statistics

### RNA-seq Data
- **Total ENCFF RNA-seq files**: 7 (across 7 cell lines)
- **Total ENCSR RNA-seq experiments**: 7 (complete mapping)
- **Cell lines with RNA-seq**: 7 (A549, K562, HepG2, GM23248, Caki2, NCI-H460, Panc1)

### WGS Data  
- **Total ENCFF WGS files**: 4 (legacy variant sets) + 12 (FASTQ files from processing)
- **Total ENCSR WGS experiments**: 5+ (2 confirmed, 7+ need lookup)
- **WGS FASTQ files identified**: 12 (6 cell lines with paired-end data)
- **Cell lines with WGS processing**: 9 (A549, Caki2, GM23248, HepG2, K562, NCI-H460, Panc1, T47D, sknmc)
- **WGS variant calling completed**: All 9 cell lines (DeepVariant + DeepSomatic)

### Overall Project Coverage
- **Primary assay types**: 2 (RNA-seq, WGS)
- **Cell lines with both RNA-seq and WGS**: 7 (A549, K562, HepG2, GM23248, Caki2, NCI-H460, Panc1)
- **Cell lines with WGS only**: 2 (T47D, sknmc)
- **Assemblies used**: GRCh38/hg38 (primary), hg19 (legacy)
- **Total variant calls generated**: ~44M variants across 9 cell lines

## Data Availability and Access

### Public Data
- **RNA-seq files**: All ENCFF files publicly available via ENCODE portal
- **Download URLs**: Direct HTTPS links provided in metadata
- **Cloud Storage**: S3 and Azure URIs available

### Generated Data (Project-specific)
- **Merged VCF**: Combined ENCODE cell line variants
- **Standardized datasets**: H5AD files with harmonized metadata
- **Croissant metadata**: JSON-LD formatted dataset descriptions

This comprehensive inventory provides the foundation for proper attribution, reproducibility, and future integration of ENCODE datasets in the project.