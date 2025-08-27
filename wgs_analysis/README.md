# WGS Variant Analysis Pipeline

This directory contains scripts for analyzing variants from whole genome sequencing (WGS) data, particularly focused on comparing GTEx samples with cell line data.

## Key Scripts

- `bcftools_bychromosome.py`: Processes VCF files chromosome by chromosome
- `bcftools_stats.py`: Runs bcftools stats on multiple VCF files with parallel processing
- `count_variants.sh`: Counts variants in a single sample
- `direct_count.py`: Performs direct counting of variants from VCF files
- `gtex_extraction.py`: Extracts information from GTEx VCF files
- `process_gtex_batches.sh`: Processes GTEx samples in batches
- `process_gtex_batches_newvcf.sh`: Processes samples from the unphased VCF file
- `process_gtex_samples_simple.sh`: Processes multiple GTEx samples in parallel
- `targeted_analysis.py`: Performs targeted analysis on specific regions

## Results

The analysis revealed key differences between GTEx samples and cell lines:

- Cell lines show an average of 4.33M SNPs and 896K indels
- GTEx samples show an average of 3.32M SNPs and 246K indels
- The most striking difference is in indel counts, with cell lines showing 3.6x more indels

## Analysis Process

1. Located and validated the correct GTEx VCF file
2. Developed a framework for counting variants per sample and chromosome
3. Created parallel processing scripts for efficient analysis
4. Compared results between GTEx samples and cell lines
5. Generated detailed statistics and visualizations

## Usage

To process a single GTEx sample:
```bash
./count_variants.sh GTEX-SAMPLE-ID To process multiple samples in parallel:
bash./process_gtex_samples_simple.sh NUM_SAMPLES NUM_PARALLEL_PROCESSES
To process samples in batches:
bash./process_gtex_batches.sh BATCH_SIZE NUM_PARALLEL_PROCESSES
