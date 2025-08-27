# Methods

## Whole Genome Variant Calling Pipeline for ENCODE Cell Lines

### Sample Data and Preprocessing

We analyzed whole genome sequencing data from nine ENCODE cell lines representing diverse cancer types: A549 (lung adenocarcinoma), Caki2 (renal cell carcinoma), GM23248 (lymphoblastoid cell line), HepG2 (hepatocellular carcinoma), K562 (chronic myelogenous leukemia), NCI-H460 (lung large cell carcinoma), Panc1 (pancreatic adenocarcinoma), T47D (breast ductal carcinoma), and SK-N-MC (neuroepithelioma). Raw sequencing data were obtained from the ENCODE Data Coordination Center and processed using standardized quality control procedures.

Quality assessment was performed using FastQC v0.11.9^1^ with parallel processing (16 threads, 16GB memory per process) to evaluate read quality, adapter content, and sequence duplication levels. Reads were aligned to the GRCh38 reference genome (GCA_000001405.15) using BWA-MEM v0.7.17^2^ configured for optimal performance (64 threads, 128GB memory) with default parameters optimized for whole genome sequencing. Post-alignment processing included duplicate marking using Picard MarkDuplicates v2.27.5^3^ (64GB memory allocation) and quality score recalibration.

### Dual-Caller Variant Detection Strategy

We implemented a comprehensive variant calling approach combining germline and somatic variant detection to capture the full spectrum of genetic variation in cancer cell lines. This dual-caller strategy addresses the unique characteristics of immortalized cell lines, which may harbor both inherited germline variants and acquired somatic mutations.

#### DeepVariant Germline Calling

Germline variants were identified using DeepVariant v1.6.0^4^ with the WGS model trained on whole genome sequencing data. DeepVariant employs a deep neural network architecture to convert genomic regions into image-like representations, enabling accurate variant calling through computer vision techniques. The pipeline was configured with 64 parallel shards distributed across make_examples, call_variants, and postprocess_variants stages to optimize computational efficiency while maintaining accuracy.

Only variants with PASS filter status were retained, relying on DeepVariant's internal quality assessment and filtering mechanisms. No additional depth or quality score thresholds were applied beyond the caller's built-in filtering.

#### DeepSomatic Somatic Variant Detection

Somatic variants were detected using DeepSomatic v1.6.0^5^ operating in tumor-only mode (WGS_TUMOR_ONLY model) with 64 parallel shards matching the DeepVariant configuration. DeepSomatic is specifically designed to identify somatic mutations in cancer samples by modeling the unique mutational patterns and allele frequency distributions characteristic of malignant cells. The algorithm accounts for potential contamination from normal cells and variable tumor purity across different cell line preparations.

DeepSomatic processing included both somatic variant identification and germline variant labeling, allowing for comprehensive characterization of the mutational landscape. Variants were classified as either somatic (PASS filter) or germline-labeled based on algorithmic assessment of allele frequency patterns and genomic context. Only variants with PASS filter status were retained for somatic calls, relying on DeepSomatic's internal quality assessment.

### Variant Integration and Annotation

To create unified variant sets for each cell line, we developed a systematic integration procedure combining non-overlapping calls from both callers. DeepSomatic PASS variants were prioritized for somatic positions, while DeepVariant calls were retained for positions not identified as somatic. This approach ensures comprehensive coverage while avoiding double-counting of variants at identical genomic coordinates.

All variants were annotated with origin type (germline or somatic) using custom VCF INFO fields. Sample names were standardized across datasets, and variants were assigned unique identifiers following the format CHROM_POS_REF_ALT. Integration was performed using BCFtools v1.15^6^ with custom annotation scripts to maintain data provenance and quality metrics.

### Quality Control and Validation

Comprehensive quality control measures were implemented throughout the pipeline to ensure high-confidence variant calls. Coverage analysis, mapping quality metrics, and read pairing statistics were assessed using standard genomics quality control tools. 

Variant quality was evaluated using multiple metrics including transition/transversion (Ts/Tv) ratios, heterozygote/homozygote ratios, and allele frequency distributions. Variant density and distribution patterns were analyzed across chromosomes to identify potential systematic biases or technical artifacts.

### Statistical Analysis and Output Generation

Final variant sets were characterized using comprehensive summary statistics including total variant counts, variant type distributions (SNVs, indels), and functional annotation categories. All analyses were performed using Python 3.9 with pandas^7^, NumPy^8^, and custom bioinformatics libraries. Statistical comparisons were conducted using appropriate non-parametric tests, with significance assessed at α = 0.05. Visualizations were generated using matplotlib^9^ and seaborn^10^ following publication guidelines.

### Computational Infrastructure

Analyses were performed on high-performance computing clusters with 64-96 CPU cores and 128-256 GB RAM per node. Total computational time was approximately 8-12 hours per cell line, with parallel processing reducing wall-clock time to 2-4 hours. Storage requirements included ~50 GB per sample for intermediate files and ~5 GB for final compressed variant files.

### Data Availability and Reproducibility

Analysis pipelines are implemented in modular Python and bash scripts with comprehensive logging and error handling. All parameters and software versions are documented in configuration files to ensure reproducibility. Raw sequencing data are available through the ENCODE Data Coordination Center under appropriate data access agreements. Processed variant files and analysis results are available upon reasonable request.

---

## References

1. Andrews, S. FastQC: a quality control tool for high throughput sequence data. *Babraham Bioinformatics* (2010).

2. Li, H. & Durbin, R. Fast and accurate short read alignment with Burrows-Wheeler transform. *Bioinformatics* **25**, 1754–1760 (2009).

3. Van der Auwera, G. A. et al. From FastQ data to high confidence variant calls: the Genome Analysis Toolkit best practices pipeline. *Curr. Protoc. Bioinformatics* **43**, 11.10.1–11.10.33 (2013).

4. Poplin, R. et al. A universal SNP and small-indel variant caller using deep neural networks. *Nat. Biotechnol.* **36**, 983–987 (2018).

5. Huang, Y. et al. DeepSomatic: accurate somatic small variant calling with deep learning. *bioRxiv* 2023.12.09.571027 (2023).

6. Danecek, P. et al. Twelve years of SAMtools and BCFtools. *Gigascience* **10**, giab008 (2021).

7. McKinney, W. Data structures for statistical computing in Python. *Proc. 9th Python Sci. Conf.* 56–61 (2010).

8. Harris, C. R. et al. Array programming with NumPy. *Nature* **585**, 357–362 (2020).

9. Hunter, J. D. Matplotlib: a 2D graphics environment. *Comput. Sci. Eng.* **9**, 90–95 (2007).

10. Waskom, M. seaborn: statistical data visualization. *J. Open Source Softw.* **6**, 3021 (2021).