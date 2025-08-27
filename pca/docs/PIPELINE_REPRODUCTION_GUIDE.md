# Hybrid Ancestry Inference Pipeline - Reproduction Guide

## 🎯 **Overview**
This guide provides step-by-step instructions to reproduce the hybrid ancestry inference pipeline results, from raw genetic data to final ancestry assignments with 86.8% success rate.

## 📋 **Prerequisites**

### **Software Requirements**
- **Python 3.8+** with packages:
  - pandas, numpy, matplotlib, seaborn
  - scikit-learn (for ML components)
  - scipy (for statistical tests)
- **PLINK v1.90b6.24** or later
- **Memory**: Minimum 32GB RAM (64GB+ recommended)
- **Storage**: ~500GB for intermediate files

### **Input Data Requirements**
- **Genotype files**: PLINK binary format (.bed/.bim/.fam)
- **Sample metadata**: Donor information with ethnicity data
- **Reference populations**: MAGE 1000 Genomes samples included

## 🔄 **Complete Pipeline Execution**

### **Step 1: Data Preparation**
```bash
# Ensure you're in the WGS PCA directory
cd /path/to/wgs/pca

# Verify input files exist
ls data/plink_files/
# Expected: full_genome_clean.bed, .bim, .fam
```

### **Step 2: Stage 1 - Reference-Guided PCA Clustering**
```bash
# Run Stage 1 analysis
python scripts/analysis/plink_pca_analysis.py
```

**What this does:**
- Loads 2,884,699 ancestry-informative markers
- Performs LD pruning and quality control
- Runs PCA with PLINK v1.90b6.24
- Implements reference-guided clustering using MAGE 1000 Genomes centroids
- Assigns ancestry with confidence scoring
- **Output**: 81.0% success rate (1,888/2,331 samples)

**Key outputs:**
- `results/plink_analysis/plink_pca_results.csv`
- `results/plink_analysis/plink_pca_analysis.png`
- `results/plink_analysis/plink_pca_analysis_report.txt`

### **Step 3: Stage 2 - ML Reference Projection**
```bash
# Run Stage 2 ML inference
python scripts/analysis/reference_projection_analysis.py
```

**What this does:**
- Uses 731 MAGE samples as training references (5 ancestry groups)
- Trains ensemble classifiers (KNN + Random Forest) on 20 PCs
- Projects 443 Stage 1 uncertain samples
- Applies 70% confidence threshold for high-confidence assignments
- **Output**: Additional 117 high-confidence assignments (5.0% improvement)

**Key outputs:**
- `results/reference_projection/reference_projection_results.csv`
- `results/reference_projection/reference_projection_analysis.png`
- `results/reference_projection/reference_projection_report.txt`

### **Step 4: Final Consensus Integration**
```bash
# Run comprehensive combined analysis
python scripts/analysis/combined_ancestry_analysis.py
```

**What this does:**
- Combines Stage 1 and Stage 2 results
- Creates final hybrid consensus ancestry assignments
- Performs concordance analysis with self-reported ethnicity
- Generates comprehensive performance statistics
- **Final Output**: 86.8% overall success rate (2,024/2,331 samples)

**Key outputs:**
- `results/combined_analysis/complete_ancestry_results.csv`
- `results/combined_analysis/ancestry_method_comparison.csv`
- `results/combined_analysis/comprehensive_ancestry_analysis.png`
- `results/combined_analysis/comprehensive_ancestry_analysis_report.txt`

## 🔍 **Optional Validation Analyses**

### **ADNI Population Structure Investigation**
```bash
# Investigate ADNI vs GTEx demographic differences
python scripts/analysis/investigate_adni_ancestry.py
python scripts/analysis/compare_adni_gtex.py
```

**What this reveals:**
- ADNI contains 83.2% Latino/Hispanic participants self-identifying as "white"
- Population admixture explains low concordance (5.2% vs GTEx 96.6%)
- Validates pipeline accuracy in detecting real demographic patterns

### **Concordance Analysis**
```bash
# Detailed ethnicity vs ancestry concordance
python scripts/analysis/investigate_concordance.py
```

## 📊 **Key Result Files**

### **Main Dataset Files**
| File | Description | Location |
|------|-------------|----------|
| `complete_ancestry_results.csv` | Complete dataset with all ancestry assignments | `results/combined_analysis/` |
| `ancestry_method_comparison.csv` | Method-by-method comparison | `results/combined_analysis/` |
| `pipeline_performance_summary.csv` | Overall performance metrics | `results/combined_analysis/` |
| `donor_summary.csv` | Original donor metadata with ancestry | `results/updated_metadata/` |

### **Visualization Files**
| File | Description | Shows |
|------|-------------|--------|
| `comprehensive_ancestry_analysis.png` | Complete pipeline overview | All stages, methods, performance |
| `plink_pca_analysis.png` | Stage 1 results | PCA clustering, 81.0% success |
| `reference_projection_analysis.png` | Stage 2 results | ML projection, confidence dist. |
| `adni_gtex_comparison.png` | Population comparison | Demographic differences |
| `adni_ancestry_investigation.png` | ADNI analysis | Population structure validation |

### **Report Files**
| File | Description | Contains |
|------|-------------|----------|
| `comprehensive_ancestry_analysis_report.txt` | Complete summary | Final performance, all statistics |
| `plink_pca_analysis_report.txt` | Stage 1 report | PCA methodology, results |
| `reference_projection_report.txt` | Stage 2 report | ML performance, confidence |

## ⚡ **Quick Reproduction (Complete Pipeline)**

For users wanting to reproduce the entire pipeline in one command:

```bash
# Navigate to project directory
cd /path/to/wgs/pca

# Run complete pipeline
./scripts/run_complete_pipeline.sh
```

Or manually in sequence:
```bash
python scripts/analysis/plink_pca_analysis.py && \
python scripts/analysis/reference_projection_analysis.py && \
python scripts/analysis/combined_ancestry_analysis.py
```

## 🎯 **Expected Performance Metrics**

### **Stage 1 (Reference-Guided PCA)**
- **Success Rate**: 81.0% (1,888/2,331 samples)
- **Mean Confidence**: 0.690
- **High Confidence (>0.7)**: 1,144 samples (60.6%)
- **AFR Recovery**: 292 samples (vs 2 with fixed boundaries)

### **Stage 2 (ML Projection)**
- **Samples Processed**: 443 Stage 1 uncertain samples
- **High-Confidence Assignments**: 57 samples (confidence >0.7)
- **Medium-Confidence Assignments**: 60 samples (0.5-0.7)
- **ML Classifier Performance**: KNN 99.7%, RF 99.5% (5-fold CV)

### **Final Hybrid Consensus**
- **Overall Success**: 86.8% (2,024/2,331 samples)
- **Population Distribution**: EUR 38.4%, AMR 19.0%, AFR 16.0%, EAS 6.8%, SAS 5.8%
- **Uncertain Samples**: 307 samples (13.2%)
- **Mean Final Confidence**: 0.69

## 🔧 **Troubleshooting**

### **Common Issues**
1. **Memory errors**: Increase available RAM or use compute cluster
2. **Missing files**: Ensure all input data files are present and accessible
3. **Python package errors**: Install required packages with pip/conda
4. **PLINK not found**: Ensure PLINK is in system PATH

### **Performance Optimization**
- **Parallel processing**: Some scripts support multiprocessing
- **Memory management**: Large datasets may need chunking
- **Compute resources**: Use HPC/cloud for large-scale analyses

## 📈 **Validation Metrics**

To validate successful reproduction:

1. **Stage 1 success rate**: Should be ~81.0%
2. **AFR sample count**: Should recover ~292 AFR samples
3. **Final success rate**: Should achieve ~86.8% overall
4. **Population distribution**: Should match reported percentages
5. **Key visualizations**: Should match provided plots

## 📚 **Citation and References**

When using this pipeline, please cite:
- PLINK: Purcell et al. (2007) American Journal of Human Genetics
- 1000 Genomes Project: Nature 526:68-74 (2015)
- HANCESTRO: McMahon et al. (2021) Genome Biology

---

**Last Updated**: December 2024
**Pipeline Version**: v2024.12 (Reference-Guided Hybrid)
**Success Rate**: 86.8% final ancestry assignments