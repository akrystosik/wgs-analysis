#!/bin/bash
# cell_line_config.sh
#
# Purpose: Configuration file for the cell line variant analysis pipeline.
# This file contains all cell-line-specific settings and file paths.
# Keeping configuration separate from implementation makes it easier to:
# - Add new cell lines
# - Update file paths
# - Maintain consistency across analyses
# - Share configurations between different scripts

# Base directories for the project
BASE_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/encode_analysis"
OUTPUT_ROOT="${BASE_DIR}/data/variants/cell_line_analysis_v4"

# Define mappings for DeepSomatic files
declare -A DEEPSOMATIC_FILES
DEEPSOMATIC_FILES["Caki2"]="data/variants/deepsomatic/Caki2_somatic_full.vcf"
DEEPSOMATIC_FILES["K562"]="data/variants/deepsomatic/K562_WGS_pair1.deepsomatic.with_germline.vcf.gz"
# Add new cell lines here following the same pattern
# DEEPSOMATIC_FILES["NewCellLine"]="path/to/file"

# Define mappings for DeepVariant files
declare -A DEEPVARIANT_FILES
DEEPVARIANT_FILES["Caki2"]="data/variants/deepvariant/Caki2_WGS_pair1.vcf.gz"
DEEPVARIANT_FILES["K562"]="data/variants/deepvariant/K562_WGS_pair1.deepvariant.vcf.gz"
# Add new cell lines here following the same pattern
# DEEPVARIANT_FILES["NewCellLine"]="path/to/file"

# Validate that each cell line has both files defined
for cell_line in "${!DEEPSOMATIC_FILES[@]}"; do
    if [[ ! "${DEEPVARIANT_FILES[$cell_line]+isset}" ]]; then
        echo "Error in configuration: ${cell_line} missing DeepVariant file definition"
        exit 1
    fi
done