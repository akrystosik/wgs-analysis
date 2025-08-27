#!/bin/bash
#/mnt/czi-sci-ai/intrinsic-variation-gene-ex/wgs/scripts/process_gtex_samples_simple.sh
# Script to count non-reference variants in multiple GTEx samples
# Usage: ./process_gtex_samples_simple.sh [number_of_samples] [max_parallel_processes]

# Settings with default values
NUM_SAMPLES=${1:-10}  # Default to 10 samples if not specified
MAX_PARALLEL=${2:-4}  # Default to 4 parallel processes
RANDOM_SELECTION=true # Set to true for random sample selection, false for first N

# Define paths
VCF="/mnt/czi-sci-ai/intrinsic-variation-gene-ex-3/project_gene_regulation/data/GTEx/WGS/phg001796.v1.GTEx_v9_WGS_phased.genotype-calls-vcf.c1/GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv.SHAPEIT2_phased.vcf.gz"
OUTPUT_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/wgs/data/qc/variants/gtex_counts"
TEMP_DIR="${OUTPUT_DIR}/temp"
SUMMARY_FILE="${OUTPUT_DIR}/gtex_summary.tsv"

# Create output directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$TEMP_DIR"

echo "GTEx Variant Counter"
echo "===================="
echo "VCF file: $VCF"
echo "Output directory: $OUTPUT_DIR"
echo "Processing $NUM_SAMPLES samples with max $MAX_PARALLEL parallel processes"

# Get all sample IDs
echo "Retrieving sample IDs..."
bcftools query -l "$VCF" > "${TEMP_DIR}/all_samples.txt"
TOTAL_SAMPLES=$(wc -l < "${TEMP_DIR}/all_samples.txt")
echo "Found $TOTAL_SAMPLES total samples"

# Select samples to process
if [ "$RANDOM_SELECTION" = true ] && [ "$NUM_SAMPLES" -lt "$TOTAL_SAMPLES" ]; then
    echo "Randomly selecting $NUM_SAMPLES samples..."
    # If shuf is available
    if command -v shuf &> /dev/null; then
        shuf -n "$NUM_SAMPLES" "${TEMP_DIR}/all_samples.txt" > "${TEMP_DIR}/selected_samples.txt"
    else
        # Simple random selection if shuf is not available
        sort -R "${TEMP_DIR}/all_samples.txt" | head -n "$NUM_SAMPLES" > "${TEMP_DIR}/selected_samples.txt"
    fi
else
    echo "Selecting first $NUM_SAMPLES samples..."
    head -n "$NUM_SAMPLES" "${TEMP_DIR}/all_samples.txt" > "${TEMP_DIR}/selected_samples.txt"
fi

# Create header for summary file
echo -e "Sample\tSNPs\tIndels\tTotal" > "$SUMMARY_FILE"

# Process each sample
active_processes=0
while read -r sample; do
    # Wait if max parallel jobs reached
    while [ $active_processes -ge $MAX_PARALLEL ]; do
        active_processes=$(jobs -p | wc -l)
        sleep 1
    done
    
    # Start a background process for this sample
    (
        sample_output_dir="${OUTPUT_DIR}/${sample}"
        mkdir -p "$sample_output_dir"
        
        echo "Processing sample: $sample"
        
        # Initialize counters
        total_snps=0
        total_indels=0
        
        # Process each chromosome
        for chr in {1..22} X; do
            echo "  Processing chromosome chr${chr} for $sample..."
            
            # Count SNPs
            bcftools view -s "$sample" -r "chr${chr}" -v snps "$VCF" > "${TEMP_DIR}/${sample}_${chr}_snps.vcf"
            snp_count=$(grep -v "^#" "${TEMP_DIR}/${sample}_${chr}_snps.vcf" | grep -v "0|0" | grep -v "\\.|\\." | wc -l)
            
            # Count indels
            bcftools view -s "$sample" -r "chr${chr}" -v indels "$VCF" > "${TEMP_DIR}/${sample}_${chr}_indels.vcf"
            indel_count=$(grep -v "^#" "${TEMP_DIR}/${sample}_${chr}_indels.vcf" | grep -v "0|0" | grep -v "\\.|\\." | wc -l)
            
            # Remove temporary files
            rm -f "${TEMP_DIR}/${sample}_${chr}_snps.vcf" "${TEMP_DIR}/${sample}_${chr}_indels.vcf"
            
            # Save chromosome-specific counts
            echo -e "Chromosome\tSNPs\tIndels" > "${sample_output_dir}/chr${chr}.txt"
            echo -e "chr${chr}\t${snp_count}\t${indel_count}" >> "${sample_output_dir}/chr${chr}.txt"
            
            # Update totals
            total_snps=$((total_snps + snp_count))
            total_indels=$((total_indels + indel_count))
            
            echo "    ${sample}: chr${chr} - SNPs: $snp_count, Indels: $indel_count"
        done
        
        # Save total counts
        echo -e "Type\tCount" > "${sample_output_dir}/total.txt"
        echo -e "SNPs\t$total_snps" >> "${sample_output_dir}/total.txt"
        echo -e "Indels\t$total_indels" >> "${sample_output_dir}/total.txt"
        echo -e "Total\t$((total_snps + total_indels))" >> "${sample_output_dir}/total.txt"
        
        # Add to summary file (with file lock to prevent corruption)
        flock "$SUMMARY_FILE" bash -c "echo -e \"$sample\t$total_snps\t$total_indels\t$((total_snps + total_indels))\" >> \"$SUMMARY_FILE\""
        
        echo "Completed sample $sample: $total_snps SNPs, $total_indels indels"
    ) &
    
    # Update active process count
    active_processes=$(jobs -p | wc -l)
    echo "Current active processes: $active_processes"
    
done < "${TEMP_DIR}/selected_samples.txt"

# Wait for all background processes to complete
echo "Waiting for all processes to complete..."
wait

echo "All sample processing jobs completed."

# Calculate summary statistics
echo "Calculating summary statistics..."
if command -v awk &> /dev/null; then
    # Add summary row with averages
    awk 'BEGIN {FS="\t"; snps=0; indels=0; total=0; count=0} 
        NR>1 {snps+=$2; indels+=$3; total+=$4; count++} 
        END {if(count>0) printf "AVERAGE\t%.2f\t%.2f\t%.2f\n", snps/count, indels/count, total/count}' "$SUMMARY_FILE" >> "$SUMMARY_FILE"
    
    # Extract averages for display
    read -r avg_snps avg_indels avg_total < <(awk 'BEGIN {FS="\t"} $1=="AVERAGE" {print $2, $3, $4}' "$SUMMARY_FILE")
else
    avg_snps="N/A"
    avg_indels="N/A"
    avg_total="N/A"
fi

echo "Processing complete!"
echo "Processed $(grep -v "AVERAGE" "$SUMMARY_FILE" | wc -l) samples"
echo "Average counts: SNPs = $avg_snps, Indels = $avg_indels, Total = $avg_total"
echo "Summary file: $SUMMARY_FILE"

# Clean up temporary directory
rm -rf "$TEMP_DIR"