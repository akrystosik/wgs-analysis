#!/bin/bash

# Script to process GTEx samples in batches with full sample IDs
# Usage: ./process_gtex_batches_newvcf.sh [batch_size] [max_parallel] [start_index]

# Settings with default values
BATCH_SIZE=${1:-50}      # Number of samples per batch
MAX_PARALLEL=${2:-32}    # Number of parallel processes
START_INDEX=${3:-0}      # Start index (skip previously processed samples)

# Define paths
VCF="/mnt/czi-sci-ai/intrinsic-variation-gene-ex-3/project_gene_regulation/data/GTEx/WGS/phg001796.v1.GTEx_v9_WGS_953.genotype-calls-vcf.c1/GTEx_Analysis_2021-02-11_v9_WholeGenomeSeq_953Indiv.vcf.gz"
OUTPUT_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/wgs/data/qc/variants/gtex_counts_newvcf"
TEMP_DIR="${OUTPUT_DIR}/temp"
SUMMARY_FILE="${OUTPUT_DIR}/gtex_summary_complete.tsv"
ALL_SAMPLES_FILE="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/wgs/data/qc/variants/gtex_samples.txt"
PROCESSED_FILE="${OUTPUT_DIR}/processed_samples.txt"

# Create necessary directories
mkdir -p "$OUTPUT_DIR"
mkdir -p "$TEMP_DIR"

echo "GTEx Batch Processor (New VCF)"
echo "=========================="
echo "VCF file: $VCF"
echo "Output directory: $OUTPUT_DIR"
echo "Batch size: $BATCH_SIZE"
echo "Max parallel processes: $MAX_PARALLEL"
echo "Starting from index: $START_INDEX"

# Check if sample file exists
if [ ! -f "$ALL_SAMPLES_FILE" ]; then
    echo "Error: Sample ID file not found at $ALL_SAMPLES_FILE"
    exit 1
fi

TOTAL_SAMPLES=$(wc -l < "$ALL_SAMPLES_FILE")
echo "Found $TOTAL_SAMPLES total samples"

# Initialize or load processed samples file
if [ ! -f "$PROCESSED_FILE" ]; then
    # If we're starting fresh, create empty file
    touch "$PROCESSED_FILE"
    
    # If we're starting from a specific index, mark previous samples as processed
    if [ "$START_INDEX" -gt 0 ]; then
        head -n "$START_INDEX" "$ALL_SAMPLES_FILE" >> "$PROCESSED_FILE"
        echo "Marked first $START_INDEX samples as processed"
    fi
fi

# Count already processed samples
PROCESSED_COUNT=$(wc -l < "$PROCESSED_FILE")
echo "Already processed $PROCESSED_COUNT samples"

# Calculate remaining samples
REMAINING=$((TOTAL_SAMPLES - PROCESSED_COUNT))
echo "Remaining samples to process: $REMAINING"

# Initialize or load the complete summary file
if [ ! -f "$SUMMARY_FILE" ]; then
    # Create header
    echo -e "Sample\tSNPs\tIndels\tTotal" > "$SUMMARY_FILE"
fi

# Function to process a single batch
process_batch() {
    local batch_start=$1
    local batch_size=$2
    local batch_end=$((batch_start + batch_size - 1))
    
    if [ "$batch_end" -ge "$TOTAL_SAMPLES" ]; then
        batch_end=$((TOTAL_SAMPLES - 1))
    fi
    
    echo "Processing batch from index $batch_start to $batch_end"
    
    # Get sample IDs for this batch
    local batch_file="${TEMP_DIR}/batch_${batch_start}_${batch_end}.txt"
    sed -n "$((batch_start + 1)),$((batch_end + 1))p" "$ALL_SAMPLES_FILE" > "$batch_file"
    
    # Create a batch-specific summary file
    local batch_summary="${TEMP_DIR}/summary_${batch_start}_${batch_end}.tsv"
    echo -e "Sample\tSNPs\tIndels\tTotal" > "$batch_summary"
    
    # Process each sample in the batch
    local active_processes=0
    while read -r sample; do
        # Check if already processed
        if grep -q "^$sample$" "$PROCESSED_FILE"; then
            echo "Skipping $sample (already processed)"
            continue
        fi
        
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
            
            # Add to batch summary file (with file lock)
            flock "$batch_summary" bash -c "echo -e \"$sample\t$total_snps\t$total_indels\t$((total_snps + total_indels))\" >> \"$batch_summary\""
            
            # Mark as processed
            flock "$PROCESSED_FILE" bash -c "echo \"$sample\" >> \"$PROCESSED_FILE\""
            
            echo "Completed sample $sample: $total_snps SNPs, $total_indels indels"
        ) &
        
        # Update active process count
        active_processes=$(jobs -p | wc -l)
        echo "Current active processes: $active_processes"
        
    done < "$batch_file"
    
    # Wait for all processes in this batch to complete
    echo "Waiting for batch $batch_start-$batch_end to complete..."
    wait
    
    # Append batch results to main summary
    tail -n +2 "$batch_summary" >> "$SUMMARY_FILE"
    
    echo "Batch $batch_start-$batch_end complete"
}

# Process remaining samples in batches
for ((i = PROCESSED_COUNT; i < TOTAL_SAMPLES; i += BATCH_SIZE)); do
    process_batch $i $BATCH_SIZE
    
    # Calculate progress
    current_processed=$(wc -l < "$PROCESSED_FILE")
    progress=$((100 * current_processed / TOTAL_SAMPLES))
    echo "Overall progress: $progress% ($current_processed/$TOTAL_SAMPLES samples processed)"
done

# Calculate final statistics
echo "Calculating final statistics..."
if command -v awk &> /dev/null; then
    # Add summary row with averages
    awk 'BEGIN {FS="\t"; snps=0; indels=0; total=0; count=0} 
        NR>1 {snps+=$2; indels+=$3; total+=$4; count++} 
        END {if(count>0) printf "AVERAGE\t%.2f\t%.2f\t%.2f\n", snps/count, indels/count, total/count}' "$SUMMARY_FILE" > "${TEMP_DIR}/average.tsv"
    
    # Replace existing average line or append
    if grep -q "AVERAGE" "$SUMMARY_FILE"; then
        sed -i '/AVERAGE/d' "$SUMMARY_FILE"
    fi
    cat "${TEMP_DIR}/average.tsv" >> "$SUMMARY_FILE"
    
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
rm -rf "${TEMP_DIR}"