#!/bin/bash
#bash scripts/02_bwa_alignment_03_bam_qc.sh
# Processing Pair 2 from FASTQ to QC

# Define paths
PAIR_ID="K562_WGS_pair2"
FASTQ_R1="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/encode_analysis/data/fastq/K562/${PAIR_ID}_R1.fastq.gz"
FASTQ_R2="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/encode_analysis/data/fastq/K562/${PAIR_ID}_R2.fastq.gz"
REF="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/encode_analysis/data/reference/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
OUTPUT_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/encode_analysis/data/bam/K562"
QC_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/encode_analysis/data/qc"

# BWA Alignment
bwa mem -t 96 \
  -R "@RG\tID:${PAIR_ID}\tSM:K562\tLB:WGS\tPL:ILLUMINA" \
  ${REF} \
  ${FASTQ_R1} \
  ${FASTQ_R2} | \
  samtools view -b - > ${OUTPUT_DIR}/${PAIR_ID}.bam

# Verify BAM
if [[ ! -s ${OUTPUT_DIR}/${PAIR_ID}.bam ]]; then
  echo "Error: Alignment failed or generated empty BAM file."
  exit 1
fi

# Sort BAM
samtools sort -@ 96 -o ${OUTPUT_DIR}/${PAIR_ID}.sorted.bam ${OUTPUT_DIR}/${PAIR_ID}.bam

# Verify Sorted BAM
if [[ ! -s ${OUTPUT_DIR}/${PAIR_ID}.sorted.bam ]]; then
  echo "Error: Sorting failed or generated empty BAM file."
  exit 1
fi

# Mark Duplicates
java -Xmx128G -jar picard.jar MarkDuplicates \
  I=${OUTPUT_DIR}/${PAIR_ID}.sorted.bam \
  O=${OUTPUT_DIR}/${PAIR_ID}.marked_duplicates.bam \
  M=${OUTPUT_DIR}/${PAIR_ID}.duplicate_metrics.txt

# Verify Marked Duplicates
if [[ ! -s ${OUTPUT_DIR}/${PAIR_ID}.marked_duplicates.bam ]]; then
  echo "Error: MarkDuplicates failed or generated empty BAM file."
  exit 1
fi

# Index BAM
samtools index ${OUTPUT_DIR}/${PAIR_ID}.marked_duplicates.bam

# Flagstat QC
samtools flagstat ${OUTPUT_DIR}/${PAIR_ID}.marked_duplicates.bam > ${QC_DIR}/${PAIR_ID}.flagstat.txt

# Insert Size Metrics
java -Xmx32G -jar picard.jar CollectInsertSizeMetrics \
  I=${OUTPUT_DIR}/${PAIR_ID}.marked_duplicates.bam \
  O=${QC_DIR}/${PAIR_ID}.insert_size_metrics.txt \
  H=${QC_DIR}/${PAIR_ID}.insert_size_histogram.pdf
