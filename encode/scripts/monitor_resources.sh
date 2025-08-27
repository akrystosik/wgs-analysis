#!/bin/bash
# Monitor system resources during pipeline execution

OUTPUT_DIR="/mnt/czi-sci-ai/intrinsic-variation-gene-ex/encode_analysis/logs/resources"
mkdir -p "$OUTPUT_DIR"

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
OUTPUT_FILE="$OUTPUT_DIR/resources_${TIMESTAMP}.log"

echo "Starting resource monitoring. Log: $OUTPUT_FILE"
echo "Press Ctrl+C to stop monitoring."

echo "Timestamp,CPU_Usage(%),Memory_Used(GB),Memory_Total(GB),Disk_Used(GB),Disk_Total(GB)" > "$OUTPUT_FILE"

while true; do
    TIMESTAMP=$(date +%Y-%m-%d_%H:%M:%S)
    CPU_USAGE=$(top -bn1 | grep "Cpu(s)" | awk '{print $2}')
    MEM_TOTAL=$(free -g | grep Mem | awk '{print $2}')
    MEM_USED=$(free -g | grep Mem | awk '{print $3}')
    DISK_TOTAL=$(df -BG /mnt/czi-sci-ai | tail -1 | awk '{print $2}' | tr -d 'G')
    DISK_USED=$(df -BG /mnt/czi-sci-ai | tail -1 | awk '{print $3}' | tr -d 'G')
    
    echo "$TIMESTAMP,$CPU_USAGE,$MEM_USED,$MEM_TOTAL,$DISK_USED,$DISK_TOTAL" >> "$OUTPUT_FILE"
    sleep 60
done
