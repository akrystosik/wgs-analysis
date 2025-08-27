#!/bin/bash

# Create output directory
OUTPUT_DIR="data/variants/encode"
mkdir -p $OUTPUT_DIR
TEMP_DIR="$OUTPUT_DIR/temp"
mkdir -p $TEMP_DIR

# Array of ENCODE file accessions
ENCODE_FILES=(
    "ENCFF752OAX"
    "ENCFF785JVR"
    "ENCFF574MDJ"
    "ENCFF863MPP"
)

# Download and process function
process_encode_file() {
    accession=$1
    echo "Processing ${accession}..."
    
    # Download to temporary location
    echo "Downloading ${accession}.vcf.gz..."
    url="https://www.encodeproject.org/files/${accession}/@@download/${accession}.vcf.gz"
    wget -c -O "$TEMP_DIR/${accession}.vcf.gz" "$url"
    
    # Decompress and recompress with bgzip
    echo "Recompressing with bgzip..."
    gunzip -c "$TEMP_DIR/${accession}.vcf.gz" | bgzip > "$OUTPUT_DIR/${accession}.vcf.gz"
    
    # Index with tabix
    echo "Creating index..."
    tabix -p vcf "$OUTPUT_DIR/${accession}.vcf.gz"
    
    # Verify the file
    echo "Verifying ${accession}..."
    if bcftools view -h "$OUTPUT_DIR/${accession}.vcf.gz" >/dev/null; then
        echo "${accession} processed successfully"
    else
        echo "Error processing ${accession}"
        return 1
    fi
}

# Check for required tools
for tool in wget gunzip bgzip tabix bcftools; do
    if ! command -v $tool &> /dev/null; then
        echo "Error: $tool is required but not installed."
        exit 1
    fi
done

# Process each file
for file in "${ENCODE_FILES[@]}"; do
    process_encode_file "$file"
done

# Clean up
rm -rf $TEMP_DIR

echo "All files processed. Checking final status..."
# Verify final files
for file in "${ENCODE_FILES[@]}"; do
    if [ -f "${OUTPUT_DIR}/${file}.vcf.gz" ] && [ -f "${OUTPUT_DIR}/${file}.vcf.gz.tbi" ]; then
        if bcftools view -h "${OUTPUT_DIR}/${file}.vcf.gz" >/dev/null; then
            echo "${file}: OK"
        else
            echo "${file}: Failed verification"
        fi
    else
        echo "${file}: Missing files"
    fi
done