#!/bin/bash

# Set up archive directory with timestamp
ARCHIVE_DIR="archive/k562_analysis_$(date +%Y%m%d)"
mkdir -p $ARCHIVE_DIR

# List of stale directories to archive
STALE_DIRS=(
    "k562_comparison"
    "k562_encode_comparison"
    "k562_improved_comparison"
    "k562_validation_check"
    "test_output"
    "test_isec"
    "consensus_check"
    "encode_lifted"
)

# List of stale scripts to archive
STALE_SCRIPTS=(
    "k562-regions.sh"
    "k562-region-filter.sh"
    "download-encode-vcfs.sh"
    "k562-encode-compare.sh"
    "double_check.sh"
    "liftover_vcf.sh"
    "debug-compare.sh"
)

echo "Starting cleanup..."

# Handle NFS locks
echo "Syncing filesystem..."
sync
sleep 5

# Archive directories
for dir in "${STALE_DIRS[@]}"; do
    if [ -d "$dir" ]; then
        echo "Moving $dir to archive..."
        # Use rsync to handle NFS issues
        rsync -av --remove-source-files "$dir/" "$ARCHIVE_DIR/$dir/"
        rm -rf "$dir"
    fi
done

# Archive scripts
mkdir -p "$ARCHIVE_DIR/scripts"
for script in "${STALE_SCRIPTS[@]}"; do
    if [ -f "scripts/$script" ]; then
        echo "Moving script $script to archive..."
        mv "scripts/$script" "$ARCHIVE_DIR/scripts/"
    fi
done

# Create manifest
{
    echo "K562 Analysis Archive - $(date)"
    echo "=============================="
    echo
    echo "Archived Directories:"
    ls -R "$ARCHIVE_DIR"
    echo
    echo "Current Working Scripts:"
    ls -l scripts/
} > "$ARCHIVE_DIR/manifest.txt"

echo "Cleanup complete! Archive created at $ARCHIVE_DIR"