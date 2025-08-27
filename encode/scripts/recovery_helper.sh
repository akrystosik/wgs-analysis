#!/bin/bash
# Helper script to recover from common pipeline errors

function fix_bam_index() {
    local bam_file="$1"
    echo "Fixing BAM index for $bam_file"
    samtools index "$bam_file"
}

function fix_vcf_index() {
    local vcf_file="$1"
    echo "Fixing VCF index for $vcf_file"
    bcftools index -t "$vcf_file"
}

function fix_permissions() {
    local directory="$1"
    echo "Fixing permissions for $directory"
    chmod -R u+rw "$directory"
}
