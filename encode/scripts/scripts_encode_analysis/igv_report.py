#!/usr/bin/env python3

import os
import json
from argparse import ArgumentParser

def create_igv_report_json(bam_path, vcf_path, bed_path, fasta_path, output_dir):
    """Create JSON configuration for IGV report"""
    config = {
        "reference": {
            "id": "hg38",
            "name": "Human (GRCh38/hg38)",
            "fastaURL": fasta_path,
            "indexed": True
        },
        "tracks": [
            {
                "name": "K562 BAM",
                "type": "alignment",
                "format": "bam",
                "url": bam_path,
                "indexed": True,
                "height": 500
            },
            {
                "name": "DeepSomatic Variants",
                "type": "variant",
                "format": "vcf",
                "url": vcf_path,
                "indexed": True,
                "height": 50
            },
            {
                "name": "Regions of Interest",
                "type": "annotation",
                "format": "bed",
                "url": bed_path,
                "indexed": False,
                "displayMode": "EXPANDED",
                "height": 50
            }
        ]
    }
    
    # Write config
    output_file = os.path.join(output_dir, "igv_report_config.json")
    with open(output_file, 'w') as f:
        json.dump(config, f, indent=2)
    
    return output_file

def main():
    parser = ArgumentParser(description='Generate IGV report configuration')
    parser.add_argument('--bam', required=True, help='Path to BAM file')
    parser.add_argument('--vcf', required=True, help='Path to VCF file')
    parser.add_argument('--bed', required=True, help='Path to BED file with regions')
    parser.add_argument('--fasta', required=True, help='Path to reference FASTA')
    parser.add_argument('--output', required=True, help='Output directory')
    
    args = parser.parse_args()
    
    # Create output directory if it doesn't exist
    os.makedirs(args.output, exist_ok=True)
    
    # Create IGV report config
    config_file = create_igv_report_json(
        args.bam,
        args.vcf,
        args.bed,
        args.fasta,
        args.output
    )
    
    print(f"Created IGV report configuration: {config_file}")
    print("\nTo create the IGV report, run:")
    print(f"create_report {args.fasta} {config_file} {os.path.join(args.output, 'k562_report.html')} -t {args.bed}")

if __name__ == '__main__':
    main()