#!/usr/bin/env python3
"""
Validate VCF Croissant JSON-LD metadata files.

This script validates the generated VCF Croissant metadata files against the 
official schema using the mlcroissant library, with special handling for 
both open access and protected dataset metadata.
"""

import mlcroissant as mlc
from pathlib import Path
from typing import List
import json
import sys

def validate_vcf_croissant_file(file_path: Path) -> bool:
    """Validate a single VCF Croissant metadata file."""
    
    print(f"\nValidating: {file_path.name}")
    print(f"{'='*60}")
    
    try:
        # Load and validate the metadata
        dataset = mlc.Dataset(str(file_path))
        
        # If we get here, basic validation passed
        print("✓ Basic JSON-LD structure is valid")
        
        # Check metadata fields
        metadata = dataset.metadata
        print(f"✓ Dataset name: {metadata.name}")
        print(f"✓ Dataset description length: {len(metadata.description)} characters")
        
        if hasattr(metadata, 'citation') and metadata.citation:
            print(f"✓ Citations: {len(metadata.citation)} found")
        
        if hasattr(metadata, 'license'):
            print(f"✓ License: {metadata.license}")
        
        # Check distribution files
        if hasattr(metadata, 'distribution') and metadata.distribution:
            print(f"✓ Distribution files: {len(metadata.distribution)} found")
            
            # Check for VCF-specific fields
            for dist in metadata.distribution:
                if hasattr(dist, 'encoding_format'):
                    print(f"  - Encoding format: {dist.encoding_format}")
        
        # Check for genomics-specific fields by loading raw JSON
        with open(file_path, 'r') as f:
            raw_data = json.load(f)
        
        genomics_fields = [
            'sc:organism', 'sc:assayType', 'sc:dataType', 'sc:accessLevel'
        ]
        
        for field in genomics_fields:
            if field in raw_data:
                print(f"✓ {field}: {raw_data[field]}")
        
        # Check access level specific fields
        access_level = raw_data.get('sc:accessLevel', 'unknown')
        if access_level == 'protected':
            print("✓ Protected dataset - source-only metadata")
            if 'sc:dataAccess' in raw_data:
                print(f"  - Data access URL: {raw_data['sc:dataAccess']}")
        elif access_level == 'open':
            print("✓ Open access dataset - file metadata included")
            if 'sc:datasetStatistics' in raw_data:
                stats = raw_data['sc:datasetStatistics']
                if 'numberOfVariants' in stats:
                    print(f"  - Total variants: {stats['numberOfVariants']:,}")
                if 'numberOfSamples' in stats:
                    print(f"  - Total samples: {stats['numberOfSamples']}")
        
        print("✓ VCF Croissant validation successful!")
        return True
        
    except Exception as e:
        print(f"✗ Validation failed: {str(e)}")
        
        # Try to load as raw JSON to check for basic JSON validity
        try:
            with open(file_path, 'r') as f:
                json.load(f)
            print("✓ File is valid JSON")
        except json.JSONDecodeError as json_e:
            print(f"✗ JSON parsing error: {str(json_e)}")
        
        return False

def analyze_dataset_coverage(croissant_files: List[Path]) -> None:
    """Analyze coverage of different dataset types."""
    
    print(f"\n{'='*60}")
    print("Dataset Coverage Analysis")
    print(f"{'='*60}")
    
    datasets = {'open': [], 'protected': []}
    
    for file_path in croissant_files:
        try:
            with open(file_path, 'r') as f:
                data = json.load(f)
            
            access_level = data.get('sc:accessLevel', 'unknown')
            dataset_name = data.get('name', 'Unknown')
            
            if access_level in datasets:
                datasets[access_level].append({
                    'name': dataset_name,
                    'file': file_path.name,
                    'organism': data.get('sc:organism', 'Unknown'),
                    'assay': data.get('sc:assayType', 'Unknown')
                })
        except Exception as e:
            print(f"Error analyzing {file_path.name}: {e}")
    
    print(f"Open Access Datasets ({len(datasets['open'])}):")
    for dataset in datasets['open']:
        print(f"  ✓ {dataset['name']}")
        print(f"    - File: {dataset['file']}")
        print(f"    - Organism: {dataset['organism']}")
        print(f"    - Assay: {dataset['assay']}")
    
    print(f"\nProtected Datasets ({len(datasets['protected'])}):")
    for dataset in datasets['protected']:
        print(f"  🔒 {dataset['name']}")
        print(f"    - File: {dataset['file']} (source-only)")
        print(f"    - Organism: {dataset['organism']}")
        print(f"    - Assay: {dataset['assay']}")

def main():
    """Main function to validate all VCF Croissant files."""
    
    # Find all VCF Croissant files
    metadata_dir = Path('../croissant_metadata')
    if not metadata_dir.exists():
        metadata_dir = Path('./croissant_metadata')
    if not metadata_dir.exists():
        metadata_dir = Path('.')
    
    croissant_files = list(metadata_dir.glob('*_wgs_variants_croissant.jsonld'))
    
    if not croissant_files:
        print("No VCF Croissant metadata files found")
        return
    
    print(f"Found {len(croissant_files)} VCF Croissant metadata files to validate")
    
    valid_files = 0
    total_files = len(croissant_files)
    
    for file_path in sorted(croissant_files):
        if validate_vcf_croissant_file(file_path):
            valid_files += 1
    
    # Analyze dataset coverage
    analyze_dataset_coverage(croissant_files)
    
    print(f"\n{'='*60}")
    print(f"VCF Validation Summary: {valid_files}/{total_files} files passed validation")
    
    if valid_files == total_files:
        print("🎉 All VCF Croissant metadata files are valid!")
        sys.exit(0)
    else:
        print(f"⚠️  {total_files - valid_files} files failed validation")
        sys.exit(1)

if __name__ == "__main__":
    main()