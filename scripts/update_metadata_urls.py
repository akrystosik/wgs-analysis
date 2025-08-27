#!/usr/bin/env python3
"""
Update Croissant metadata files with clean GitHub repository URLs.

This script updates the multi-omics WGS Croissant JSON-LD metadata files
to use clean, short URLs for the wgs-analysis GitHub repository.

Version: 2.0 - Updated for clean repository separation by data type
"""

import json
import re
from pathlib import Path
from datetime import datetime

def create_clean_dataset_id(filename: str) -> str:
    """Create clean, short dataset identifier."""
    
    # Clean name mapping for multi-omics datasets
    name_map = {
        'adni_wgs_variants_croissant.jsonld': 'adni_wgs_variants',
        'encode_paired_wgs_variants_croissant.jsonld': 'encode_paired_wgs_variants',
        'gtex_wgs_variants_croissant.jsonld': 'gtex_wgs_variants',
        'mage_wgs_variants_croissant.jsonld': 'mage_wgs_variants'
    }
    
    return name_map.get(filename, filename.replace('_croissant.jsonld', ''))

def update_metadata_urls(file_path: Path, github_repo_url: str) -> bool:
    """Update URLs in a single metadata file for clean repository structure."""
    
    print(f"Updating: {file_path.name}")
    
    try:
        # Load the metadata file
        with open(file_path, 'r') as f:
            data = json.load(f)
        
        # Create clean dataset identifier
        clean_id = create_clean_dataset_id(file_path.name)
        
        # Update @id with clean URL structure
        data['@id'] = f"{github_repo_url}/metadata/{clean_id}"
        
        # Update contentUrl for distribution files
        if 'distribution' in data:
            for dist in data['distribution']:
                if 'contentUrl' in dist:
                    current_url = dist['contentUrl']
                    
                    # For open access datasets, update to GitHub raw content
                    if data.get('sc:accessLevel') == 'open':
                        # Extract filename from current URL
                        filename_match = re.search(r'([^/]+\.vcf\.gz)$', current_url)
                        if filename_match:
                            filename = filename_match.group(1)
                            dist['contentUrl'] = f"{github_repo_url}/raw/main/data/{filename}"
                    
                    # For protected datasets, keep the official repository URLs
                    # (ADNI, GTEx, MAGE should keep their official access portals)
        
        # Update sc:dataAccess if it's using example.com or old URLs
        if 'sc:dataAccess' in data and ('example.com' in data['sc:dataAccess'] or 'encode_analysis' in data['sc:dataAccess']):
            if data.get('sc:accessLevel') == 'open':
                data['sc:dataAccess'] = f"{github_repo_url}/data/"
        
        # Update datePublished to current timestamp
        data['datePublished'] = datetime.now().isoformat()
        
        # Add repository information for better discoverability
        data['sc:repository'] = {
            '@type': 'SoftwareSourceCode',
            'name': 'wgs-analysis',
            'url': github_repo_url,
            'description': 'Multi-omics WGS analysis pipeline with paired RNA-seq datasets'
        }
        
        # Write back the updated file
        with open(file_path, 'w') as f:
            json.dump(data, f, indent=2, ensure_ascii=False)
        
        print(f"  ✓ Updated @id: {data['@id']}")
        return True
        
    except Exception as e:
        print(f"  ✗ Error updating {file_path.name}: {e}")
        return False

def main():
    """Main function to update all multi-omics WGS metadata files."""
    
    # Configuration for clean wgs-analysis repository
    github_repo_url = "https://github.com/akrystosik/wgs-analysis"
    metadata_dir = Path("metadata/croissant_metadata")
    
    if not metadata_dir.exists():
        print(f"Error: Metadata directory not found: {metadata_dir}")
        return
    
    # Find multi-omics WGS metadata files
    metadata_files = list(metadata_dir.glob("*_wgs_*croissant.jsonld"))
    
    if not metadata_files:
        print("No WGS Croissant metadata files found")
        return
    
    print("🧬 Multi-omics WGS Metadata URL Update (v2.0)")
    print("=" * 60)
    print(f"Target repository: {github_repo_url}")
    print(f"Found {len(metadata_files)} multi-omics datasets")
    print("-" * 60)
    
    updated_count = 0
    
    for file_path in sorted(metadata_files):
        if update_metadata_urls(file_path, github_repo_url):
            updated_count += 1
    
    print("-" * 60)
    print(f"Updated {updated_count}/{len(metadata_files)} files")
    
    if updated_count == len(metadata_files):
        print("\n🎉 All multi-omics WGS metadata files updated!")
        print("\nClean URLs created:")
        for file_path in sorted(metadata_files):
            clean_id = create_clean_dataset_id(file_path.name)
            print(f"  • {github_repo_url}/metadata/{clean_id}.jsonld")
        
        print(f"\nNext steps:")
        print("1. Validate updated metadata files")
        print("2. Create wgs-analysis GitHub repository") 
        print("3. Push clean, organized content")
        print("4. Update documentation references")
    else:
        print(f"⚠️ {len(metadata_files) - updated_count} files had issues")

if __name__ == "__main__":
    main()