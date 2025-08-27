#!/usr/bin/env python3
"""
Code Organization Plan
Creates a proper organization plan for the repository based on actual usage
"""

import os
import shutil
from pathlib import Path

def create_organization_plan():
    """Create a comprehensive organization plan"""
    
    plan = {
        'keep_active': {
            'description': 'Core files that are actively used and essential for the pipeline',
            'files': [
                # Main analysis scripts (our completed work)
                'pca/scripts/analysis/plink_pca_analysis.py',
                'pca/scripts/analysis/reference_projection_analysis.py', 
                'pca/scripts/utils/metadata_schema_updater.py',
                'pca/scripts/testing/comprehensive_sanity_checks.py',
                'pca/scripts/testing/final_validation.py',
                
                # Main runner
                'pca/run_pca_analysis.py',
                
                # GitHub preparation
                'prepare_for_github.py',
                'code_review_analysis.py',
            ]
        },
        
        'keep_supporting': {
            'description': 'Supporting files that provide utility functions',
            'files': [
                'pca/scripts/utils/enhanced_ancestry_loader.py',
                'pca/scripts/core/real_vcf_pca_analysis.py',  # Used by main analysis
            ]
        },
        
        'archive_experimental': {
            'description': 'Experimental/development files - archive for reference',
            'files': [
                'pca/scripts/experimental/quick_vcf_pca.py',
                'pca/scripts/experimental/real_vcf_direct_pca.py',
                'pca/scripts/experimental/simple_vcf_pca.py',
                'pca/scripts/analysis/enhanced_pca_analysis.py',  # Superseded by PLINK analysis
                'pca/scripts/core/vcf_pca_analysis.py',
                'pca/scripts/core/vcf_pca_analysis_optimized.py',
            ]
        },
        
        'archive_legacy_pipeline': {
            'description': 'Legacy variant calling pipeline - archive as it\'s complete',
            'files': [
                'pipeline/',  # Entire pipeline directory
            ]
        },
        
        'archive_analysis_scripts': {
            'description': 'Analysis scripts from variant calling phase - archive',
            'files': [
                'scripts/analyze_filtered_variants.py',
                'scripts/bcftools_bychromosome.py',
                'scripts/bcftools_stats.py',
                'scripts/compare_matched_variants.py',
                'scripts/compare_variants.py',
                'scripts/deepmap.py',
                'scripts/direct_count.py',
                'scripts/gtex_extraction.py',
                'scripts/igv_report.py',
                'scripts/investigate_high_af_discord.py',
                'scripts/investigate_variant_discordance.py',
                'scripts/summarize_data_by_type.py',
                'scripts/targeted_analysis.py',
            ]
        },
        
        'archive_utils': {
            'description': 'Utility scripts - archive as they were development helpers',
            'files': [
                'pca/scripts/utils/analyze_full_vcf.py',
                'pca/scripts/utils/extract_real_vcf_sample.py',
                'pca/scripts/utils/full_genome_extractor.py',
                'pca/scripts/utils/improved_ancestry_loader.py',
                'pca/scripts/utils/indexed_multi_chr_extractor.py',
                'pca/scripts/utils/optimized_genome_extractor.py',
                'pca/scripts/utils/quick_multi_chr_extractor.py',
                'pca/scripts/utils/test_ancestry_loading.py',
            ]
        },
        
        'delete': {
            'description': 'Files to delete completely',
            'files': [
                'pipeline/core/__pycache__/',
                'pipeline/pipeline_io/__pycache__/',
                'pipeline/qc/__pycache__/',
                'pipeline/utils/__pycache__/',
                'code_review_analysis.json',  # Generated file
                'project_structure.json',     # Generated file
            ]
        }
    }
    
    return plan

def print_organization_plan():
    """Print the organization plan for review"""
    plan = create_organization_plan()
    
    print("Code Organization Plan for GitHub")
    print("=" * 50)
    
    total_keep = len(plan['keep_active']['files']) + len(plan['keep_supporting']['files'])
    total_archive = sum(len(category['files']) for key, category in plan.items() if key.startswith('archive'))
    total_delete = len(plan['delete']['files'])
    
    print(f"\n📊 SUMMARY:")
    print(f"Files to keep: {total_keep}")
    print(f"Files to archive: {total_archive}")
    print(f"Files to delete: {total_delete}")
    
    for category_key, category in plan.items():
        print(f"\n{category_key.upper().replace('_', ' ')} ({len(category['files'])} items):")
        print(f"Description: {category['description']}")
        
        for file_path in category['files']:
            if Path(file_path).exists():
                if Path(file_path).is_dir():
                    file_count = len(list(Path(file_path).rglob('*.py')))
                    print(f"  📁 {file_path} ({file_count} Python files)")
                else:
                    file_size = Path(file_path).stat().st_size
                    print(f"  📄 {file_path} ({file_size:,} bytes)")
            else:
                print(f"  ❓ {file_path} (not found)")
    
    return plan

def create_archive_structure():
    """Create archive directory structure"""
    print(f"\n🗂️  RECOMMENDED ARCHIVE STRUCTURE:")
    print(f"archive/")
    print(f"├── experimental/          # Development and experimental scripts")
    print(f"├── legacy_pipeline/       # Original variant calling pipeline") 
    print(f"├── analysis_scripts/      # Variant analysis scripts")
    print(f"├── utils/                 # Development utility scripts")
    print(f"└── README.md             # Archive documentation")
    
def main():
    """Main function"""
    plan = print_organization_plan()
    create_archive_structure()
    
    print(f"\n💡 RECOMMENDATIONS:")
    print(f"1. Create 'archive/' directory for historical code")
    print(f"2. Keep only the {len(plan['keep_active']['files']) + len(plan['keep_supporting']['files'])} essential files for GitHub")
    print(f"3. Archive {sum(len(cat['files']) for key, cat in plan.items() if key.startswith('archive'))} files")
    print(f"4. Delete {len(plan['delete']['files'])} cache/temporary files")
    print(f"\n✅ This will create a clean, focused repository for the ancestry inference pipeline")
    
    return plan

if __name__ == "__main__":
    main()