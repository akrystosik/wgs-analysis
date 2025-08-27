#!/usr/bin/env python3
"""
Code Review Analysis Script
Analyzes all Python files to determine their usage status and relevance
"""

import os
import re
import ast
from pathlib import Path
from collections import defaultdict
import json

def analyze_file_usage(file_path):
    """Analyze a Python file to determine its usage status"""
    try:
        with open(file_path, 'r', encoding='utf-8') as f:
            content = f.read()
        
        # Check file size
        file_size = len(content)
        
        # Check for main execution
        has_main = '__main__' in content
        
        # Check for class definitions
        class_count = len(re.findall(r'class\s+\w+', content))
        
        # Check for function definitions
        func_count = len(re.findall(r'def\s+\w+', content))
        
        # Check for imports
        import_count = len(re.findall(r'(?:from\s+\w+\s+)?import\s+', content))
        
        # Check for docstrings
        has_docstring = '"""' in content or "'''" in content
        
        # Check for recent modifications (based on content analysis)
        has_recent_dates = any(year in content for year in ['2025', '2024'])
        
        # Check for TODO/FIXME/NOTE comments
        todo_count = len(re.findall(r'#\s*(?:TODO|FIXME|NOTE|HACK)', content, re.IGNORECASE))
        
        # Check for test functions
        test_count = len(re.findall(r'def\s+test_', content))
        
        # Check for specific patterns that indicate usage
        patterns = {
            'pca_analysis': 'pca.*analysis',
            'ancestry': 'ancestry',
            'reference_projection': 'reference.*projection',
            'metadata': 'metadata',
            'validation': 'validation|test',
            'utils': 'util',
            'experimental': 'experimental',
            'deprecated': 'deprecated|old|legacy',
            'pipeline': 'pipeline',
            'variant_calling': 'variant.*call'
        }
        
        pattern_matches = {}
        for pattern_name, pattern in patterns.items():
            pattern_matches[pattern_name] = len(re.findall(pattern, content, re.IGNORECASE))
        
        return {
            'file_size': file_size,
            'has_main': has_main,
            'class_count': class_count,
            'func_count': func_count,
            'import_count': import_count,
            'has_docstring': has_docstring,
            'has_recent_dates': has_recent_dates,
            'todo_count': todo_count,
            'test_count': test_count,
            'pattern_matches': pattern_matches
        }
        
    except Exception as e:
        return {'error': str(e)}

def categorize_files():
    """Categorize Python files by their purpose and usage"""
    
    # Get all Python files
    python_files = []
    for root, dirs, files in os.walk('.'):
        if '.git' in dirs:
            dirs.remove('.git')
        for file in files:
            if file.endswith('.py'):
                python_files.append(Path(root) / file)
    
    # Analyze each file
    file_analysis = {}
    for file_path in python_files:
        analysis = analyze_file_usage(file_path)
        file_analysis[str(file_path)] = analysis
    
    # Categorize files
    categories = {
        'active_core': [],          # Currently used core functionality
        'active_analysis': [],      # Currently used analysis scripts
        'active_utils': [],         # Currently used utilities
        'active_testing': [],       # Currently used testing scripts
        'experimental': [],         # Experimental/development scripts
        'legacy_pipeline': [],      # Legacy pipeline components
        'legacy_analysis': [],      # Legacy analysis scripts
        'deprecated': [],           # Deprecated/unused scripts
        'cache_files': [],          # __pycache__ files
        'small_utility': []         # Small utility scripts
    }
    
    for file_path, analysis in file_analysis.items():
        if 'error' in analysis:
            continue
            
        file_path_obj = Path(file_path)
        
        # Cache files
        if '__pycache__' in str(file_path):
            categories['cache_files'].append(file_path)
            continue
        
        # Recently created/modified files (our new analysis pipeline)
        if analysis['has_recent_dates'] and any(keyword in file_path.lower() for keyword in ['plink', 'reference', 'projection', 'metadata', 'testing', 'validation']):
            if 'testing' in file_path or 'validation' in file_path:
                categories['active_testing'].append(file_path)
            elif 'analysis' in file_path:
                categories['active_analysis'].append(file_path)
            elif 'utils' in file_path:
                categories['active_utils'].append(file_path)
            else:
                categories['active_core'].append(file_path)
            continue
        
        # Experimental scripts
        if 'experimental' in file_path or 'quick' in file_path or 'simple' in file_path:
            categories['experimental'].append(file_path)
            continue
        
        # Legacy pipeline components
        if 'pipeline' in file_path and not analysis['has_recent_dates']:
            categories['legacy_pipeline'].append(file_path)
            continue
        
        # Legacy analysis scripts
        if any(keyword in file_path for keyword in ['scripts/', 'core/', 'enhanced_pca']) and not analysis['has_recent_dates']:
            categories['legacy_analysis'].append(file_path)
            continue
        
        # Small utility files
        if analysis['file_size'] < 1000 and analysis['func_count'] < 3:
            categories['small_utility'].append(file_path)
            continue
        
        # Everything else - analyze by content
        if analysis['pattern_matches']['pca_analysis'] > 0:
            categories['active_analysis'].append(file_path)
        elif analysis['pattern_matches']['validation'] > 0:
            categories['active_testing'].append(file_path)
        elif analysis['pattern_matches']['utils'] > 0:
            categories['active_utils'].append(file_path)
        elif analysis['pattern_matches']['deprecated'] > 0:
            categories['deprecated'].append(file_path)
        else:
            categories['active_core'].append(file_path)
    
    return categories, file_analysis

def generate_recommendations(categories, file_analysis):
    """Generate recommendations for file organization"""
    
    recommendations = {
        'keep_active': [],
        'archive_experimental': [],
        'archive_legacy': [],
        'delete_cache': [],
        'review_needed': []
    }
    
    # Files to definitely keep (active core functionality)
    for category in ['active_core', 'active_analysis', 'active_utils', 'active_testing']:
        recommendations['keep_active'].extend(categories[category])
    
    # Files to archive (experimental/legacy but might be useful)
    for category in ['experimental', 'legacy_analysis']:
        recommendations['archive_experimental'].extend(categories[category])
    
    # Files to archive (legacy pipeline)
    recommendations['archive_legacy'].extend(categories['legacy_pipeline'])
    
    # Files to delete (cache)
    recommendations['delete_cache'].extend(categories['cache_files'])
    
    # Files that need review
    for category in ['deprecated', 'small_utility']:
        recommendations['review_needed'].extend(categories[category])
    
    return recommendations

def main():
    """Main analysis function"""
    print("Code Review Analysis")
    print("=" * 50)
    
    # Analyze all files
    categories, file_analysis = categorize_files()
    
    # Generate recommendations
    recommendations = generate_recommendations(categories, file_analysis)
    
    # Print analysis results
    print("\n=== FILE CATEGORIZATION ===")
    for category, files in categories.items():
        if files:
            print(f"\n{category.upper()} ({len(files)} files):")
            for file in sorted(files):
                size = file_analysis.get(file, {}).get('file_size', 0)
                print(f"  {file} ({size:,} bytes)")
    
    print("\n=== RECOMMENDATIONS ===")
    
    print(f"\n✅ KEEP ACTIVE ({len(recommendations['keep_active'])} files):")
    for file in sorted(recommendations['keep_active']):
        print(f"  {file}")
    
    print(f"\n📦 ARCHIVE EXPERIMENTAL ({len(recommendations['archive_experimental'])} files):")
    for file in sorted(recommendations['archive_experimental']):
        print(f"  {file}")
    
    print(f"\n📦 ARCHIVE LEGACY ({len(recommendations['archive_legacy'])} files):")
    for file in sorted(recommendations['archive_legacy']):
        print(f"  {file}")
    
    print(f"\n🗑️  DELETE CACHE ({len(recommendations['delete_cache'])} files):")
    for file in sorted(recommendations['delete_cache']):
        print(f"  {file}")
    
    print(f"\n🔍 REVIEW NEEDED ({len(recommendations['review_needed'])} files):")
    for file in sorted(recommendations['review_needed']):
        print(f"  {file}")
    
    # Summary statistics
    total_files = sum(len(files) for files in categories.values())
    active_files = len(recommendations['keep_active'])
    
    print(f"\n=== SUMMARY ===")
    print(f"Total Python files: {total_files}")
    print(f"Active files to keep: {active_files}")
    print(f"Files to archive: {len(recommendations['archive_experimental']) + len(recommendations['archive_legacy'])}")
    print(f"Files to delete: {len(recommendations['delete_cache'])}")
    print(f"Files needing review: {len(recommendations['review_needed'])}")
    
    # Save detailed analysis
    with open('code_review_analysis.json', 'w') as f:
        json.dump({
            'categories': categories,
            'recommendations': recommendations,
            'file_analysis': file_analysis
        }, f, indent=2)
    
    print(f"\n📄 Detailed analysis saved to: code_review_analysis.json")
    
    return recommendations

if __name__ == "__main__":
    main()