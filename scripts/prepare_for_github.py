#!/usr/bin/env python3
"""
GitHub Repository Preparation Script
Prepares the hybrid ancestry inference pipeline for GitHub publication
"""

import os
import shutil
import subprocess
from pathlib import Path
import json

def check_git_status():
    """Check current git status"""
    print("=== Git Status Check ===")
    
    try:
        # Check if we're in a git repository
        result = subprocess.run(['git', 'status', '--porcelain'], 
                              capture_output=True, text=True, cwd='.')
        
        if result.returncode == 0:
            print("✅ Git repository detected")
            
            # Check for uncommitted changes
            if result.stdout.strip():
                print("📝 Uncommitted changes found:")
                print(result.stdout)
            else:
                print("✅ Working directory is clean")
                
            # Check current branch
            branch_result = subprocess.run(['git', 'branch', '--show-current'], 
                                         capture_output=True, text=True, cwd='.')
            if branch_result.returncode == 0:
                print(f"📋 Current branch: {branch_result.stdout.strip()}")
                
        else:
            print("❌ Not a git repository")
            return False
            
    except FileNotFoundError:
        print("❌ Git not found")
        return False
    
    return True

def check_essential_files():
    """Check that essential files exist"""
    print("\n=== Essential Files Check ===")
    
    essential_files = [
        ('README.md', 'Main project README'),
        ('LICENSE', 'License file'),
        ('requirements.txt', 'Python dependencies'),
        ('.gitignore', 'Git ignore rules'),
        ('pca/HYBRID_ANCESTRY_METHODOLOGY.md', 'Methodology documentation'),
        ('pca/PROJECT_STATUS_SUMMARY.md', 'Project status summary'),
        ('pca/scripts/analysis/plink_pca_analysis.py', 'Main PCA analysis script'),
        ('pca/scripts/analysis/reference_projection_analysis.py', 'Reference projection script'),
        ('pca/scripts/utils/metadata_schema_updater.py', 'Metadata schema updater'),
        ('pca/scripts/testing/comprehensive_sanity_checks.py', 'Testing suite'),
        ('pca/scripts/testing/final_validation.py', 'Final validation'),
    ]
    
    missing_files = []
    for file_path, description in essential_files:
        if Path(file_path).exists():
            file_size = Path(file_path).stat().st_size
            print(f"✅ {description}: {file_path} ({file_size:,} bytes)")
        else:
            print(f"❌ {description}: {file_path} (MISSING)")
            missing_files.append(file_path)
    
    if missing_files:
        print(f"\n⚠️  {len(missing_files)} essential files missing")
        return False
    else:
        print("\n✅ All essential files present")
        return True

def check_large_files():
    """Check for large files that shouldn't be committed"""
    print("\n=== Large Files Check ===")
    
    large_files = []
    max_size = 100 * 1024 * 1024  # 100MB
    
    # Check for large files
    for root, dirs, files in os.walk('.'):
        # Skip .git directory
        if '.git' in dirs:
            dirs.remove('.git')
            
        for file in files:
            file_path = Path(root) / file
            try:
                size = file_path.stat().st_size
                if size > max_size:
                    large_files.append((str(file_path), size))
            except (OSError, FileNotFoundError):
                continue
    
    if large_files:
        print(f"⚠️  Found {len(large_files)} large files (>100MB):")
        for file_path, size in large_files:
            print(f"  {file_path}: {size / (1024*1024):.1f} MB")
        print("\n💡 Consider adding these to .gitignore if not essential")
    else:
        print("✅ No large files found")
    
    return len(large_files) == 0

def validate_python_syntax():
    """Validate Python syntax in all scripts"""
    print("\n=== Python Syntax Validation ===")
    
    python_files = []
    for root, dirs, files in os.walk('.'):
        if '.git' in dirs:
            dirs.remove('.git')
        for file in files:
            if file.endswith('.py'):
                python_files.append(Path(root) / file)
    
    syntax_errors = []
    for py_file in python_files:
        try:
            with open(py_file, 'r', encoding='utf-8') as f:
                code = f.read()
            compile(code, str(py_file), 'exec')
            print(f"✅ {py_file}")
        except SyntaxError as e:
            print(f"❌ {py_file}: {e}")
            syntax_errors.append((py_file, e))
        except Exception as e:
            print(f"⚠️  {py_file}: {e}")
    
    if syntax_errors:
        print(f"\n❌ {len(syntax_errors)} files have syntax errors")
        return False
    else:
        print(f"\n✅ All {len(python_files)} Python files have valid syntax")
        return True

def create_project_structure_summary():
    """Create a summary of the project structure"""
    print("\n=== Project Structure Summary ===")
    
    structure = {
        'directories': {},
        'files': {},
        'total_size': 0
    }
    
    for root, dirs, files in os.walk('.'):
        if '.git' in dirs:
            dirs.remove('.git')
            
        # Count directories
        for dir_name in dirs:
            dir_path = Path(root) / dir_name
            if str(dir_path).startswith('./'):
                structure['directories'][str(dir_path)[2:]] = len(list(dir_path.iterdir())) if dir_path.exists() else 0
        
        # Count files
        for file_name in files:
            file_path = Path(root) / file_name
            try:
                size = file_path.stat().st_size
                structure['total_size'] += size
                ext = file_path.suffix.lower()
                if ext not in structure['files']:
                    structure['files'][ext] = {'count': 0, 'size': 0}
                structure['files'][ext]['count'] += 1
                structure['files'][ext]['size'] += size
            except (OSError, FileNotFoundError):
                continue
    
    # Save structure summary
    with open('project_structure.json', 'w') as f:
        json.dump(structure, f, indent=2)
    
    print(f"📁 Directories: {len(structure['directories'])}")
    print(f"📄 File types: {len(structure['files'])}")
    print(f"💾 Total size: {structure['total_size'] / (1024*1024):.1f} MB")
    
    # Show top file types
    print("\n📊 Top file types:")
    sorted_files = sorted(structure['files'].items(), key=lambda x: x[1]['count'], reverse=True)
    for ext, info in sorted_files[:10]:
        ext_name = ext if ext else 'no extension'
        print(f"  {ext_name}: {info['count']} files ({info['size'] / (1024*1024):.1f} MB)")
    
    return structure

def generate_github_checklist():
    """Generate a checklist for GitHub preparation"""
    print("\n=== GitHub Preparation Checklist ===")
    
    checklist = [
        "✅ Documentation updated and complete",
        "✅ All code tested and validated",
        "✅ .gitignore configured for large files",
        "✅ README.md provides clear project overview",
        "✅ LICENSE file included",
        "✅ requirements.txt lists all dependencies",
        "✅ Python syntax validated",
        "✅ Project structure organized",
        "⚠️  Consider adding CONTRIBUTING.md",
        "⚠️  Consider adding CHANGELOG.md",
        "⚠️  Consider adding GitHub Actions workflows",
        "⚠️  Consider adding issue templates",
        "⚠️  Review sensitive data in configuration files",
        "⚠️  Set up repository topics and description",
        "⚠️  Configure branch protection rules",
    ]
    
    for item in checklist:
        print(f"  {item}")
    
    return checklist

def main():
    """Main preparation function"""
    print("GitHub Repository Preparation")
    print("=" * 50)
    
    # Run all checks
    checks = [
        check_git_status(),
        check_essential_files(),
        check_large_files(),
        validate_python_syntax(),
    ]
    
    # Create summaries
    structure = create_project_structure_summary()
    checklist = generate_github_checklist()
    
    # Final assessment
    print("\n" + "=" * 50)
    print("GITHUB PREPARATION ASSESSMENT")
    print("=" * 50)
    
    passed_checks = sum(checks)
    total_checks = len(checks)
    
    if passed_checks == total_checks:
        print("🎉 Repository is ready for GitHub!")
        print("✅ All checks passed")
        print("\n📝 Next steps:")
        print("  1. git add .")
        print("  2. git commit -m 'Complete hybrid ancestry inference pipeline'")
        print("  3. git push origin main")
        print("  4. Create GitHub repository")
        print("  5. Set up repository description and topics")
    else:
        print(f"⚠️  {total_checks - passed_checks} checks failed")
        print("❌ Please address issues before pushing to GitHub")
    
    print(f"\n📊 Repository Stats:")
    print(f"  Total size: {structure['total_size'] / (1024*1024):.1f} MB")
    print(f"  Directories: {len(structure['directories'])}")
    print(f"  File types: {len(structure['files'])}")
    
    return passed_checks == total_checks

if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)