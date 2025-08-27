# encode_analysis/pipeline/utils/naming.py

from typing import Dict

# encode_analysis/pipeline/utils/naming.py

def get_standardized_names(pair_id: str) -> Dict[str, str]:
    """
    Creates standardized file names matching our established directory structure.
    
    Given a pair_id like 'A549', this generates paths and names that match:
    /data/bam/A549/A549_WGS_pair1.marked_duplicates.bam
    
    The function extracts the base sample name (everything before _WGS_pair1)
    and uses it consistently for directory structure and file naming.
    
    Args:
        pair_id: Base sample identifier (e.g., 'A549')
        
    Returns:
        Dictionary containing standardized names for all pipeline stages
    """
    # Extract base name (remove _WGS_pair1 if present)
    base_name = pair_id.replace('_WGS_pair1', '')
    
    return {
        # Directory and identification fields
        'sample_dir': base_name,
        'sample_name': base_name,
        
        # BAM processing files
        'sorted_bam': f"{base_name}_WGS_pair1.sorted.bam",
        'marked_bam': f"{base_name}_WGS_pair1.marked_duplicates.bam",
        'metrics_file': f"{base_name}_WGS_pair1.marked_dup_metrics.txt",
        
        # Variant calling outputs
        'vcf_output': f"{base_name}_WGS_pair1.deepvariant.vcf.gz",
        'gvcf_output': f"{base_name}_WGS_pair1.deepvariant.g.vcf.gz"
    }
    
def validate_pipeline_names(pair_id: str) -> None:
    """
    Validate that the provided pair_id follows pipeline naming conventions.
    
    This function checks:
    1. Pair ID format and characters
    2. Presence/absence of _WGS_pair1 suffix
    3. Basic sample name requirements
    
    Args:
        pair_id: The sample identifier to validate
        
    Raises:
        ValueError: If the pair_id doesn't meet naming requirements
    """
    # Check for invalid characters
    invalid_chars = set('<>:"/\\|?*')
    if any(char in pair_id for char in invalid_chars):
        raise ValueError(f"pair_id contains invalid characters: {pair_id}")
        
    # If it has _WGS_pair1, validate base name
    if pair_id.endswith('_WGS_pair1'):
        base_name = pair_id.replace('_WGS_pair1', '')
        if not base_name:
            raise ValueError("pair_id has _WGS_pair1 suffix but no base name")
    else:
        # If it doesn't have suffix, validate it could form a valid full name
        if not pair_id:
            raise ValueError("pair_id cannot be empty")
            
