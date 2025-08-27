import requests
import logging

logger = logging.getLogger(__name__)

def get_encode_metadata(file_id):
    """Get metadata from ENCODE API for a file."""
    url = f"https://www.encodeproject.org/files/{file_id}/?format=json"
    headers = {'Accept': 'application/json'}
    try:
        response = requests.get(url, headers=headers)
        if response.ok:
            return response.json()
        else:
            logger.error(f"Failed to get metadata for {file_id}: {response.status_code}")
            return None
    except Exception as e:
        logger.error(f"Error fetching metadata for {file_id}: {str(e)}")
        return None

def get_file_id_from_path(filepath):
    """Extract ENCODE file ID from filepath or URL."""
    filepath = str(filepath)
    filename = filepath.split('/')[-1]
    return filename.split('.')[0]

def extract_file_id_from_paired_with(paired_with):
    """Extract the file ID from the 'paired_with' field (handles full paths)."""
    if paired_with:
        # Sometimes paired_with is a full path like /files/ENCFFXXXXXX/, extract the ID
        return paired_with.strip('/').split('/')[-1]
    return None

def verify_fastq_pairs(fastq_paths):
    """
    Verify and organize FASTQ pairs using ENCODE metadata.
    
    Args:
        fastq_paths: List of FASTQ file paths
        
    Returns:
        dict: {'read1': path_to_read1, 'read2': path_to_read2}
    """
    if len(fastq_paths) != 2:
        raise ValueError(f"Expected exactly 2 FASTQ files, got {len(fastq_paths)}")

    # Get file IDs from paths
    file_ids = [get_file_id_from_path(path) for path in fastq_paths]
    logger.info(f"Verifying FASTQ pair: {file_ids[0]} - {file_ids[1]}")
    
    # Get metadata from ENCODE API
    metadata1 = get_encode_metadata(file_ids[0])
    metadata2 = get_encode_metadata(file_ids[1])
    
    if not metadata1 or not metadata2:
        raise ValueError("Failed to retrieve ENCODE metadata for one or both files.")
    
    # Verify they're from the same experiment
    exp1 = metadata1.get('dataset', '').split('/')[-2]
    exp2 = metadata2.get('dataset', '').split('/')[-2]
    logger.info(f"File 1 from experiment: {exp1}, File 2 from experiment: {exp2}")
    
    if exp1 != exp2:
        raise ValueError(f"Files are from different experiments: {exp1} vs {exp2}")

    # Verify they're properly paired using the paired_with field
    paired1 = extract_file_id_from_paired_with(metadata1.get('paired_with', None))
    paired2 = extract_file_id_from_paired_with(metadata2.get('paired_with', None))
    
    if not paired1 or not paired2:
        raise ValueError("One or both files are missing the 'paired_with' field in the metadata.")
    
    logger.info(f"File 1 is paired with: {paired1}, File 2 is paired with: {paired2}")
    
    if paired1 != file_ids[1] or paired2 != file_ids[0]:
        raise ValueError("Files are not properly paired in ENCODE metadata")
    
    # Determine which file is read1 and read2 based on the paired_end field
    pe1 = str(metadata1.get('paired_end', ''))
    pe2 = str(metadata2.get('paired_end', ''))
    
    if pe1 == '1':
        logger.info(f"Determined read1: {fastq_paths[0]}, read2: {fastq_paths[1]}")
        return {'read1': fastq_paths[0], 'read2': fastq_paths[1]}
    elif pe2 == '1':
        logger.info(f"Determined read1: {fastq_paths[1]}, read2: {fastq_paths[0]}")
        return {'read1': fastq_paths[1], 'read2': fastq_paths[0]}
    else:
        raise ValueError("Could not determine read1/read2 from metadata")
