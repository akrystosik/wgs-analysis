#!/usr/bin/env python3
"""
Script to lookup ENCODE ENCSR experiment accessions from ENCFF file accessions
Using ENCODE Portal API to complete the experiment mapping.
"""

import requests
import json
import time
from typing import Dict, List, Optional

class EncodePortalLookup:
    def __init__(self):
        self.base_url = "https://www.encodeproject.org"
        self.headers = {'Accept': 'application/json'}
        
    def lookup_file_experiment(self, encff_accession: str) -> Optional[str]:
        """
        Look up the ENCSR experiment accession from an ENCFF file accession.
        """
        try:
            url = f"{self.base_url}/{encff_accession}/?format=json"
            print(f"Querying: {url}")
            
            response = requests.get(url, headers=self.headers)
            response.raise_for_status()
            
            data = response.json()
            
            # The experiment field should contain the ENCSR accession
            if 'dataset' in data:
                return data['dataset'].split('/')[-2]  # Extract ENCSR from path
            elif 'experiment' in data:
                return data['experiment'].split('/')[-2]
            else:
                print(f"No experiment/dataset field found for {encff_accession}")
                return None
                
        except requests.exceptions.RequestException as e:
            print(f"Error querying {encff_accession}: {e}")
            return None
        except KeyError as e:
            print(f"Unexpected response format for {encff_accession}: {e}")
            return None

    def lookup_experiment_details(self, encsr_accession: str) -> Dict:
        """
        Get detailed information about an ENCSR experiment.
        """
        try:
            url = f"{self.base_url}/{encsr_accession}/?format=json"
            print(f"Querying experiment details: {url}")
            
            response = requests.get(url, headers=self.headers)
            response.raise_for_status()
            
            data = response.json()
            
            return {
                'accession': encsr_accession,
                'assay_title': data.get('assay_title', 'Unknown'),
                'biosample_term_name': data.get('biosample_term_name', 'Unknown'),
                'description': data.get('description', ''),
                'lab': data.get('lab', {}).get('title', 'Unknown'),
                'status': data.get('status', 'Unknown')
            }
            
        except requests.exceptions.RequestException as e:
            print(f"Error querying experiment {encsr_accession}: {e}")
            return {}

def main():
    """
    Main function to lookup all missing ENCSR mappings
    """
    lookup = EncodePortalLookup()
    
    # Cell lines and their FASTQ accessions that need ENCSR lookup
    cell_line_fastqs = {
        'A549': ['ENCFF122NPY', 'ENCFF846WHK'],
        'NCI-H460': ['ENCFF022XPK', 'ENCFF534EUU'],
        'Panc1': ['ENCFF477JTA', 'ENCFF896PZG'],
        'HepG2': ['ENCFF320KMG', 'ENCFF045JFV'],
        'sknmc': ['ENCFF212SDN', 'ENCFF675TKC']
    }
    
    results = {}
    
    for cell_line, fastq_files in cell_line_fastqs.items():
        print(f"\n=== Processing {cell_line} ===")
        cell_results = {'cell_line': cell_line, 'fastq_files': fastq_files, 'encsr_accessions': set()}
        
        for fastq in fastq_files:
            print(f"Looking up {fastq}...")
            encsr = lookup.lookup_file_experiment(fastq)
            if encsr:
                cell_results['encsr_accessions'].add(encsr)
                print(f"  Found ENCSR: {encsr}")
            else:
                print(f"  No ENCSR found for {fastq}")
            
            time.sleep(0.5)  # Be nice to the API
        
        # Convert set to list for JSON serialization
        cell_results['encsr_accessions'] = list(cell_results['encsr_accessions'])
        results[cell_line] = cell_results
        
        # Get experiment details
        for encsr in cell_results['encsr_accessions']:
            details = lookup.lookup_experiment_details(encsr)
            if details:
                cell_results[f'{encsr}_details'] = details
            time.sleep(0.5)
    
    # Save results
    output_file = 'encode_encsr_lookup_results.json'
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"\n=== Results Summary ===")
    print(f"Results saved to: {output_file}")
    
    for cell_line, data in results.items():
        encsr_list = data['encsr_accessions']
        if encsr_list:
            print(f"{cell_line}: {', '.join(encsr_list)}")
        else:
            print(f"{cell_line}: No ENCSR found")
    
    return results

if __name__ == "__main__":
    main()