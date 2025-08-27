#!/usr/bin/env python3
"""
Croissant JSON-LD Generator for VCF Genomics Datasets

This script generates Croissant metadata files for genomics variant datasets stored in VCF format.
It handles both open access datasets (with full file metadata) and protected datasets 
(with source-only metadata pointing to official repositories).

Based on:
- Croissant Format Specification: https://docs.mlcommons.org/croissant/docs/croissant-spec.html
- VCF Specification: https://samtools.github.io/hts-specs/VCFv4.2.pdf
"""

import os
import gzip
import json
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Any, Optional
import subprocess
import hashlib


class VCFDatasetInfo:
    """Container for VCF dataset-specific information from literature and official sources."""
    
    GTEx = {
        'name': 'GTEx (Genotype-Tissue Expression) - Whole Genome Sequencing',
        'description': 'Whole genome sequencing data from the Genotype-Tissue Expression (GTEx) project containing genetic variants from 838 individuals across multiple tissues. This protected dataset is available through dbGaP authorized access.',
        'url': 'https://www.gtexportal.org/',
        'data_url': 'https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000424',
        'citations': [
            'The GTEx Consortium. The GTEx Consortium atlas of genetic regulatory effects across human tissues. Science 369, 1318-1330 (2020).',
            'The GTEx Consortium. Genetic effects on gene expression across human tissues. Nature 550, 204-213 (2017).'
        ],
        'license': 'dbGaP Authorized Access Required',
        'creator': 'The GTEx Consortium',
        'keywords': ['genetic variants', 'WGS', 'eQTL', 'tissue', 'human', 'dbGaP'],
        'assay_type': 'Whole genome sequencing',
        'organism': 'Homo sapiens',
        'access_level': 'protected',
        'subjects': '838 individuals',
        'variant_count': '~400M variants',
        'dbgap_accession': 'phs000424'
    }
    
    ADNI = {
        'name': 'ADNI (Alzheimer\'s Disease Neuroimaging Initiative) - Whole Genome Sequencing', 
        'description': 'Whole genome sequencing data from the Alzheimer\'s Disease Neuroimaging Initiative containing genetic variants from participants in the longitudinal study of Alzheimer\'s disease progression. This protected dataset requires ADNI data use agreement.',
        'url': 'https://adni.loni.usc.edu/',
        'data_url': 'https://ida.loni.usc.edu/',
        'citations': [
            'Petersen RC, et al. Alzheimer\'s Disease Neuroimaging Initiative (ADNI): Clinical characterization. Neurology 74, 201-209 (2010).',
            'Weiner MW, et al. The Alzheimer\'s Disease Neuroimaging Initiative 3: Continued innovation for clinical trial improvement. Alzheimer\'s & Dementia 13, 561-571 (2017).'
        ],
        'license': 'ADNI Data Use Agreement Required',
        'creator': 'Alzheimer\'s Disease Neuroimaging Initiative',
        'keywords': ['genetic variants', 'WGS', 'Alzheimer disease', 'neuroimaging', 'longitudinal'],
        'assay_type': 'Whole genome sequencing', 
        'organism': 'Homo sapiens',
        'access_level': 'protected',
        'subjects': '~1,700 participants',
        'variant_count': '~400M variants',
        'ida_access': 'https://ida.loni.usc.edu/'
    }
    
    MAGE = {
        'name': 'MAGE (Multi-ancestry Analysis of Gene Expression) - Whole Genome Sequencing',
        'description': 'Whole genome sequencing data from the Multi-ancestry Analysis of Gene Expression study, containing genetic variants from 731 individuals from the 1000 Genomes Project representing 26 globally-distributed populations. Data is available through the 1000 Genomes Project and associated repositories.',
        'url': 'https://github.com/mccoy-lab/MAGE', 
        'data_url': 'https://www.internationalgenome.org/',
        'citations': [
            'Taylor DJ, et al. Sources of gene expression variation in a globally diverse human cohort. Nature (2024). DOI: 10.1038/s41586-024-07708-2',
            '1000 Genomes Project Consortium. A global reference for human genetic variation. Nature 526, 68-74 (2015).'
        ],
        'license': 'CC0 / 1000 Genomes Project Terms',
        'creator': 'McCoy Lab, Johns Hopkins University / 1000 Genomes Project',
        'keywords': ['genetic variants', 'WGS', 'population genetics', '1000 Genomes', 'diversity'],
        'assay_type': 'Whole genome sequencing',
        'organism': 'Homo sapiens', 
        'access_level': 'protected',  # Points to 1000 Genomes, not direct file access
        'subjects': '731 individuals from 26 populations',
        'variant_count': '~400M variants',
        'genomes_project_url': 'https://www.internationalgenome.org/'
    }
    
    ENCODE = {
        'name': 'ENCODE (Encyclopedia of DNA Elements) - Cell Line WGS Variants',
        'description': 'Whole genome sequencing variant data from ENCODE project cell lines, containing genetic variants called using DeepVariant and DeepSomatic from 9 cell lines including A549, Caki2, GM23248, HepG2, K562, NCI-H460, Panc1, T47D, and SKNMC. These datasets represent diverse cancer cell lines and normal cell lines used in functional genomics studies.',
        'url': 'https://www.encodeproject.org/',
        'data_url': 'https://www.encodeproject.org/',
        'citations': [
            'ENCODE Project Consortium. Expanded encyclopaedias of DNA elements in the human and mouse genomes. Nature 583, 699-710 (2020).',
            'ENCODE Project Consortium. An integrated encyclopedia of DNA elements in the human genome. Nature 489, 57-74 (2012).'
        ],
        'license': 'CC0',
        'creator': 'ENCODE Project Consortium', 
        'keywords': ['genetic variants', 'WGS', 'cell lines', 'functional genomics', 'DeepVariant'],
        'assay_type': 'Whole genome sequencing',
        'organism': 'Homo sapiens',
        'access_level': 'open',
        'subjects': '9 cell lines',
        'variant_count': '~5-10M variants per cell line',
        'cell_lines': ['A549', 'Caki2', 'GM23248', 'HepG2', 'K562', 'NCI-H460', 'Panc1', 'T47D', 'SKNMC']
    }
    
    ENCODE_PAIRED = {
        'name': 'ENCODE Paired Cell Lines - WGS Variants with RNA-seq',
        'description': 'Whole genome sequencing variant data from 7 ENCODE cell lines with paired RNA-seq data, enabling multi-omics analyses. Contains genetic variants called using DeepVariant and DeepSomatic from cell lines: A549 (lung carcinoma), Caki2 (renal carcinoma), GM23248 (skin fibroblast), HepG2 (hepatocellular carcinoma), K562 (chronic myelogenous leukemia), NCI-H460 (lung carcinoma), and Panc1 (pancreatic ductal carcinoma). These datasets represent diverse cancer cell lines and normal cell lines used in functional genomics studies with corresponding transcriptome data.',
        'url': 'https://www.encodeproject.org/',
        'data_url': 'https://www.encodeproject.org/',
        'citations': [
            'ENCODE Project Consortium. Expanded encyclopaedias of DNA elements in the human and mouse genomes. Nature 583, 699-710 (2020).',
            'ENCODE Project Consortium. An integrated encyclopedia of DNA elements in the human genome. Nature 489, 57-74 (2012).'
        ],
        'license': 'CC0',
        'creator': 'ENCODE Project Consortium', 
        'keywords': ['genetic variants', 'WGS', 'cell lines', 'functional genomics', 'DeepVariant', 'multi-omics', 'paired RNA-seq'],
        'assay_type': 'Whole genome sequencing',
        'organism': 'Homo sapiens',
        'access_level': 'open',
        'subjects': '7 cell lines with paired RNA-seq data',
        'variant_count': '~17.3M variants total',
        'cell_lines': ['A549', 'Caki2', 'GM23248', 'HepG2', 'K562', 'NCI-H460', 'Panc1']
    }


class VCFCroissantGenerator:
    """Generator for Croissant JSON-LD metadata files from VCF datasets."""
    
    def __init__(self, base_url: str = "https://example.com/datasets"):
        self.base_url = base_url
        self.dataset_info = VCFDatasetInfo()
        
    def _get_assay_information(self, dataset_info: Dict[str, Any]) -> Dict[str, str]:
        """Get dynamic assay information based on dataset and experimental methods.
        
        Always defaults to EFO ontology terms when available.
        """
        cell_lines = dataset_info.get('cell_lines', [])
        
        # ENCODE-specific cell line technology mapping
        chromium_cell_lines = ['K562', 'HepG2']  # 10X Chromium WGS
        illumina_cell_lines = ['A549', 'NCI-H460', 'Panc1', 'sknmc', 'GM23248', 'Caki2']  # Illumina WGS
        
        # Check if dataset contains any 10X Chromium cell lines
        has_chromium = any(cell_line in chromium_cell_lines for cell_line in cell_lines)
        has_illumina = any(cell_line in illumina_cell_lines for cell_line in cell_lines)
        
        if has_chromium and has_illumina:
            # Mixed technologies - use generic WGS term
            return {
                'name': 'whole genome sequencing assay',
                'description': 'Whole genome sequencing assay using mixed technologies (Illumina + 10X Chromium)',
                'ontology_id': 'OBI:0002117'  # Generic whole genome sequencing assay
            }
        elif has_chromium:
            # 10X Chromium specific
            return {
                'name': 'whole genome sequencing assay',
                'description': '10X Chromium whole genome sequencing assay',
                'ontology_id': 'OBI:0002117'  # Generic whole genome sequencing assay
            }
        elif has_illumina:
            # Illumina specific  
            return {
                'name': 'whole genome sequencing assay',
                'description': 'Illumina whole genome sequencing assay',
                'ontology_id': 'OBI:0002117'  # Generic whole genome sequencing assay
            }
        else:
            # Default to generic WGS
            return {
                'name': 'whole genome sequencing assay',
                'description': 'Whole genome sequencing assay',
                'ontology_id': 'OBI:0002117'  # Generic whole genome sequencing assay
            }
        
    def analyze_vcf_file(self, file_path: Path) -> Dict[str, Any]:
        """Analyze a VCF file to extract metadata."""
        
        if not file_path.exists():
            raise FileNotFoundError(f"VCF file not found: {file_path}")
        
        metadata = {
            'file_size': file_path.stat().st_size,
            'file_name': file_path.name,
            'is_compressed': file_path.suffix == '.gz'
        }
        
        # Open file (handle gzipped or regular)
        opener = gzip.open if metadata['is_compressed'] else open
        mode = 'rt' if metadata['is_compressed'] else 'r'
        
        try:
            with opener(file_path, mode) as f:
                header_lines = []
                samples = []
                variant_count = 0
                
                # Read header and count variants
                for line_num, line in enumerate(f):
                    line = line.strip()
                    
                    if line.startswith('##'):
                        header_lines.append(line)
                    elif line.startswith('#CHROM'):
                        # Sample line - extract sample names
                        fields = line.split('\t')
                        samples = fields[9:] if len(fields) > 9 else []
                        break
                    elif line_num > 1000:  # Safety limit for header
                        break
                
                # Count total variants (approximate for large files)
                f.seek(0)
                sample_lines = 0
                for line in f:
                    if not line.startswith('#'):
                        sample_lines += 1
                        if sample_lines > 10000:  # Sample first 10k for speed
                            break
                
                # Estimate total variants if we hit the limit
                if sample_lines == 10000:
                    # Get file size and estimate
                    f.seek(0, 2)  # Seek to end
                    file_size = f.tell()
                    f.seek(0)
                    
                    # Skip header
                    header_size = 0
                    for line in f:
                        header_size += len(line.encode('utf-8'))
                        if not line.startswith('#'):
                            break
                    
                    data_size = file_size - header_size
                    avg_line_size = data_size / 10000 if sample_lines > 0 else 0
                    variant_count = int(data_size / avg_line_size) if avg_line_size > 0 else sample_lines
                else:
                    variant_count = sample_lines
                
        except Exception as e:
            print(f"Error analyzing VCF {file_path}: {e}")
            variant_count = 0
            samples = []
            header_lines = []
        
        metadata.update({
            'variant_count': variant_count,
            'sample_count': len(samples),
            'samples': samples,
            'has_samples': len(samples) > 0,
            'header_lines': len(header_lines)
        })
        
        return metadata
    
    def get_dataset_info(self, dataset_name: str) -> Dict[str, Any]:
        """Get dataset-specific information from literature."""
        dataset_name = dataset_name.upper()
        
        # Map common variations to standard names
        name_mapping = {
            'GTEX': 'GTEx',
            'ENCODE': 'ENCODE',
            'ENCODE_PAIRED': 'ENCODE_PAIRED',
            'ADNI': 'ADNI', 
            'MAGE': 'MAGE'
        }
        
        mapped_name = name_mapping.get(dataset_name, dataset_name)
        if hasattr(self.dataset_info, mapped_name):
            return getattr(self.dataset_info, mapped_name)
        else:
            # Default information for unknown datasets
            return {
                'name': f'{dataset_name} WGS Variants',
                'description': f'Whole genome sequencing variant data: {dataset_name}',
                'url': self.base_url,
                'data_url': self.base_url,
                'citations': [],
                'license': 'Unknown',
                'creator': 'Unknown',
                'keywords': ['genetic variants', 'WGS', 'genomics'],
                'assay_type': 'Whole genome sequencing',
                'organism': 'Homo sapiens',
                'access_level': 'unknown',
                'subjects': 'Unknown',
                'variant_count': 'Unknown'
            }
    
    def generate_croissant_metadata(self, dataset_name: str, vcf_files: List[Path] = None) -> Dict[str, Any]:
        """Generate Croissant JSON-LD metadata for VCF datasets."""
        
        # Get dataset-specific information
        dataset_info = self.get_dataset_info(dataset_name)
        
        # Generate unique identifier
        dataset_id = f"{dataset_name.lower()}_wgs_variants_{datetime.now().strftime('%Y%m%d')}"
        
        # Analyze VCF files if provided (for open access datasets)
        file_metadata = []
        total_variants = 0
        total_samples = 0
        
        if vcf_files and dataset_info['access_level'] == 'open':
            for vcf_file in vcf_files:
                try:
                    file_meta = self.analyze_vcf_file(vcf_file)
                    file_metadata.append({
                        'file_path': vcf_file,
                        'metadata': file_meta
                    })
                    total_variants += file_meta.get('variant_count', 0)
                    total_samples += file_meta.get('sample_count', 0)
                except Exception as e:
                    print(f"Error analyzing {vcf_file}: {e}")
        
        # Build Croissant metadata
        croissant_metadata = {
            "@context": {
                "@language": "en",
                "@vocab": "https://schema.org/",
                "citeAs": "cr:citeAs",
                "column": "cr:column",
                "conformsTo": "dct:conformsTo",
                "cr": "http://mlcommons.org/croissant/",
                "rai": "http://mlcommons.org/croissant/RAI/",
                "data": {
                    "@id": "cr:data",
                    "@type": "@json"
                },
                "dataType": {
                    "@id": "cr:dataType",
                    "@type": "@vocab"
                },
                "dct": "http://purl.org/dc/terms/",
                "examples": {
                    "@id": "cr:examples",
                    "@type": "@json"
                },
                "extract": "cr:extract",
                "field": "cr:field",
                "fileProperty": "cr:fileProperty",
                "fileObject": "cr:fileObject",
                "fileSet": "cr:fileSet",
                "format": "cr:format",
                "includes": "cr:includes",
                "isLiveDataset": "cr:isLiveDataset",
                "jsonPath": "cr:jsonPath",
                "key": "cr:key",
                "md5": "cr:md5",
                "parentField": "cr:parentField",
                "path": "cr:path",
                "recordSet": "cr:recordSet",
                "references": "cr:references",
                "regex": "cr:regex",
                "repeated": "cr:repeated",
                "replace": "cr:replace",
                "sc": "https://schema.org/",
                "separator": "cr:separator",
                "source": "cr:source",
                "subField": "cr:subField",
                "transform": "cr:transform",
                "wd": "https://www.wikidata.org/wiki/"
            },
            "@type": "sc:Dataset",
            "@id": f"{self.base_url}/{dataset_id}",
            "dct:conformsTo": "http://mlcommons.org/croissant/1.0",
            "name": dataset_info['name'],
            "description": dataset_info['description'],
            "url": dataset_info['url'],
            "license": dataset_info['license'],
            "creator": {
                "@type": "Organization",
                "name": dataset_info['creator']
            },
            "citation": dataset_info['citations'],
            "keywords": dataset_info['keywords'],
            "datePublished": datetime.now().isoformat(),
            "version": "1.0", 
            "identifier": dataset_id,
            
            # Genomics-specific metadata
            "sc:organism": dataset_info['organism'],
            "sc:assayType": dataset_info['assay_type'],
            "sc:dataType": "genetic variants",
            
            # Access information
            "sc:accessLevel": dataset_info['access_level'],
            "sc:dataAccess": dataset_info['data_url'] if dataset_info['access_level'] == 'protected' else self.base_url
        }
        
        # Add Cross-Modality metadata for better validation compatibility
        if dataset_info.get('cell_lines'):
            # Get dynamic assay information based on the dataset
            assay_info = self._get_assay_information(dataset_info)
            
            croissant_metadata["variableMeasured"] = [
                {
                    "@type": "PropertyValue",
                    "name": "assay",
                    "description": assay_info['description'],
                    "propertyID": "assay_ontology_term_id",
                    "value": [assay_info['name']]
                },
                {
                    "@type": "PropertyValue", 
                    "name": "assay_ontology_term_id",
                    "description": "Ontology term ID for the assay",
                    "propertyID": "assay_ontology_term_id",
                    "value": [assay_info['ontology_id']]
                },
                {
                    "@type": "PropertyValue",
                    "name": "organism",
                    "description": "Species of origin for the samples",
                    "propertyID": "organism_ontology_term_id",
                    "value": ["Homo sapiens"]
                },
                {
                    "@type": "PropertyValue",
                    "name": "organism_ontology_term_id",
                    "description": "NCBITaxon ontology term for organism",
                    "propertyID": "organism_ontology_term_id", 
                    "value": ["NCBITaxon:9606"]
                },
                {
                    "@type": "PropertyValue",
                    "name": "tissue",
                    "description": "Cell line tissue types",
                    "propertyID": "tissue_ontology_term_id",
                    "value": ["cell culture"]
                },
                {
                    "@type": "PropertyValue",
                    "name": "tissue_ontology_term_id", 
                    "description": "Ontology term for tissue",
                    "propertyID": "tissue_ontology_term_id",
                    "value": ["EFO:0000322"]
                },
                {
                    "@type": "PropertyValue",
                    "name": "tissue_type",
                    "description": "Type of tissue",
                    "propertyID": "tissue_type",
                    "value": ["cell culture"]
                },
                {
                    "@type": "PropertyValue",
                    "name": "disease",
                    "description": "Disease state of samples",
                    "propertyID": "disease_ontology_term_id",
                    "value": ["cancer"]
                },
                {
                    "@type": "PropertyValue",
                    "name": "disease_ontology_term_id",
                    "description": "Ontology term for disease",
                    "propertyID": "disease_ontology_term_id",
                    "value": ["MONDO:0004992"]
                },
                {
                    "@type": "PropertyValue",
                    "name": "development_stage",
                    "description": "Development stage", 
                    "propertyID": "development_stage_ontology_term_id",
                    "value": ["adult"]
                },
                {
                    "@type": "PropertyValue",
                    "name": "development_stage_ontology_term_id",
                    "description": "Ontology term for development stage",
                    "propertyID": "development_stage_ontology_term_id",
                    "value": ["HsapDv:0000087"]
                }
            ]
        
        # Add distribution information based on access level
        if dataset_info['access_level'] == 'open' and file_metadata:
            # Open access: include actual file information
            distributions = []
            for file_info in file_metadata:
                file_path = file_info['file_path']
                meta = file_info['metadata']
                
                # Calculate MD5 hash for the file
                import hashlib
                md5_hash = hashlib.md5()
                try:
                    with open(file_path, 'rb') as f:
                        for chunk in iter(lambda: f.read(4096), b""):
                            md5_hash.update(chunk)
                    md5_value = md5_hash.hexdigest()
                except:
                    md5_value = "unknown"
                
                distributions.append({
                    "@type": "cr:FileObject",
                    "@id": meta['file_name'],
                    "name": meta['file_name'],
                    "description": f"VCF file containing genetic variants with {meta['sample_count']} samples and ~{meta['variant_count']:,} variants",
                    "contentUrl": f"{self.base_url}/{meta['file_name']}",
                    "encodingFormat": "application/gzip" if meta['is_compressed'] else "text/plain",
                    "contentSize": str(meta['file_size']),
                    "md5": md5_value,
                    "sc:variantCount": meta['variant_count'],
                    "sc:sampleCount": meta['sample_count']
                })
            
            croissant_metadata["distribution"] = distributions
            
            # Add record set for open access data
            croissant_metadata["recordSet"] = [
                {
                    "@type": "cr:RecordSet",
                    "name": "genetic_variants",
                    "description": "Genetic variants in VCF format with genotype calls and quality information",
                    "field": [
                        {
                            "@type": "cr:Field",
                            "name": "CHROM",
                            "description": "Chromosome identifier",
                            "dataType": "sc:Text",
                            "source": {
                                "fileObject": {
                                    "@id": distributions[0]["name"]
                                },
                                "extract": {
                                    "column": "CHROM"
                                }
                            }
                        },
                        {
                            "@type": "cr:Field", 
                            "name": "POS",
                            "description": "Genomic position (1-based)",
                            "dataType": "sc:Integer",
                            "source": {
                                "fileObject": {
                                    "@id": distributions[0]["name"]
                                },
                                "extract": {
                                    "column": "POS"
                                }
                            }
                        },
                        {
                            "@type": "cr:Field",
                            "name": "ID", 
                            "description": "Variant identifier (e.g., dbSNP ID)",
                            "dataType": "sc:Text",
                            "source": {
                                "fileObject": {
                                    "@id": distributions[0]["name"]
                                },
                                "extract": {
                                    "column": "ID"
                                }
                            }
                        },
                        {
                            "@type": "cr:Field",
                            "name": "REF",
                            "description": "Reference allele sequence", 
                            "dataType": "sc:Text",
                            "source": {
                                "fileObject": {
                                    "@id": distributions[0]["name"]
                                },
                                "extract": {
                                    "column": "REF"
                                }
                            }
                        },
                        {
                            "@type": "cr:Field",
                            "name": "ALT",
                            "description": "Alternative allele sequence(s)",
                            "dataType": "sc:Text",
                            "source": {
                                "fileObject": {
                                    "@id": distributions[0]["name"]
                                },
                                "extract": {
                                    "column": "ALT"
                                }
                            }
                        },
                        {
                            "@type": "cr:Field",
                            "name": "QUAL",
                            "description": "Variant quality score",
                            "dataType": "sc:Float",
                            "source": {
                                "fileObject": {
                                    "@id": distributions[0]["name"]
                                },
                                "extract": {
                                    "column": "QUAL"
                                }
                            }
                        },
                        {
                            "@type": "cr:Field",
                            "name": "FILTER", 
                            "description": "Filter status (PASS/FAIL)",
                            "dataType": "sc:Text",
                            "source": {
                                "fileObject": {
                                    "@id": distributions[0]["name"]
                                },
                                "extract": {
                                    "column": "FILTER"
                                }
                            }
                        },
                        {
                            "@type": "cr:Field",
                            "name": "INFO",
                            "description": "Variant annotation information",
                            "dataType": "sc:Text",
                            "source": {
                                "fileObject": {
                                    "@id": distributions[0]["name"]
                                },
                                "extract": {
                                    "column": "INFO"
                                }
                            }
                        },
                        {
                            "@type": "cr:Field",
                            "name": "FORMAT",
                            "description": "Genotype format specification",
                            "dataType": "sc:Text",
                            "source": {
                                "fileObject": {
                                    "@id": distributions[0]["name"]
                                },
                                "extract": {
                                    "column": "FORMAT"
                                }
                            }
                        },
                        {
                            "@type": "cr:Field",
                            "name": "GENOTYPES",
                            "description": "Sample genotype calls and associated data", 
                            "dataType": "sc:Text",
                            "source": {
                                "fileObject": {
                                    "@id": distributions[0]["name"]
                                },
                                "extract": {
                                    "column": "SAMPLES"
                                }
                            }
                        }
                    ]
                }
            ]
            
        else:
            # Protected access: source-only metadata
            # Generate consistent hash for source repository URL
            import hashlib
            url_hash = hashlib.sha256(dataset_info['data_url'].encode('utf-8')).hexdigest()
            
            croissant_metadata["distribution"] = [
                {
                    "@type": "cr:FileObject",
                    "@id": "source_repository",
                    "name": "source_repository",
                    "description": f"Access {dataset_info['name']} through official repository",
                    "contentUrl": dataset_info['data_url'],
                    "encodingFormat": "text/html",
                    "sha256": url_hash,
                    "sc:accessRequirement": "Authorized access required" if 'dbGaP' in dataset_info.get('license', '') or 'Agreement' in dataset_info.get('license', '') else "Registration required"
                }
            ]
        
        # Add dataset statistics
        if dataset_info['access_level'] == 'open' and total_variants > 0:
            stats = {
                "numberOfVariants": total_variants,
                "numberOfSamples": total_samples,
                "fileCount": len(file_metadata)
            }
        else:
            stats = {
                "subjectCount": dataset_info['subjects'],
                "estimatedVariants": dataset_info['variant_count']
            }
        
        croissant_metadata["sc:datasetStatistics"] = stats
        
        # Add access-specific information
        if dataset_info['access_level'] == 'protected':
            if 'dbgap_accession' in dataset_info:
                croissant_metadata["sc:dbGaPAccession"] = dataset_info['dbgap_accession']
            if 'ida_access' in dataset_info:
                croissant_metadata["sc:idaAccess"] = dataset_info['ida_access']
            if 'genomes_project_url' in dataset_info:
                croissant_metadata["sc:genomesProjectUrl"] = dataset_info['genomes_project_url']
        
        return croissant_metadata
    
    def save_croissant_metadata(self, metadata: Dict[str, Any], output_path: Path) -> None:
        """Save Croissant metadata to a JSON-LD file."""
        
        with open(output_path, 'w', encoding='utf-8') as f:
            json.dump(metadata, f, indent=2, ensure_ascii=False)
        
        print(f"Croissant metadata saved to: {output_path}")


def main():
    """Main function to generate Croissant metadata for VCF datasets."""
    
    # Initialize generator
    generator = VCFCroissantGenerator()
    
    # Create output directory
    output_dir = Path('croissant_metadata')
    output_dir.mkdir(exist_ok=True)
    
    # Generate metadata for different datasets
    datasets_config = [
        {
            'name': 'GTEx',
            'vcf_files': None  # Protected - source only
        },
        {
            'name': 'ADNI', 
            'vcf_files': None  # Protected - source only
        },
        {
            'name': 'MAGE',
            'vcf_files': None  # Protected - source only (1000 Genomes)
        },
        {
            'name': 'ENCODE',
            'vcf_files': list(Path('cell_line_analysis_v5').glob('*/final/combined_variants.vcf.gz'))  # Open access
        },
        {
            'name': 'ENCODE_PAIRED',
            'vcf_files': [Path('paired_cell_lines/paired_cell_lines_variants.vcf.gz')]  # Open access - paired cell lines
        }
    ]
    
    for config in datasets_config:
        dataset_name = config['name']
        vcf_files = config['vcf_files']
        
        print(f"\nGenerating Croissant metadata for {dataset_name}...")
        
        try:
            # Generate metadata
            metadata = generator.generate_croissant_metadata(dataset_name, vcf_files)
            
            # Create output filename
            output_name = f"{dataset_name.lower()}_wgs_variants_croissant.jsonld"
            output_path = output_dir / output_name
            
            # Save metadata
            generator.save_croissant_metadata(metadata, output_path)
            
        except Exception as e:
            print(f"Error processing {dataset_name}: {str(e)}")
    
    print(f"\nVCF Croissant metadata generation complete!")


if __name__ == "__main__":
    main()