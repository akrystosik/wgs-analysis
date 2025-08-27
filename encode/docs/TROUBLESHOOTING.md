# Troubleshooting Guide

## Common Issues and Solutions

### Installation Issues

#### DeepVariant/DeepSomatic Not Found
**Error:** `Command not found: deepvariant`
**Solution:**
```bash
# Ensure DeepVariant is in PATH
export PATH=$PATH:/path/to/deepvariant/bin
# Or use full path in config
```

#### Missing BAM Index
**Error:** `[E::hts_idx_load3] Could not load index`
**Solution:**
```bash
samtools index your_file.bam
```

#### Picard Memory Errors
**Error:** `OutOfMemoryError: Java heap space`
**Solution:**
```bash
# Increase Java heap size
export JAVA_OPTS="-Xmx64g"
# Or modify pipeline config memory settings
```

### Pipeline Issues

#### Variant Calling Fails
**Error:** `make_examples failed with exit code 1`
**Solution:**
1. Check BAM file integrity: `samtools quickcheck file.bam`
2. Verify reference genome index exists
3. Ensure sufficient disk space
4. Check log files in `logs/` directory

#### Permission Denied
**Error:** `Permission denied: /tmp/deepvariant_tmp`
**Solution:**
```bash
# Create writable temp directory
mkdir -p /tmp/deepvariant_tmp
chmod 755 /tmp/deepvariant_tmp
```

#### Resource Exhaustion
**Error:** `Cannot allocate memory` or process killed
**Solution:**
1. Reduce parallelization: edit `num_shards` in config
2. Increase system memory allocation
3. Use disk-based temporary files

### Configuration Issues

#### Invalid YAML Config
**Error:** `yaml.scanner.ScannerError`
**Solution:**
```bash
# Validate YAML syntax
python -c "import yaml; yaml.safe_load(open('config/pipeline_config.yaml'))"
```

#### Path Not Found
**Error:** `FileNotFoundError: [Errno 2] No such file or directory`
**Solution:**
1. Use absolute paths in configuration
2. Verify all input files exist
3. Check file permissions

### Performance Issues

#### Slow Processing
**Symptoms:** Pipeline takes >24 hours per sample
**Solutions:**
1. Increase CPU allocation (`num_threads`)
2. Optimize shard count (`num_shards`)
3. Use SSD storage for temporary files
4. Monitor system resources: `htop`, `iostat`

#### High Memory Usage
**Symptoms:** System becomes unresponsive
**Solutions:**
1. Reduce memory per process
2. Process samples sequentially instead of parallel
3. Implement disk-based sorting

### Quality Control Issues

#### Low Coverage Warnings
**Warning:** `Mean coverage < 20x`
**Solutions:**
1. Check FASTQ quality with FastQC
2. Verify alignment parameters
3. Examine duplicate rates

#### High Duplicate Rates
**Warning:** `Duplicate rate > 30%`
**Investigation:**
1. Check library preparation protocol
2. Examine PCR amplification cycles
3. Consider using different duplicate marking strategy

### Data Issues

#### VCF Format Errors
**Error:** `Invalid VCF header` or `Malformed VCF`
**Solution:**
```bash
# Validate VCF format
bcftools view -h file.vcf.gz | head -20
# Fix with bcftools
bcftools norm -f reference.fa input.vcf > fixed.vcf
```

#### Missing Chromosomes
**Error:** `Chromosome 'chr1' not found in reference`
**Solution:**
1. Check chromosome naming convention (chr1 vs 1)
2. Ensure BAM and reference match
3. Use consistent reference genome version

## Getting Help

### Log File Analysis
Always check these log files first:
- `logs/pipeline_[timestamp]/pipeline_execution.log`
- `logs/pipeline_[timestamp]/pipeline_errors.log`
- Individual tool logs in respective directories

### Debug Mode
Enable verbose logging:
```bash
python -m pipeline.pipeline --debug --verbose
```

### System Resource Monitoring
```bash
# Monitor during pipeline execution
htop
iostat -x 1
df -h
```

### Contact and Support
- GitHub Issues: https://github.com/akrystosik/encode_analysis/issues
- Include log files and system information in issue reports