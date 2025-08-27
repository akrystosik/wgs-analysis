# DeepVariant Setup for ENCODE WGS Analysis

## Testing DeepVariant Access
1. Created test job to verify container access
   - Image: google/deepvariant:1.6.0
   - Location: jobs/deepvariant-test.yaml
   - Result: Successfully pulled and ran despite policy warning

## Production Setup
1. Job Configuration
   - Location: jobs/deepvariant-wgs.yaml
   - Resources: 16 CPU (32 max), 32GB RAM (64GB max)
   - Input: GM23248_WGS.markdup.bam
   - Output: data/variants/GM23248.vcf.gz

## Workflow
1. Wait for BAM processing to complete:
   - Pair1 conversion
   - Pair2 processing
   - Mark duplicates
2. Run DeepVariant job from local Mac
3. Monitor variant calling progress

## Commands
```yaml
# Key DeepVariant parameters
--model_type=WGS
--num_shards=16
--ref=/path/to/reference
--reads=/path/to/bam
