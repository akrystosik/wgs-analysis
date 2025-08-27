#!/usr/bin/env python3
import argparse
import subprocess
import sys

def main():
    parser = argparse.ArgumentParser(
        description="Wrapper script to run DeepSomatic (an extension of DeepVariant)."
    )
    parser.add_argument("--model_type", default="WGS_TUMOR_ONLY",
                        help="DeepSomatic model type (e.g., WGS, WES, PACBIO, ONT, WGS_TUMOR_ONLY, etc.).")
    parser.add_argument("--ref", required=True,
                        help="Path to reference FASTA file.")
    parser.add_argument("--reads_tumor", required=True,
                        help="Path to tumor BAM file.")
    parser.add_argument("--reads_normal",
                        help="Path to normal BAM file (if tumor-normal). Omit for tumor-only.")
    parser.add_argument("--output_vcf", required=True,
                        help="Path to output VCF file.")
    parser.add_argument("--output_gvcf",
                        help="Path to output GVCF file.")
    parser.add_argument("--sample_name_tumor", default="tumor",
                        help="Sample name for tumor in the VCF.")
    parser.add_argument("--sample_name_normal",
                        help="Sample name for normal in the VCF, if tumor-normal.")
    parser.add_argument("--num_shards", type=int, default=64,
                        help="Number of shards (threads). Default = 64.")
    parser.add_argument("--logging_dir", default="/tmp/deepsomatic/logs",
                        help="Directory for logs.")
    parser.add_argument("--intermediate_results_dir", default="/tmp/deepsomatic/intermediate",
                        help="Directory for intermediate files.")
    parser.add_argument("--regions",
                        help="Genome regions to process (e.g., 'chr20' or a BED file). If omitted, runs on whole genome.")
    parser.add_argument("--use_default_pon_filtering", default="false",
                        help="Set to 'true' for default Panel of Normals filtering in tumor-only. Default = 'false'.")
    parser.add_argument("--dry_run", default="false",
                        help="Set to 'true' to print commands without running them.")
    
    args = parser.parse_args()

    # Construct the command
    cmd = [
        "run_deepsomatic",
        f"--model_type={args.model_type}",
        f"--ref={args.ref}",
        f"--reads_tumor={args.reads_tumor}",
        f"--output_vcf={args.output_vcf}",
        f"--num_shards={args.num_shards}",
        f"--logging_dir={args.logging_dir}",
        f"--intermediate_results_dir={args.intermediate_results_dir}",
        f"--sample_name_tumor={args.sample_name_tumor}",
        f"--use_default_pon_filtering={args.use_default_pon_filtering}",
        f"--dry_run={args.dry_run}",
    ]

    # Add normal BAM if provided (tumor-normal)
    if args.reads_normal:
        cmd.append(f"--reads_normal={args.reads_normal}")
        if args.sample_name_normal:
            cmd.append(f"--sample_name_normal={args.sample_name_normal}")

    # Add output gVCF if provided
    if args.output_gvcf:
        cmd.append(f"--output_gvcf={args.output_gvcf}")

    # Restrict regions, if specified
    if args.regions:
        cmd.append(f"--regions={args.regions}")

    # Show the final command line
    print("DeepSomatic command:\n" + " ".join(cmd), file=sys.stderr)

    # Execute
    subprocess.run(cmd, check=True)

if __name__ == "__main__":
    main()
