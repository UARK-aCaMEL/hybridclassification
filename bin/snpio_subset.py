#!/usr/bin/env python3
import argparse
import os
from snpio import VCFReader


def get_prefix_from_vcf_path(vcf_path):
    basename = os.path.basename(vcf_path)
    for ext in [".vcf.gz", ".vcf"]:
        if basename.endswith(ext):
            return basename[: -len(ext)]
    return basename  # fallback if no known extension


def main():
    parser = argparse.ArgumentParser(
        description="Run SNPio to filter individuals and generate missingness reports"
    )
    parser.add_argument("--vcf", required=True, help="Path to the VCF file.")
    parser.add_argument("--popmap", required=True, help="Path to the population map file.")
    parser.add_argument("--include", action="append", required=True,
                        help="Population to include (use multiple times for multiple pops)")
    parser.add_argument("--prefix", type=str, default=None,
                        help="Prefix for output files (default: derived from VCF name)")

    args = parser.parse_args()

    # Use user-specified prefix or fallback to VCF-derived
    prefix = args.prefix or get_prefix_from_vcf_path(args.vcf)

    # read data
    gd = VCFReader(
        filename=args.vcf,
        popmapfile=args.popmap,
        force_popmap=True,
        verbose=True,
        plot_format="png",
        plot_fontsize=8,
        plot_dpi=300,
        prefix=prefix,
        include_pops=args.include
    )

    output_vcf = f"{prefix}.subset.vcf"
    gd.write_vcf(output_vcf)


if __name__ == "__main__":
    main()
