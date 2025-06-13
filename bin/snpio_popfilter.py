#!/usr/bin/env python3
import argparse
import os
from snpio import NRemover2, VCFReader


def get_prefix_from_vcf_path(vcf_path):
    basename = os.path.basename(vcf_path)
    for ext in [".vcf.gz", ".vcf"]:
        if basename.endswith(ext):
            return basename[: -len(ext)]
    return basename


def main():
    parser = argparse.ArgumentParser(
        description="Filter a VCF by population-level missingness using SNPio"
    )
    parser.add_argument(
        "--vcf", required=True,
        help="Path to the input VCF file (vcf or vcf.gz)"
    )
    parser.add_argument(
        "--popmap", required=True,
        help="Path to the population map file (tab-delimited: sample\tpop)"
    )
    parser.add_argument(
        "--pop_cov", type=float, default=0.5,
        help="Maximum allowed proportion of missing data per population (default: 0.5)"
    )
    parser.add_argument(
        "--prefix", type=str, default=None,
        help="Prefix for output files (default: derived from VCF name)"
    )

    args = parser.parse_args()

    prefix = args.prefix or get_prefix_from_vcf_path(args.vcf)

    # Load VCF and popmap
    gd = VCFReader(
        filename=args.vcf,
        popmapfile=args.popmap,
        force_popmap=True,
        verbose=True,
        plot_format="png",
        plot_fontsize=8,
        plot_dpi=300,
        prefix=prefix,
    )

    # Filter by population missingness
    nrm = NRemover2(gd)
    gd_filt = nrm.filter_missing_pop(args.pop_cov).resolve()

    # Write filtered VCF
    output_vcf = f"{prefix}.filtered.vcf"
    gd_filt.write_vcf(output_vcf)

    print(f"Filtered VCF written to: {output_vcf}")


if __name__ == "__main__":
    main()
