#!/usr/bin/env python3
import argparse
import os
import sys
from collections import Counter

from snpio import NRemover2, VCFReader


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
    parser.add_argument(
        "--popmap", required=True, help="Path to the population map file."
    )
    parser.add_argument(
        "--ind_cov",
        type=float,
        default=0.75,
        help="Maximum allowed missingness per individual (default: 0.75)",
    )
    parser.add_argument(
        "--flank_dist",
        type=int,
        default=75,
        help="Maximum allowed distance between SNPs (default: 75)",
    )
    parser.add_argument(
        "--min_maf",
        type=float,
        default=0.05,
        help="Maximum MAF to retain a SNP (default: 0.05)",
    )
    parser.add_argument(
        "--snp_cov",
        type=float,
        default=0.9,
        help="Maximum missing data to retain a SNP (default: 0.9)",
    )
    parser.add_argument(
        "--pop_cov",
        type=float,
        default=0.5,
        help="Maximum missing data to retain a SNP, within each species in --popmap (default: 0.5)",
    )
    parser.add_argument(
        "--min_samples_per_pop",
        type=int,
        default=5,
        help="Stop if fewer samples of any species remain after filtering (default: 5)",
    )
    parser.add_argument("--prefix", type=str, default=None,
        help="Prefix for output files (default: derived from VCF name)")
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for thinning (which SNP is kept in each window)",
    )
    parser.add_argument(
        "--save_plots",
        action="store_true",
        help="Also write SNPio's static plot images (not used by the report)",
    )

    args = parser.parse_args()

    # extract prefix from VCF filename
    prefix = args.prefix or get_prefix_from_vcf_path(args.vcf)

    # read data (static plots are skipped unless --save_plots is given)
    gd = VCFReader(
        filename=args.vcf,
        popmapfile=args.popmap,
        force_popmap=True,
        verbose=True,
        plot_format="png",
        plot_fontsize=8,
        plot_dpi=300,
        prefix=prefix,
        save_plots=args.save_plots,
    )

    # generate missingness reports
    gd.missingness_reports()

    # Filter VCF. Per-species missingness comes first: filtered on the pooled
    # samples alone, the retained loci can be genotyped mostly in the larger
    # species, and the per-sample filter then removes the other species
    nrm = NRemover2(gd)
    gd_filt = (
        nrm.filter_missing_pop(args.pop_cov)
        .filter_monomorphic(exclude_heterozygous=False)
        .filter_biallelic(exclude_heterozygous=False)
        .filter_missing(args.snp_cov)
        .filter_maf(args.min_maf)
        .filter_missing_sample(args.ind_cov)
        .thin_loci(remove_all=False, size=args.flank_dist, seed=args.seed)
        .resolve()
    )

    # Stop if either species of the pair (nearly) disappeared: the analysis
    # would otherwise go ahead without it
    popmap = dict(line.split()[:2] for line in open(args.popmap) if line.strip())
    before = Counter(popmap[s] for s in gd.samples)
    after = Counter(popmap[s] for s in gd_filt.samples)
    short = {p: (after.get(p, 0), n) for p, n in before.items() if after.get(p, 0) < args.min_samples_per_pop}
    if short:
        detail = ", ".join(f"{p} {k} of {n}" for p, (k, n) in sorted(short.items()))
        sys.exit(
            f"ERROR: too few samples left after filtering ({detail}; minimum "
            f"{args.min_samples_per_pop}, set by --min_species_samples). Relax "
            "--ind_cov, --snp_cov or --pop_cov, or remove this pair from --combinations."
        )

    nrm.plot_sankey_filtering_report()

    # Write the filtered VCF using the modified prefix
    output_vcf = f"{prefix}.filter.vcf"
    gd_filt.write_vcf(output_vcf)


if __name__ == "__main__":
    main()
