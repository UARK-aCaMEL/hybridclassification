#!/usr/bin/env python3
import sys
import re
import pandas as pd


def parse_vcf(vcf_path):
    """
    Parse VCF, return sample order, locus names, and genotype codes per sample.
    Genepop codes: 0/0->0101, 0/1 or 1/0->0102, 1/1->0202, missing->0000
    """
    samples = []
    loci = []
    genotypes = {}
    gt_pattern = re.compile(r"[/|]")

    with open(vcf_path) as fh:
        for line in fh:
            if line.startswith('#CHROM'):
                cols = line.strip().split('\t')
                samples = cols[9:]
                genotypes = {s: [] for s in samples}
                break
        # read variant lines
        for line in fh:
            if line.startswith('#'):
                continue
            cols = line.strip().split('\t')
            chrom, pos, vid, ref, alt = cols[0], cols[1], cols[2], cols[3], cols[4]
            loci.append(vid if vid != '.' else f"{chrom}_{pos}")
            fmt = cols[8].split(':')
            try:
                gt_idx = fmt.index('GT')
            except ValueError:
                gt_idx = 0
            for s, sf in zip(samples, cols[9:]):
                fields = sf.split(':')
                gt = fields[gt_idx]
                if gt in ('./.', '.|.'):
                    code = '0000'
                else:
                    alleles = re.split(gt_pattern, gt)
                    parts = []
                    for a in alleles:
                        if a == '0': parts.append('01')
                        elif a == '1': parts.append('02')
                        else: parts.append('00')
                    if len(parts) == 1:
                        parts = [parts[0], parts[0]]
                    code = ''.join(parts)
                genotypes[s].append(code)
    return samples, loci, genotypes


def write_genepop(path, title, loci, groups, genotypes):
    """
    Write a GENEPOP file:
    - title: first line
    - loci names, comma-separated
    - for each population: 'Pop' line, then 'sample, code code ...'
    groups: ordered list of (pop_label, list_of_samples)
    """
    with open(path, 'w') as out:
        out.write(title + '\n')
        out.write(','.join(loci) + '\n')
        for pop_label, samples in groups:
            if not samples:
                continue
            out.write('Pop\n')
            for s in samples:
                geno_str = ' '.join(genotypes[s])
                out.write(f"{s}, {geno_str}\n")


def main(vcf_path, popmap_path, prefix):
    # load popmap: two-column, tab-delim, no header: SAMPLE \t POP
    popdf = pd.read_csv(popmap_path, sep='\t', header=None, names=['sample','pop'], dtype=str)
    popdf = popdf.set_index('sample')
    # determine two base populations
    unique_pops = popdf['pop'].unique().tolist()
    if len(unique_pops) < 2:
        sys.exit('Error: need at least two populations in popmap')
    p0, p1 = unique_pops[:2]

    samples, loci, genotypes = parse_vcf(vcf_path)

    # assign samples
    pure0 = [s for s in samples if s in popdf.index and popdf.at[s,'pop'] == p0]
    pure1 = [s for s in samples if s in popdf.index and popdf.at[s,'pop'] == p1]
    admix = [s for s in samples if s not in pure0 + pure1]

    # write pure two-population file
    pure_path = f"{prefix}_pure.txt"
    write_genepop(
        pure_path,
        title=f"{prefix} two-population ({p0} vs {p1})",
        loci=loci,
        groups=[('P0', pure0), ('P1', pure1)],
        genotypes=genotypes
    )
    print(f"Wrote pure file: {pure_path}")

    # write admixed file (single population 'ADMIX')
    admix_path = f"{prefix}_admix.txt"
    if admix:
        write_genepop(
            admix_path,
            title=f"{prefix} admixed individuals",
            loci=loci,
            groups=[('ADMIX', admix)],
            genotypes=genotypes
        )
        print(f"Wrote admixed file: {admix_path}")
    else:
        # create empty file
        open(admix_path, 'w').close()
        print(f"No admixed individuals; created empty file {admix_path}")


if __name__ == '__main__':
    if len(sys.argv) != 4:
        sys.exit('Usage: vcf_to_genepop.py <input.vcf> <popmap.txt> <prefix>')
    main(sys.argv[1], sys.argv[2], sys.argv[3])
