#!/usr/bin/env python3
"""
vcf_to_genepop.py

Convert one (or two, via --simulation) VCFs + popmap into a single
NewHybrids‐formatted file, with user-specified reference labels, and
produce an index-to-sample mapping file.

Usage:
    vcf_to_genepop.py \
      --vcf MAIN_VCF \
      --popmap POPMAP.tsv \
      --prefix PREFIX \
      --p0-label POP0_NAME \
      --p1-label POP1_NAME \
      [--simulation SIM_VCF]

Produces:
    <PREFIX>_NewHybrids.txt
    <PREFIX>_NH_index_map.tsv
"""
import sys
import re
import argparse
import pandas as pd

def parse_vcf(vcf_path):
    samples = []
    loci = []
    gen6 = {}
    gt_split = re.compile(r"[\/|]")

    with open(vcf_path) as fh:
        for line in fh:
            if line.startswith('#CHROM'):
                cols = line.strip().split('\t')
                samples = cols[9:]
                gen6 = {s: [] for s in samples}
                break

        for line in fh:
            if line.startswith('#'): continue
            cols = line.strip().split('\t')
            chrom, pos, vid = cols[0], cols[1], cols[2]
            loci.append(vid if vid != '.' else f"{chrom}_{pos}")

            fmt = cols[8].split(':')
            gt_idx = fmt.index('GT') if 'GT' in fmt else 0

            for s, sf in zip(samples, cols[9:]):
                gt = sf.split(':')[gt_idx]
                if gt in ('./.', '.|.'):
                    code = '000000'
                else:
                    parts = []
                    for a in re.split(gt_split, gt):
                        if a == '0': parts.append('001')
                        elif a == '1': parts.append('002')
                        else: parts.append('000')
                    if len(parts) == 1:
                        parts *= 2
                    code = ''.join(parts)
                gen6[s].append(code)

    return samples, loci, gen6

def generate_lumped(gen6):
    lumped = {}
    for s, arr in gen6.items():
        two = []
        for c in arr:
            a1, a2 = c[:3], c[3:]
            def d(x):
                if x == '001': return '1'
                if x == '002': return '2'
                return '0'
            two.append(d(a1) + d(a2))
        lumped[s] = two
    return lumped

def write_newhybrids(path, loci, pure0, pure1, unknown, sim, lumped):
    with open(path, 'w') as f:
        inds = []
        ztag = []
        for s in pure0:
            inds.append(s); ztag.append('z0')
        for s in pure1:
            inds.append(s); ztag.append('z1')
        for s in unknown:
            inds.append(s); ztag.append('')
        for s in sim:
            inds.append(s); ztag.append('')

        f.write(f"NumIndivs {len(inds)}\n")
        f.write(f"NumLoci {len(loci)}\n")
        f.write("Digits 1\nFormat Lumped\n\n")
        f.write("LocusNames " + " ".join(loci) + "\n")
        for i, (s, z) in enumerate(zip(inds, ztag), start=1):
            prefix = f"{i} {z} " if z else f"{i} "
            f.write(prefix + " ".join(lumped[s]) + "\n")
    return inds

def main(vcf_path, popmap_path, prefix, p0_label, p1_label, sim_vcf):
    # Load popmap
    popdf = pd.read_csv(
        popmap_path, sep='\t', header=None,
        names=['sample','pop'], dtype=str
    ).set_index('sample')

    # Verify user-specified labels
    pops = popdf['pop'].unique().tolist()
    if p0_label not in pops or p1_label not in pops:
        sys.exit(f"Error: --p0-label '{p0_label}' and --p1-label '{p1_label}' must both be present in the popmap populations {pops}")
    if p0_label == p1_label:
        sys.exit("Error: --p0-label and --p1-label must be different")

    # Parse main VCF
    main_samples, loci, gen6 = parse_vcf(vcf_path)

    # Parse simulation VCF if provided
    if sim_vcf:
        sim_samples, loci2, sim6 = parse_vcf(sim_vcf)
        if loci2 != loci:
            sys.exit("Error: loci in simulation VCF do not match main VCF")
        gen6.update(sim6)
    else:
        sim_samples = []

    # Assign samples based on labels
    pure0   = [s for s in main_samples if popdf.at[s,'pop'] == p0_label]
    pure1   = [s for s in main_samples if popdf.at[s,'pop'] == p1_label]
    unknown = [s for s in main_samples if s not in pure0 + pure1]

    # Generate lumped genotypes
    lumped = generate_lumped(gen6)

    # Write NewHybrids file and capture ordered inds
    out_path = f"{prefix}_NewHybrids.txt"
    inds = write_newhybrids(out_path, loci, pure0, pure1, unknown, sim_samples, lumped)
    print(f"Wrote NewHybrids file: {out_path}")

    # Write index-to-sample mapping
    map_path = f"{prefix}_NH_index_map.tsv"
    with open(map_path, 'w') as mf:
        mf.write("Index\tSample\n")
        for idx, s in enumerate(inds, start=1):
            mf.write(f"{idx}\t{s}\n")
    print(f"Wrote NewHybrids index map: {map_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="VCF → single NewHybrids file with custom reference labels"
    )
    parser.add_argument(
        '--vcf', required=True, help="Path to the main VCF file"
    )
    parser.add_argument(
        '--popmap', required=True, help="Tab-delimited popmap: SAMPLE<TAB>POP"
    )
    parser.add_argument(
        '--prefix', required=True, help="Output prefix"
    )
    parser.add_argument(
        '--p0-label', required=True, help="Population name to use as P0 in the popmap"
    )
    parser.add_argument(
        '--p1-label', required=True, help="Population name to use as P1 in the popmap"
    )
    parser.add_argument(
        '--simulation', help="Optional VCF of simulated admixed samples"
    )
    args = parser.parse_args()

    main(
        vcf_path    = args.vcf,
        popmap_path = args.popmap,
        prefix      = args.prefix,
        p0_label    = args.p0_label,
        p1_label    = args.p1_label,
        sim_vcf     = args.simulation
    )
