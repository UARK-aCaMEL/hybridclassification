#!/usr/bin/env python3
"""
Simulate hybrids (Pure, F1, F2, BC1, BC2) via two strategies:
  - freq: sample genotypes by allele frequencies (default)
  - sampled: draw two parents and combine alleles at random

Usage:
    python hybrid_simulator.py \
      --vcf input.vcf.gz \
      --popmap ref_popmap.tsv \
      --p0 POP0_ID \
      --p1 POP1_ID \
      --num_reps 3 \
      --size_pure 200 \
      --size_f1 200 \
      --size_f2 200 \
      --size_bc 200 \
      --strategy freq|sampled \
      --out_prefix sim
"""
import argparse
import numpy as np
import pandas as pd
import pysam
import sys
import random

def parse_args():
    p = argparse.ArgumentParser(description="Simulate hybrids → VCF + popmap")
    p.add_argument("--vcf",       required=True, help="bgzipped/tabix'd VCF")
    p.add_argument("--popmap",    required=True, help="TSV (no header): SAMPLE<TAB>POP")
    p.add_argument("--p0",        required=True, help="ref-pop0 ID")
    p.add_argument("--p1",        required=True, help="ref-pop1 ID")
    p.add_argument("--num_reps",  type=int, default=1,   help="number of replicates")
    p.add_argument("--size_pure", type=int, default=200, help="pure-pop sample size")
    p.add_argument("--size_f1",   type=int, default=200, help="F1 sample size")
    p.add_argument("--size_f2",   type=int, default=200, help="F2 sample size")
    p.add_argument("--size_bc",   type=int, default=200, help="backcross sample size")
    p.add_argument("--strategy",  choices=["freq","sampled"], default="freq",
                   help="Simulation strategy: 'freq' (allele-frequency) or 'sampled' (resample parents)")
    p.add_argument("--out_prefix",default="sim",          help="output prefix")
    return p.parse_args()

def compute_ref_af(vcf_path, samples, popmap, pop_id):
    vcf = pysam.VariantFile(vcf_path)
    idxs = [i for i,s in enumerate(samples) if popmap.get(s)==pop_id]
    freqs = []
    for rec in vcf.fetch():
        ac = 0
        called = 0
        for i in idxs:
            gt = rec.samples[samples[i]]["GT"]
            if gt is None or None in gt or gt[0]<0 or gt[1]<0:
                continue
            ac     += gt[0] + gt[1]
            called += 2
        freqs.append(ac/called if called>0 else 0.0)
    return np.array(freqs, dtype=float)

def sample_genotypes_from_af(p, n):
    L = p.size
    draws = (np.random.rand(n, L, 2) < p[np.newaxis,:,None])
    return draws.sum(axis=2).astype(np.int8)

def write_vcf_manually(template_vcf, sample_ids, geno_mat, out_vcf_path):
    raw_hdr = str(template_vcf.header).splitlines()
    meta_lines = [l for l in raw_hdr if l.startswith("##")]
    old_cols = [l for l in raw_hdr if l.startswith("#CHROM")][0].split("\t")[:9]
    new_hdr_line = "\t".join(old_cols + sample_ids)
    with open(out_vcf_path, "w") as fout:
        for l in meta_lines:
            fout.write(l + "\n")
        fout.write(new_hdr_line + "\n")
        L = geno_mat.shape[1]
        for li, rec in enumerate(template_vcf.fetch()):
            if li >= L:
                break
            chrom = rec.contig
            pos   = str(rec.pos)
            rid   = rec.id or "."
            ref   = rec.ref
            alt   = ",".join(rec.alts)
            qual  = str(rec.qual) if rec.qual is not None else "."
            filt  = ";".join(rec.filter.keys()) if rec.filter.keys() else "PASS"
            info_items = rec.info.items()
            info = "." if not info_items else ";".join(f"{k}={v}" for k,v in info_items)
            fmt  = "GT"
            fields = [chrom, pos, rid, ref, alt, qual, filt, info, fmt]
            for g in geno_mat[:, li]:
                if   g == 0: gt_str = "0/0"
                elif g == 1: gt_str = "0/1"
                else:         gt_str = "1/1"
                fields.append(gt_str)
            fout.write("\t".join(fields) + "\n")

def main():
    args = parse_args()
    pm_df = pd.read_csv(args.popmap, sep="\t", header=None, names=["sample","pop"] )
    popmap = dict(zip(pm_df["sample"], pm_df["pop"]))
    vcf_in = pysam.VariantFile(args.vcf)
    samples = list(vcf_in.header.samples)

    if args.strategy == "freq":
        af0 = compute_ref_af(args.vcf, samples, popmap, args.p0)
        af1 = compute_ref_af(args.vcf, samples, popmap, args.p1)
        L   = af0.size
        all_genos   = []
        popmap_rows = []
        sample_ids  = []
        for rep in range(1, args.num_reps+1):
            pure0 = sample_genotypes_from_af(af0, args.size_pure)
            pure1 = sample_genotypes_from_af(af1, args.size_pure)
            pF0 = pure0.sum(axis=0)/(2*args.size_pure)
            pF1 = pure1.sum(axis=0)/(2*args.size_pure)
            f1 = (np.random.rand(args.size_f1, L)<af0).astype(int) + \
                 (np.random.rand(args.size_f1, L)<af1).astype(int)
            pF1_hat = f1.sum(axis=0)/(2*args.size_f1)
            f2 = sample_genotypes_from_af(pF1_hat, args.size_f2)
            bc0 = (np.random.rand(args.size_bc, L)<af0).astype(int) + \
                  (np.random.rand(args.size_bc, L)<pF1_hat).astype(int)
            bc1 = (np.random.rand(args.size_bc, L)<af1).astype(int) + \
                  (np.random.rand(args.size_bc, L)<pF1_hat).astype(int)
            cats = [
                (pure0, f"{rep}_Pure_{args.p0}", f"Pure_{args.p0}"),
                (pure1, f"{rep}_Pure_{args.p1}", f"Pure_{args.p1}"),
                (f1,    f"{rep}_F1",              "F1"),
                (f2,    f"{rep}_F2",              "F2"),
                (bc0,   f"{rep}_BC_{args.p0}",     f"BC_{args.p0}"),
                (bc1,   f"{rep}_BC_{args.p1}",     f"BC_{args.p1}"),
            ]
            for mat, sid_pref, pop_lbl in cats:
                for i in range(mat.shape[0]):
                    sid = f"{sid_pref}_{i+1}"
                    sample_ids.append(sid)
                    popmap_rows.append({"ID":sid, "POP": pop_lbl})
                    all_genos.append(mat[i])
        geno_array = np.vstack(all_genos)

    else:  # sampled strategy
        records = [r for r in vcf_in.fetch()]
        L = len(records)
        # map each locus to {sample: (allele1,allele2)}
        geno_haps = [ {s: records[i].samples[s]["GT"] for s in samples} for i in range(L) ]
        p0_samples = [s for s in samples if popmap.get(s)==args.p0]
        p1_samples = [s for s in samples if popmap.get(s)==args.p1]
        all_haps = []
        popmap_rows = []
        sample_ids = []
        for rep in range(1, args.num_reps+1):
            # Pure p0
            for i in range(args.size_pure):
                hap = np.zeros((L,2), dtype=np.int8)
                parent = random.choice(p0_samples)
                for li in range(L):
                    gt = geno_haps[li].get(parent)
                    a0, a1 = (gt if gt and None not in gt else (0,0))
                    hap[li,0], hap[li,1] = a0, a1
                all_haps.append(hap)
                sid = f"{rep}_Pure_{args.p0}_{i+1}"
                sample_ids.append(sid)
                popmap_rows.append({"ID":sid, "POP": f"Pure_{args.p0}"})
            # Pure p1
            for i in range(args.size_pure):
                hap = np.zeros((L,2), dtype=np.int8)
                parent = random.choice(p1_samples)
                for li in range(L):
                    gt = geno_haps[li].get(parent)
                    a0, a1 = (gt if gt and None not in gt else (0,0))
                    hap[li,0], hap[li,1] = a0, a1
                all_haps.append(hap)
                sid = f"{rep}_Pure_{args.p1}_{i+1}"
                sample_ids.append(sid)
                popmap_rows.append({"ID":sid, "POP": f"Pure_{args.p1}"})
            # F1: cross one random P0 and one random P1
            f1_haps = []
            for i in range(args.size_f1):
                hap = np.zeros((L,2), dtype=np.int8)
                p0 = random.choice(p0_samples)
                p1 = random.choice(p1_samples)
                for li in range(L):
                    gt0 = geno_haps[li].get(p0)
                    gt1 = geno_haps[li].get(p1)
                    allele0 = gt0[random.choice([0,1])] if gt0 and None not in gt0 else 0
                    allele1 = gt1[random.choice([0,1])] if gt1 and None not in gt1 else 0
                    hap[li,0], hap[li,1] = allele0, allele1
                f1_haps.append(hap)
                all_haps.append(hap)
                sid = f"{rep}_F1_{i+1}"
                sample_ids.append(sid)
                popmap_rows.append({"ID":sid, "POP": "F1"})
            # F2: cross two random F1
            for i in range(args.size_f2):
                hap = np.zeros((L,2), dtype=np.int8)
                pa, pb = random.choice(f1_haps), random.choice(f1_haps)
                for li in range(L):
                    allele0 = pa[li, random.choice([0,1])]
                    allele1 = pb[li, random.choice([0,1])]
                    hap[li,0], hap[li,1] = allele0, allele1
                all_haps.append(hap)
                sid = f"{rep}_F2_{i+1}"
                sample_ids.append(sid)
                popmap_rows.append({"ID":sid, "POP": "F2"})
            # Backcross to P0 (BC_{p0})
            for i in range(args.size_bc):
                hap = np.zeros((L,2), dtype=np.int8)
                f1p = random.choice(f1_haps)
                for li in range(L):
                    a1 = f1p[li, random.choice([0,1])]
                    gt0 = geno_haps[li].get(random.choice(p0_samples))
                    a0 = gt0[random.choice([0,1])] if gt0 and None not in gt0 else 0
                    hap[li,0], hap[li,1] = a0, a1
                all_haps.append(hap)
                sid = f"{rep}_BC_{args.p0}_{i+1}"
                sample_ids.append(sid)
                popmap_rows.append({"ID":sid, "POP": f"BC_{args.p0}"})
            # Backcross to P1 (BC_{p1})
            for i in range(args.size_bc):
                hap = np.zeros((L,2), dtype=np.int8)
                f1p = random.choice(f1_haps)
                for li in range(L):
                    a1 = f1p[li, random.choice([0,1])]
                    gt1 = geno_haps[li].get(random.choice(p1_samples))
                    a0 = gt1[random.choice([0,1])] if gt1 and None not in gt1 else 0
                    hap[li,0], hap[li,1] = a0, a1
                all_haps.append(hap)
                sid = f"{rep}_BC_{args.p1}_{i+1}"
                sample_ids.append(sid)
                popmap_rows.append({"ID":sid, "POP": f"BC_{args.p1}"})
        # convert to genotype counts
        geno_array = np.array([h.sum(axis=1) for h in all_haps], dtype=np.int8)

    # write popmap TSV
    df_popmap = pd.DataFrame(popmap_rows)
    popmap_file = f"{args.out_prefix}_simulation.tsv"
    df_popmap.to_csv(popmap_file, sep="\t", index=False)
    print(f"Wrote population map: {popmap_file}")

    # write VCF manually
    write_vcf_manually(vcf_in, sample_ids, geno_array, f"{args.out_prefix}_simulation.vcf")

if __name__ == "__main__":
    main()
