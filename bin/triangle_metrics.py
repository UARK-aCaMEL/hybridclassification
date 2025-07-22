#!/usr/bin/env python3
"""
triangle_metrics.py
Compute hybrid index (h) & inter‑class heterozygosity (H).
Produces:
  <prefix>_hindex.tsv             (header #Num_loci=…)
  <prefix>_hindex_fixed.tsv       (header #Num_loci=…)
  <prefix>_classification_popmap.tsv
"""

from __future__ import annotations
import argparse, re, sys
from pathlib import Path
import numpy as np
import pandas as pd
import pysam

# ───────── helpers ─────────
def load_popmap(path: Path) -> dict[str, str]:
    df = pd.read_csv(path, sep="\t", header=None,
                     names=["sample", "pop"], dtype=str)
    return dict(zip(df["sample"], df["pop"]))      # <- FIXED

def recode(rec):
    v = np.full(len(rec.samples), np.nan, dtype=np.float32)
    for i, call in enumerate(rec.samples.values()):
        gt = call["GT"]
        if gt and None not in gt and min(gt) >= 0:
            v[i] = gt[0] + gt[1]
    return v

def vcf_matrix(vcf):
    fh = pysam.VariantFile(vcf)
    samp, loci, g = list(fh.header.samples), [], []
    for r in fh.fetch():
        loci.append(r.id if r.id and r.id!="." else f"{r.contig}_{r.pos}")
        g.append(recode(r))
    if not g:
        sys.exit(f"{vcf} is empty")
    return samp, loci, np.vstack(g).T

def build_N(g, af0, af1):
    N = np.full_like(g, np.nan, dtype=np.float32)
    for i in range(g.shape[1]):
        if np.isnan(af0[i]) or np.isnan(af1[i]): continue
        if af0[i] > af1[i]:
            N[g[:,i]==2,i] = 0; N[g[:,i]==0,i] = 2
        else:
            N[g[:,i]==0,i] = 0; N[g[:,i]==2,i] = 2
        N[g[:,i]==1,i] = 1
    return N

def h_H(N, g):
    S = np.isfinite(N).sum(1).astype(float);  S[S==0]=np.nan
    h = np.nansum(N,1)/(2*S);  H = np.nansum(g==1,1)/S
    pmiss = 1 - (S/g.shape[1])
    return h, H, pmiss

def infer_sim(s):
    for t in s.split("_"):
        if re.fullmatch(r"(F1|F2|BC.*|Pure.*)", t, re.I): return t
    return "Sim"

# ───────── main ─────────
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True, type=Path)
    ap.add_argument("--popmap", required=True, type=Path)
    ap.add_argument("--p0", required=True); ap.add_argument("--p1", required=True)
    ap.add_argument("--out_prefix", required=True)
    ap.add_argument("--simulation", type=Path)
    ap.add_argument("--af_diff_min", type=float, default=0.0)
    args = ap.parse_args()

    # ensure output dir exists
    Path(args.out_prefix).parent.mkdir(parents=True, exist_ok=True)

    popmap = load_popmap(args.popmap)
    if args.p0==args.p1 or args.p0 not in popmap.values() or args.p1 not in popmap.values():
        sys.exit("Invalid parental labels")

    samp, loci, g_main = vcf_matrix(args.vcf)
    if args.simulation:
        sim_samp, loci2, g_sim = vcf_matrix(args.simulation)
        if loci2!=loci: sys.exit("Sim VCF loci/order differ")
    else:
        sim_samp, g_sim = [], np.empty((0,len(loci)), dtype=np.float32)

    mask_p0 = np.array([popmap[s]==args.p0 for s in samp])
    mask_p1 = np.array([popmap[s]==args.p1 for s in samp])

    het=(g_main==1).astype(float); hom2=(g_main==2).astype(float)
    af0=(het[mask_p0].sum(0)+2*hom2[mask_p0].sum(0))/(2*np.isfinite(g_main[mask_p0]).sum(0))
    af1=(het[mask_p1].sum(0)+2*hom2[mask_p1].sum(0))/(2*np.isfinite(g_main[mask_p1]).sum(0))

    diff=np.abs(af0-af1)
    keep_user  = diff>=args.af_diff_min
    keep_fixed = diff>=0.999

    unknown=[s for s in samp if popmap[s] not in (args.p0,args.p1)]
    g_class=np.concatenate((g_main[~(mask_p0|mask_p1)], g_sim), axis=0)
    class_samp=unknown+sim_samp

    # ----- all loci -----
    N=build_N(g_class[:,keep_user], af0[keep_user], af1[keep_user])
    h,H,pm=h_H(N, g_class[:,keep_user])
    df=pd.DataFrame({"Sample":class_samp,"HybridIndex":h,"Heterozygosity":H,"PercMissing":pm})
    with open(f"{args.out_prefix}_hindex.tsv","w") as f:
        f.write(f"#Num_loci={keep_user.sum()}\n"); df.to_csv(f,sep="\t",index=False)

    # ----- fixed -----
    if keep_fixed.sum():
        Nf=build_N(g_class[:,keep_fixed], af0[keep_fixed], af1[keep_fixed])
        hf,Hf,pf=h_H(Nf, g_class[:,keep_fixed])
        dff=pd.DataFrame({"Sample":class_samp,"HybridIndex":hf,"Heterozygosity":Hf,"PercMissing":pf})
    else:
        dff=pd.DataFrame({"Sample":class_samp,"HybridIndex":np.nan,"Heterozygosity":np.nan,"PercMissing":np.nan})
    with open(f"{args.out_prefix}_hindex_fixed.tsv","w") as f:
        f.write(f"#Num_loci={keep_fixed.sum()}\n"); dff.to_csv(f,sep="\t",index=False)

    # ----- classification popmap -----
    groups=[popmap.get(s, infer_sim(s)) for s in class_samp]
    pd.DataFrame({"Sample":class_samp,"Group":groups}).to_csv(
        f"{args.out_prefix}_classification_popmap.tsv", sep="\t", index=False)

if __name__=="__main__":
    main()
