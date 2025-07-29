#!/usr/bin/env python3
"""
mask_samples.py
---------------
Outlier‐mask empirical individuals whose (HybridIndex, Heterozygosity)
fall outside the Mahalanobis ellipse of their simulated distribution
for their assigned NewHybrids category—excluding pure (P0/P1) classes.
Also outputs a detailed mask table and diagnostic plots per class
showing the 95% Mahalanobis ellipse.
"""
import sys
import argparse
from pathlib import Path

import pandas as pd
import numpy as np
from numpy.linalg import inv

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

# SciPy for chi2; fallback if missing
try:
    from scipy.stats import chi2
    def chi2_cutoff(alpha):
        return chi2.ppf(1 - alpha, df=2)
    def chi2_pvalue(d2):
        return chi2.sf(d2, df=2)
except ImportError:
    def chi2_cutoff(alpha):
        if abs(alpha - 0.05) < 1e-6:
            return 5.991
        raise RuntimeError("SciPy required for other alpha levels")
    def chi2_pvalue(d2):
        return np.nan

def read_hindex(path):
    df = pd.read_csv(path, sep='\t', comment='#')
    for c in ['Sample','HybridIndex','Heterozygosity']:
        if c not in df.columns:
            sys.exit(f"ERROR: Missing column '{c}' in {path}")
    return df

def read_nh(path):
    cats = ["P0","P1","F1","F2","Bx0","Bx1"]
    names = ["Index","Individual"] + cats
    return pd.read_csv(path, sep=r'\s+', skiprows=1, header=None, names=names)

def load_index_map(path):
    df = pd.read_csv(path, sep='\t', header=0)
    if 'Sample' in df.columns:
        df = df.rename(columns={'Sample':'Individual'})
    if not {'Index','Individual'}.issubset(df.columns):
        sys.exit("ERROR: index map must have 'Index' and 'Sample/Individual'")
    return df[['Index','Individual']]

def plot_diagnostics(cls, sim_pts, emp_pts, mu, cov, cutoff, prefix):
    fig, ax = plt.subplots(figsize=(6,6))
    ax.scatter(sim_pts[:,0], sim_pts[:,1], alpha=0.5, label='Simulation')
    ax.scatter(emp_pts[:,0], emp_pts[:,1], alpha=0.8, label='Empirical')
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:,order]
    angle = np.degrees(np.arctan2(*vecs[:,0][::-1]))
    width, height = 2 * np.sqrt(vals * chi2_cutoff(0.05))
    ellipse = Ellipse(
        xy=mu, width=width, height=height, angle=angle,
        edgecolor='black', facecolor='none', lw=1.5, label='95% Mahalanobis'
    )
    ax.add_patch(ellipse)
    ax.set_xlabel('HybridIndex')
    ax.set_ylabel('Heterozygosity')
    ax.set_title(f'Diagnostics: {cls}')
    ax.legend()
    out_png = f"{prefix}_{cls}_diag.png"
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    print(f"✅ Diagnostic plot: {out_png}")

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--hindex',         required=True,
                   help="TSV with Sample, HybridIndex, Heterozygosity")
    p.add_argument('--nh_results',     required=True,
                   help="NewHybrids posterior probabilities file")
    p.add_argument('--nh_index',       required=True,
                   help="TSV mapping Index→Individual")
    p.add_argument('--hindex_popmap',  required=True,
                   help="TSV with Sample→SimClass for simulated data")
    p.add_argument('--prob_threshold', type=float, required=True,
                   help="Minimum posterior probability to include empirical")
    p.add_argument('--out_prefix',     required=True,
                   help="Prefix for output files")
    p.add_argument('--alpha',          type=float, default=0.05,
                   help="Type I error rate for Mahalanobis cutoff")
    args = p.parse_args()

    # Read data
    df_h   = read_hindex(args.hindex)
    df_nh  = read_nh(args.nh_results)
    df_map = load_index_map(args.nh_index)

    # Merge NewHybrids assignments
    cats = ["P0","P1","F1","F2","Bx0","Bx1"]
    df_cls = df_nh.drop(columns='Individual').merge(df_map, on='Index')
    df_cls['MaxP']      = df_cls[cats].max(axis=1)
    df_cls['Assigned']  = df_cls[cats].idxmax(axis=1)
    df_emp = df_cls[df_cls['MaxP'] > args.prob_threshold][['Individual','Assigned']]
    df_emp = df_emp[~df_emp['Assigned'].isin(['P0','P1'])]
    if df_emp.empty:
        sys.exit("No non-pure empirical samples above threshold")

    # Load simulation mapping
    df_simmap = (
        pd.read_csv(args.hindex_popmap, sep='\t')
          .rename(columns={'Group':'SimClass'})
    )
    df_simmap = df_simmap[df_simmap['Sample'].str.startswith('SIM_')].copy()
    def map_cls(r):
        raw, s = r['SimClass'], r['Sample']
        if raw=='Pure': return 'P0' if '_P0_' in s else 'P1'
        if raw in ['F1','F2']: return raw
        if raw.startswith('BC-P0'): return 'Bx0'
        if raw.startswith('BC-P1'): return 'Bx1'
        return raw
    df_simmap['SimClass'] = df_simmap.apply(map_cls, axis=1)
    df_simmap = df_simmap[df_simmap['SimClass'].isin(cats)]

    # Merge simulation points
    df_sim = df_h.merge(
        df_simmap[['Sample','SimClass']],
        on='Sample'
    )

    # Prepare output records
    results = []
    cutoff  = chi2_cutoff(args.alpha)

    for cls in sorted(df_emp['Assigned'].unique()):
        sim_sub = df_sim[df_sim['SimClass'] == cls]
        sim_pts = sim_sub[['HybridIndex','Heterozygosity']].values
        if len(sim_pts) < 3:
            continue

        # Compute ellipse parameters
        mu     = sim_pts.mean(axis=0)
        cov    = np.cov(sim_pts, rowvar=False)
        invcov = inv(cov)

        # Empirical points for this class
        emp_sub = df_emp[df_emp['Assigned'] == cls]
        emp_data = emp_sub.merge(df_h, left_on='Individual', right_on='Sample')

        emp_pts = []
        for _, r in emp_data.iterrows():
            x  = np.array([r['HybridIndex'], r['Heterozygosity']])
            d2 = float((x - mu) @ invcov @ (x - mu))
            pval = chi2_pvalue(d2)
            inside = d2 <= cutoff
            results.append({
                'Sample': r['Individual'],
                'Class': cls,
                'HybridIndex': r['HybridIndex'],
                'Heterozygosity': r['Heterozygosity'],
                'MahalanobisD2': d2,
                'p_value': pval,
                'In_Mahalanobis': inside
            })
            emp_pts.append(x)

        # Plot diagnostics ellipse
        plot_diagnostics(cls, sim_pts, np.vstack(emp_pts), mu, cov, cutoff, args.out_prefix)

    # Save outputs
    df_res    = pd.DataFrame(results)
    masked    = df_res[~df_res['In_Mahalanobis']]['Sample'].tolist()
    mask_file = f"{args.out_prefix}_masked_samples.txt"
    Path(mask_file).write_text("\n".join(masked))
    print(f"✅ Masked {len(masked)} samples → {mask_file}")

    table_file = f"{args.out_prefix}_mask_table.tsv"
    df_res.to_csv(table_file, sep='\t', index=False)
    print(f"✅ Detailed table → {table_file}")

if __name__ == '__main__':
    main()
