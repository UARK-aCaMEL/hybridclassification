#!/usr/bin/env python3
"""
mask_samples.py
---------------
Outlier‐mask empirical individuals whose (HybridIndex,Heterozygosity)
fall outside the Mahalanobis ellipse of their simulated distribution
for their assigned NewHybrids category—excluding pure (P0/P1) classes.
Also outputs a detailed mask table and diagnostic plots per class,
including KDE density surface plots with overlaid points.
"""
import sys
import argparse
from pathlib import Path

import pandas as pd
import numpy as np
from numpy.linalg import inv

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from sklearn.neighbors import KernelDensity

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

# Utility functions

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


def plot_diagnostics(cls, sim_pts, emp_pts, mu, cov, cutoff, kde, kde_cut, prefix):
    fig, ax = plt.subplots(figsize=(6,6))
    ax.scatter(sim_pts[:,0], sim_pts[:,1], alpha=0.5, label='Simulation')
    ax.scatter(emp_pts[:,0], emp_pts[:,1], alpha=0.8, label='Empirical')
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:,order]
    angle = np.degrees(np.arctan2(*vecs[:,0][::-1]))
    width, height = 2 * np.sqrt(vals * chi2_cutoff(0.05))
    ellipse = Ellipse(xy=mu, width=width, height=height, angle=angle,
                      edgecolor='black', facecolor='none', lw=1.5, label='95% Mahalanobis')
    ax.add_patch(ellipse)
    xx, yy, zz = _compute_kde_grid(sim_pts, kde)
    ax.contour(xx, yy, zz, levels=[kde_cut], linestyles='--', colors='green')
    ax.set_xlabel('HybridIndex')
    ax.set_ylabel('Heterozygosity')
    ax.set_title(f'Diagnostics: {cls}')
    ax.legend()
    out_png = f"{prefix}_{cls}_diag.png"
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    print(f"✅ Diagnostic plot: {out_png}")


def plot_kde_surface(cls, sim_pts, emp_pts, kde, kde_cut, prefix):
    xx, yy, zz = _compute_kde_grid(sim_pts, kde)
    fig, ax = plt.subplots(figsize=(6,6))
    cf = ax.contourf(xx, yy, zz, levels=20, cmap='viridis')
    ax.contour(xx, yy, zz, levels=[kde_cut], colors='white', linestyles='--')
    # overlay points
    ax.scatter(sim_pts[:,0], sim_pts[:,1], c='blue', alpha=0.4, label='Simulation')
    ax.scatter(emp_pts[:,0], emp_pts[:,1], c='red', alpha=0.8, label='Empirical')
    cbar = fig.colorbar(cf, ax=ax, label='Density')
    ax.set_xlabel('HybridIndex')
    ax.set_ylabel('Heterozygosity')
    ax.set_title(f'KDE density surface: {cls}')
    ax.legend()
    out_png = f"{prefix}_{cls}_kde.png"
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    print(f"✅ KDE surface plot: {out_png}")


def _compute_kde_grid(sim_pts, kde):
    x_min, x_max = sim_pts[:,0].min(), sim_pts[:,0].max()
    y_min, y_max = sim_pts[:,1].min(), sim_pts[:,1].max()
    xx, yy = np.mgrid[x_min:x_max:100j, y_min:y_max:100j]
    grid = np.vstack([xx.ravel(), yy.ravel()]).T
    zz = np.exp(kde.score_samples(grid)).reshape(xx.shape)
    return xx, yy, zz


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--hindex',        required=True)
    p.add_argument('--nh_results',    required=True)
    p.add_argument('--nh_index',      required=True)
    p.add_argument('--hindex_popmap', required=True)
    p.add_argument('--prob_threshold',type=float, required=True)
    p.add_argument('--out_prefix',    required=True)
    p.add_argument('--alpha',         type=float, default=0.05)
    p.add_argument('--kde_percentile',type=float, default=5.0)
    args = p.parse_args()

    df_h = read_hindex(args.hindex)
    df_nh = read_nh(args.nh_results)
    df_map = load_index_map(args.nh_index)
    df_cls = df_nh.drop(columns='Individual').merge(df_map, on='Index')
    cats = ["P0","P1","F1","F2","Bx0","Bx1"]
    df_cls['MaxP'] = df_cls[cats].max(axis=1)
    df_cls['Assigned'] = df_cls[cats].idxmax(axis=1)
    df_emp = df_cls[df_cls['MaxP']>args.prob_threshold][['Individual','Assigned']]
    df_emp = df_emp[~df_emp['Assigned'].isin(['P0','P1'])]
    if df_emp.empty:
        sys.exit("No non-pure empirical samples above threshold")
    df_simmap = pd.read_csv(args.hindex_popmap, sep='\t').rename(columns={'Group':'SimClass'})
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
    df_sim = df_h.merge(df_simmap[['Sample','SimClass']], on='Sample')

    results = []
    cutoff = chi2_cutoff(args.alpha)

    for cls in sorted(df_emp['Assigned'].unique()):
        sim_sub = df_sim[df_sim['SimClass']==cls]
        sim_pts = sim_sub[['HybridIndex','Heterozygosity']].values
        if len(sim_pts)<3: continue
        mu = sim_pts.mean(0)
        cov = np.cov(sim_pts, rowvar=False)
        invcov = inv(cov)
        kde = KernelDensity().fit(sim_pts)
        logdens = kde.score_samples(sim_pts)
        kde_cut = np.percentile(logdens, args.kde_percentile)
        emp_sub = df_emp[df_emp['Assigned']==cls]
        emp_data = emp_sub.merge(df_h, left_on='Individual', right_on='Sample')
        emp_pts = []
        for _,r in emp_data.iterrows():
            x = np.array([r['HybridIndex'], r['Heterozygosity']])
            d2 = float((x-mu)@invcov@(x-mu))
            pval = chi2_pvalue(d2)
            logd = float(kde.score_samples(x.reshape(1,-1)))
            in_kde = logd>=kde_cut
            out = d2>cutoff
            results.append({
                'Sample':r['Individual'],'Class':cls,
                'HybridIndex':r['HybridIndex'],'Heterozygosity':r['Heterozygosity'],
                'MahalanobisD2':d2,'p_value':pval,'In_Mahalanobis':out,
                'KDE_log_density':logd,'In_KDE':in_kde
            })
            emp_pts.append(x)
        plot_diagnostics(cls, sim_pts, np.vstack(emp_pts), mu, cov,
                         cutoff, kde, kde_cut, args.out_prefix)
        plot_kde_surface(cls, sim_pts, np.vstack(emp_pts), kde, kde_cut, args.out_prefix)

    df_res = pd.DataFrame(results)
    mask_list = df_res[df_res['In_Mahalanobis']]['Sample'].tolist()
    Path(f"{args.out_prefix}_masked_samples.txt").write_text("\n".join(mask_list))
    print(f"✅ Masked {len(mask_list)} samples → {args.out_prefix}_masked_samples.txt")
    df_res.to_csv(f"{args.out_prefix}_mask_table.tsv", sep='\t', index=False)
    print(f"✅ Detailed table → {args.out_prefix}_mask_table.tsv")

if __name__=='__main__':
    main()