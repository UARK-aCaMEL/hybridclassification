#!/usr/bin/env python3
"""
plot_nh_outliers.py
------------------
Creates a Plotly HTML with one subplot per hybrid class, showing:
  - Simulation points (HybridIndex vs Heterozygosity)
  - Empirical points (with NewHybrids posteriors in hover)
  - 95% Mahalanobis ellipse
Usage:
    plot_nh_outliers.py \
        --hindex         <hindex.tsv> \
        --nh_results     <NewHybrids posterior file> \
        --nh_index       <Index→sample TSV> \
        --hindex_popmap  <classification popmap TSV> \
        --prob_threshold <min posterior to assign empirical> \
        --alpha          <outlier α> \
        --template       <HTML header file> \
        --out            <output HTML>
"""
import sys
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import chi2
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# --- data loading ---
def read_hindex(path):
    df = pd.read_csv(path, sep="\t", comment="#")
    return df[['Sample','HybridIndex','Heterozygosity']]

def read_nh(path):
    cats = ["P0","P1","F1","F2","Bx0","Bx1"]
    names = ['Index','Individual'] + cats
    return pd.read_csv(path, sep=r"\s+", skiprows=1, header=None, names=names)

def load_index_map(path):
    df = pd.read_csv(path, sep="\t", header=0)
    if 'Individual' in df.columns and 'Sample' not in df.columns:
        df = df.rename(columns={'Individual':'Sample'})
    return df[['Index','Sample']]

def load_simmap(path):
    df = pd.read_csv(path, sep="\t", header=None, names=['Sample','Group'])
    return df[df['Sample'].str.startswith('SIM_')].copy()

# --- ellipse calculation ---
def ellipse_coords(mu, cov, alpha, ns=100):
    vals, vecs = np.linalg.eigh(cov)
    order = np.argsort(vals)[::-1]
    vals, vecs = vals[order], vecs[:,order]
    r2 = chi2.ppf(1-alpha, df=2)
    theta = np.linspace(0, 2*np.pi, ns)
    circle = np.vstack([np.cos(theta), np.sin(theta)])
    axes = np.sqrt(vals * r2)
    transform = vecs @ np.diag(axes)
    pts = (transform @ circle).T + mu
    return pts[:,0], pts[:,1]

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--hindex',         required=True)
    p.add_argument('--nh_results',     required=True)
    p.add_argument('--nh_index',       required=True)
    p.add_argument('--hindex_popmap',  required=True)
    p.add_argument('--prob_threshold', type=float, required=True)
    p.add_argument('--alpha',          type=float, default=0.05)
    p.add_argument('--template',       required=True)
    p.add_argument('--out',            required=True)
    args = p.parse_args()

    # 1) load H‐index (coords)
    df_h   = read_hindex(args.hindex)         # Sample, HybridIndex, Heterozygosity

    # 2) load NewHybrids posteriors + map Index→Sample
    df_nh  = read_nh(args.nh_results)         # Index, Individual, P0…Bx1
    df_map = load_index_map(args.nh_index)    # Index, Sample

    # join posterior → Sample, then bring in coords
    df_nh2  = df_nh.drop(columns=['Individual'])
    df_hcls = df_nh2.merge(df_map, on='Index')\
                    .merge(df_h, on='Sample')
    cats = ['P0','P1','F1','F2','Bx0','Bx1']
    df_hcls['MaxP']     = df_hcls[cats].max(axis=1)
    df_hcls['Assigned'] = df_hcls[cats].idxmax(axis=1)

    # 3) empirical hybrids above threshold
    df_emp = df_hcls[df_hcls['MaxP'] > args.prob_threshold]
    df_emp = df_emp[~df_emp['Assigned'].isin(['P0','P1'])]
    if df_emp.empty:
        sys.exit("No empirical hybrids above threshold")

    # 4) load & classify simulation samples
    simmap = load_simmap(args.hindex_popmap)
    def map_group(samp, grp):
        if grp=='Pure':       return 'P0' if '_P0_' in samp else 'P1'
        if grp in ['F1','F2']: return grp
        if grp.startswith('BC-P0'): return 'Bx0'
        if grp.startswith('BC-P1'): return 'Bx1'
        return grp
    simmap['Class'] = [map_group(s,g) for s,g in zip(simmap['Sample'], simmap['Group'])]
    # join coords → Sample,Class
    df_sim = df_h.merge(simmap[['Sample','Class']], on='Sample')

    # 5) layout subplots
    classes = sorted(df_emp['Assigned'].unique())
    n      = len(classes)
    cols   = n if n < 4 else 2
    rows   = int(np.ceil(n/cols))
    fig    = make_subplots(
                rows=rows, cols=cols,
                subplot_titles=classes,
                horizontal_spacing=0.1,
                vertical_spacing=0.2
             )

    # hover templates
    sim_hover = "Sample: %{text}<extra></extra>"
    emp_hover = (
        "Sample: %{text}<br>"
        "P0: %{customdata[0]:.2f}<br>"
        "P1: %{customdata[1]:.2f}<br>"
        "F1: %{customdata[2]:.2f}<br>"
        "F2: %{customdata[3]:.2f}<br>"
        "Bx0: %{customdata[4]:.2f}<br>"
        "Bx1: %{customdata[5]:.2f}<extra></extra>"
    )

    # 6) plot each class
    for idx, cls in enumerate(classes):
        r = idx // cols + 1
        c = idx % cols + 1

        sim_df = df_sim[df_sim['Class'] == cls]
        emp_df = df_emp[df_emp['Assigned'] == cls]

        # ellipse on sim points
        pts = sim_df[['HybridIndex','Heterozygosity']].values
        mu, cov = pts.mean(axis=0), np.cov(pts, rowvar=False)
        xe, ye  = ellipse_coords(mu, cov, args.alpha)

        # simulation scatter
        fig.add_trace(go.Scatter(
            x=sim_df['HybridIndex'], y=sim_df['Heterozygosity'],
            mode='markers',
            marker=dict(color='blue', opacity=0.5),
            name='Simulation', showlegend=(idx==0),
            text=sim_df['Sample'],
            hovertemplate=sim_hover
        ), row=r, col=c)

        # empirical scatter (with posteriors)
        fig.add_trace(go.Scatter(
            x=emp_df['HybridIndex'], y=emp_df['Heterozygosity'],
            mode='markers',
            marker=dict(color='red', symbol='x', size=8),
            name='Empirical', showlegend=(idx==0),
            text=emp_df['Sample'],
            customdata=emp_df[cats].values,
            hovertemplate=emp_hover
        ), row=r, col=c)

        # ellipse line
        fig.add_trace(go.Scatter(
            x=xe, y=ye, mode='lines',
            line=dict(color='black'),
            name='95% ellipse', showlegend=(idx==0)
        ), row=r, col=c)

        fig.update_xaxes(title_text='HybridIndex', row=r, col=c)
        fig.update_yaxes(title_text='Heterozygosity', row=r, col=c)

    # 7) finalize & write
    fig.update_layout(
        height=300*rows, width=400*cols,
        title_text='NH Outlier Diagnostics',
        template='simple_white',
        margin=dict(t=100, b=80)
    )
    header = Path(args.template).read_text()
    body   = fig.to_html(full_html=True, include_plotlyjs='cdn')
    Path(args.out).write_text(header + body)
    print(f"✅ Plot saved to {args.out}")

if __name__=='__main__':
    main()
