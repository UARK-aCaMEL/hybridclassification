#!/usr/bin/env python3
"""
plot_triangle.py
────────────────
Interactive side‑by‑side triangle plots (all loci vs fixed differences)
with Show/Hide toggle for simulated individuals.
"""

from __future__ import annotations
import argparse, re
from pathlib import Path
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px

def load_hindex(path: Path):
    with open(path) as fh:
        first = fh.readline()
    m = re.match(r"#Num_loci=(\d+)", first.strip())
    n_loci = int(m.group(1)) if m else np.nan
    return n_loci, pd.read_csv(path, sep="\t", comment="#")

def load_popmap(path: Path):
    pm = pd.read_csv(path, sep="\t", header=None,
                     names=["Sample", "Pop"], dtype=str)
    return dict(zip(pm.Sample, pm.Pop))

def triangle_outline():
    x = np.linspace(0, 1, 200)
    return [
        go.Scatter(
            x=[0, .5, 1, 0], y=[0, 1, 0, 0],
            mode="lines", line=dict(color="black"),
            showlegend=False, legendgroup="outline"
        ),
        go.Scatter(
            x=x, y=2 * x * (1 - x),
            mode="lines", line=dict(color="black"),
            showlegend=False, legendgroup="outline"
        )
    ]

def make_scatter(sub: pd.DataFrame, name: str, color: str,
                 legend_grp: str, showlegend: bool, visible: bool):
    return go.Scatter(
        x=sub.HybridIndex, y=sub.Heterozygosity,
        mode="markers",
        marker=dict(color=color, size=6, line=dict(width=.3, color="black")),
        name=name, legendgroup=legend_grp, showlegend=showlegend,
        visible=visible,
        customdata=sub[["Sample","Pop","PercMissing"]].to_numpy(),
        hovertemplate=(
            "Sample=%{customdata[0]}<br>"
            "Pop=%{customdata[1]}<br>"
            "Missing=%{customdata[2]:.1%}<br>"
            "h=%{x:.2f}, H=%{y:.2f}<extra></extra>"
        )
    )

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--result",        required=True)
    ap.add_argument("--result_fixed",  required=True)
    ap.add_argument("--popmap",        required=True)
    ap.add_argument("--triangle_map",  required=True)
    ap.add_argument("--template")
    ap.add_argument("--palette",       default="Set2")
    ap.add_argument("--out",           required=True)
    args = ap.parse_args()

    n_all,   df_all   = load_hindex(Path(args.result))
    n_fixed, df_fixed = load_hindex(Path(args.result_fixed))

    popmap = load_popmap(Path(args.popmap))
    for df in (df_all, df_fixed):
        df["Pop"] = df.Sample.map(popmap)

    tmap        = pd.read_csv(args.triangle_map, sep="\t", dtype=str)
    group_of    = tmap.set_index("Sample")["Group"].to_dict()
    sim_samples = {s for s,g in group_of.items() if g != "ADMIX"}
    emp_samples = {s for s,g in group_of.items() if g == "ADMIX"}

    pops    = sorted(df_all[df_all.Sample.isin(emp_samples)]["Pop"].dropna().unique())
    classes = sorted({group_of[s] for s in sim_samples})
    palette = px.colors.qualitative.__dict__.get(args.palette,
                                                 px.colors.qualitative.Set2)
    pop_col = {p: palette[i % len(palette)] for i,p in enumerate(pops)}
    cls_col = {c: palette[i % len(palette)] for i,c in enumerate(classes)}
    grey    = "lightgrey"

    fig = make_subplots(
        cols=2, rows=1, shared_yaxes=True, horizontal_spacing=.07,
        subplot_titles=(f"All loci (N = {n_all})",
                        f"Fixed differences (N = {n_fixed})")
    )

    # single list of (trace, row, col, visible_default, visible_sim)
    traces: list[tuple[go.Trace,int,int,bool,bool]] = []

    # 1) Triangle outlines (visible in both)
    for col in (1,2):
        for outline in triangle_outline():
            traces.append((outline, 1, col, True, True))

    # 2) Default empirical, colored by pop
    for pop, colr in pop_col.items():
        for data, col in ((df_all,1),(df_fixed,2)):
            sub = data[(data.Sample.isin(emp_samples)) & (data.Pop == pop)]
            # showlegend only on first subplot, legendgroup=pop
            tr = make_scatter(sub, pop, colr, legend_grp=pop,
                              showlegend=(col==1), visible=(True if col else False))
            traces.append((tr, 1, col, True, False))

    # 3) Simulated view: empirical in grey
    for data, col in ((df_all,1), (df_fixed,2)):
        sub = data[data.Sample.isin(emp_samples)]
        tr = make_scatter(sub, "Empirical samples", grey,
                          legend_grp="Empirical samples",
                          showlegend=(col==1), visible=False)
        traces.append((tr, 1, col, False, True))

    # 4) Simulated: classes in color
    for cls, colr in cls_col.items():
        for data, col in ((df_all,1),(df_fixed,2)):
            sub = data[(data.Sample.isin(sim_samples)) &
                       (data.Sample.map(group_of) == cls)]
            tr = make_scatter(sub, cls, colr,
                              legend_grp=cls,
                              showlegend=(col==1), visible=False)
            traces.append((tr, 1, col, False, True))

    # add traces with initial visibility = visible_default
    for tr, row, col, vd, vs in traces:
        tr.visible = vd
        fig.add_trace(tr, row=row, col=col)

    # build masks
    vis_default = [vd for _,_,_,vd,_ in traces]
    vis_sim     = [vs for _,_,_,_,vs in traces]

    # axes
    for c in (1,2):
        fig.update_xaxes(range=[-0.05,1.05], title="Hybrid Index",
                         row=1, col=c)
    fig.update_yaxes(range=[-0.05,1.05],
                     title="Inter‑class Heterozygosity",
                     row=1, col=1)

    # toggle buttons
    fig.update_layout(
        template="simple_white", height=500,
        legend_title="Group", margin=dict(t=80,b=40,l=40,r=40),
        updatemenus=[dict(
            type="buttons", direction="right",
            x=0.5, xanchor="center", y=1.20, yanchor="top",
            buttons=[
                dict(label="Hide simulated individuals",
                     method="update",
                     args=[{"visible": vis_default}]),
                dict(label="Show simulated individuals",
                     method="update",
                     args=[{"visible": vis_sim}])
            ]
        )]
    )

    html = fig.to_html(full_html=False, include_plotlyjs="cdn")
    out = Path(args.out)
    if args.template and Path(args.template).exists():
        header = Path(args.template).read_text().rstrip() + "\n"
        out.write_text(header + html)
    else:
        out.write_text(html)

    print(f"✅ Triangle plot saved → {args.out}")

if __name__ == "__main__":
    main()
