#!/usr/bin/env python3
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import argparse
import re
from pathlib import Path

def read_results(path):
    cats = ["P0","P1","F1","F2","Bx0","Bx1"]
    cols = ["Index","Individual"] + cats
    return pd.read_csv(
        path, sep=r'\s+', skiprows=1, header=None,
        names=cols, dtype={"Index": int}
    )

def read_map(path):
    df = pd.read_csv(path, sep="\t", header=0)
    df.columns = [c.strip() for c in df.columns]
    df = df.rename(
        columns={df.columns[0]: "Index", df.columns[1]: "Individual"}
    )
    df["Index"] = df["Index"].astype(int)
    return df[["Index","Individual"]]

def assemble(df_res, df_map):
    df_res["Index"] = df_res["Index"].astype(int)
    return df_res.drop(columns="Individual") \
                 .merge(df_map, on="Index")

def prepare_sim(df):
    pattern = re.compile(
        r'^(?:SIM_)?(\d+)_(?:Pure_)?'
        r'(P0|P1|F1|F2|BC-P0|BC-P1)_\d+'
    )
    df_sim = df[df["Individual"].str.match(pattern)].copy()

    def parse(x):
        m = pattern.match(x)
        rep, cat = m.group(1), m.group(2)
        if cat == "BC-P0": cat = "Bx0"
        if cat == "BC-P1": cat = "Bx1"
        return int(rep), cat

    df_sim[["Replicate","TrueCategory"]] = (
        df_sim["Individual"]
              .apply(lambda x: pd.Series(parse(x)))
    )
    cats = ["P0","P1","F1","F2","Bx0","Bx1"]
    df_sim["Assigned"]     = df_sim[cats].idxmax(axis=1)
    df_sim["MaxPosterior"] = df_sim[cats].max(axis=1)
    return df_sim

def compute_accuracy(df_sim):
    thresholds = np.round(np.linspace(0.5,1.0,51),2)
    rec = []
    for thr in thresholds:
        for (rep, cat), grp in df_sim.groupby(["Replicate","TrueCategory"]):
            n = len(grp)
            if n == 0: continue
            correct = ((grp["Assigned"] == cat) &
                       (grp[cat] >= thr)).sum()
            rec.append({
                "Replicate": rep,
                "TrueCategory": cat,
                "Threshold": thr,
                "Accuracy": correct / n
            })
    return pd.DataFrame(rec)

def compute_rates(df_sim):
    thresholds = np.round(np.linspace(0.5,1.0,51),2)
    cats = ["P0","P1","F1","F2","Bx0","Bx1"]
    rec = []
    for thr in thresholds:
        for rep, grp in df_sim.groupby("Replicate"):
            above = (grp[cats] >= thr)
            cnt_above = above.sum(axis=1)
            calls = np.where(
                cnt_above == 1,
                grp[cats].idxmax(axis=1),
                "Unclassified"
            )
            # misclassified = assigned to exactly one cat ≠ TrueCategory
            mis = ((calls != grp["TrueCategory"]) &
                   (calls != "Unclassified")).sum()
            mis_rate = mis / len(grp)
            # unclassified = zero or multiple above-threshold
            unclass = (calls == "Unclassified").sum()
            unclass_rate = unclass / len(grp)
            rec.append({
                "Replicate": rep,
                "Threshold": thr,
                "MisRate": mis_rate,
                "UnclassRate": unclass_rate
            })
    return pd.DataFrame(rec)

def choose_colors(n, palette):
    import plotly.express as px
    if hasattr(px.colors.qualitative, palette):
        return getattr(px.colors.qualitative, palette)
    if hasattr(px.colors.sequential, palette):
        scale = getattr(px.colors.sequential, palette)
        return px.colors.sample_colorscale(
            scale, [i/(n-1) for i in range(n)]
        )
    if hasattr(px.colors.diverging, palette):
        scale = getattr(px.colors.diverging, palette)
        return px.colors.sample_colorscale(
            scale, [i/(n-1) for i in range(n)]
        )
    raise ValueError(f"Palette '{palette}' not found")

def make_plot(acc_df, df_sim, out_html,
              template=None, palette="Viridis",
              threshold=None):

    # 1) per‐category accuracy ribbons + lines
    stats = acc_df.groupby(
        ["Threshold","TrueCategory"]
    )["Accuracy"].agg(mean="mean", std="std").reset_index()
    categories = ["P0","Bx0","F1","F2","Bx1","P1"]
    colors = choose_colors(len(categories), palette)

    fig = go.Figure()
    for i, cat in enumerate(categories):
        sub = stats[stats["TrueCategory"] == cat]
        t, m, s = sub["Threshold"], sub["mean"], np.nan_to_num(sub["std"])
        fig.add_trace(go.Scatter(
            x=np.concatenate([t, t[::-1]]),
            y=np.concatenate([m+s, (m-s)[::-1]]),
            fill="toself",
            fillcolor=colors[i]
                      .replace("rgb","rgba")
                      .replace(")",",0.2)"),
            line_color="rgba(0,0,0,0)",
            showlegend=False, hoverinfo="skip"
        ))
        fig.add_trace(go.Scatter(
            x=t, y=m, mode="lines", name=cat,
            line=dict(color=colors[i]),
            hovertemplate="Threshold=%{x:.2f}<br>Accuracy=%{y:.2f}"
        ))

    # 2) misclassified & unclassified mean±SD
    rates = compute_rates(df_sim)
    rate_stats = rates.groupby("Threshold").agg(
        MisMean=('MisRate','mean'),
        MisStd =('MisRate','std'),
        UncMean=('UnclassRate','mean'),
        UncStd =('UnclassRate','std')
    ).reset_index()

    # misclassified (black)
    t, m, s = rate_stats["Threshold"], rate_stats["MisMean"], rate_stats["MisStd"]
    fig.add_trace(go.Scatter(
        x=np.concatenate([t, t[::-1]]),
        y=np.concatenate([m+s, (m-s)[::-1]]),
        fill="toself",
        fillcolor="rgba(0,0,0,0.2)",
        line_color="rgba(0,0,0,0)",
        showlegend=False, hoverinfo="skip"
    ))
    fig.add_trace(go.Scatter(
        x=t, y=m, mode="lines", name="Misclassified",
        line=dict(color="black", dash="solid"),
        hovertemplate="Threshold=%{x:.2f}<br>Misclassified=%{y:.2f}"
    ))

    # unclassified (gray)
    t2, m2, s2 = rate_stats["Threshold"], rate_stats["UncMean"], rate_stats["UncStd"]
    fig.add_trace(go.Scatter(
        x=np.concatenate([t2, t2[::-1]]),
        y=np.concatenate([m2+s2, (m2-s2)[::-1]]),
        fill="toself",
        fillcolor="rgba(128,128,128,0.2)",
        line_color="rgba(0,0,0,0)",
        showlegend=False, hoverinfo="skip"
    ))
    fig.add_trace(go.Scatter(
        x=t2, y=m2, mode="lines", name="Unclassified",
        line=dict(color="gray", dash="dash"),
        hovertemplate="Threshold=%{x:.2f}<br>Unclassified=%{y:.2f}"
    ))

    # 3) layout + empirical threshold marker
    fig.update_layout(
        xaxis_title="Posterior Probability Threshold",
        yaxis_title="Accuracy / Rate",
        xaxis=dict(range=[0.5,1.0], tickmode="linear",
                   tick0=0.5, dtick=0.05),
        yaxis=dict(range=[0,1], tickformat=".2f"),
        legend_title="Category / Rate",
        template="simple_white",
        margin=dict(t=100, b=60)
    )

    if threshold is not None:
        fig.add_shape(
            type="line", x0=threshold, x1=threshold,
            y0=0, y1=1, xref="x", yref="paper",
            line=dict(color="red", dash="dash", width=2),
            layer="above"
        )
        fig.add_annotation(
            x=threshold, y=1.02, xref="x", yref="paper",
            text=f"Empirical threshold = {threshold:.2f}",
            showarrow=False, font=dict(color="red"), align="center"
        )

    body = fig.to_html(full_html=False, include_plotlyjs="cdn")
    if template and Path(template).exists():
        header = Path(template).read_text().rstrip() + "\n"
        Path(out_html).write_text(header + body)
    else:
        Path(out_html).write_text(body)

    print(f"✅ Plot saved to: {out_html}")

if __name__ == "__main__":
    p = argparse.ArgumentParser(
        description="Plot simulation assignment accuracy vs threshold."
    )
    p.add_argument("--result",   required=True,
                   help="posterior file")
    p.add_argument("--map",      required=True,
                   help="Index→sample TSV")
    p.add_argument("--template", help="HTML header to prepend")
    p.add_argument("--palette",  default="Spectral",
                   help="plotly.colors palette")
    p.add_argument("--threshold", type=float,
                   help="empirical posterior probability threshold")
    p.add_argument("--out",      required=True,
                   help="Output HTML file")
    args = p.parse_args()

    df_res = read_results(args.result)
    df_map = read_map(args.map)
    df     = assemble(df_res, df_map)
    df_sim = prepare_sim(df)
    acc_df = compute_accuracy(df_sim)

    make_plot(
        acc_df, df_sim,
        args.out, template=args.template,
        palette=args.palette,
        threshold=args.threshold
    )
