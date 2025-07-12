#!/usr/bin/env python3
import pandas as pd
import plotly.express as px
import argparse
from pathlib import Path

def load_data(qmat_file, ind_file, pop_file):
    # read Q matrix
    q_raw = pd.read_csv(qmat_file, sep=":", header=None)
    q_df  = q_raw[1].str.strip().str.split(expand=True).astype(float)
    # read metadata
    inds = pd.read_csv(ind_file, header=None, names=["Individual"])
    pops = pd.read_csv(pop_file, header=None, names=["Population"])
    if not (len(inds) == len(pops) == len(q_df)):
        raise ValueError("input length mismatch")
    q_df.columns       = [f"Cluster {i+1}" for i in range(q_df.shape[1])]
    q_df["Individual"] = inds["Individual"]
    q_df["Population"] = pops["Population"]
    return q_df

def make_plot(df, out_html, template=None,
              palette="Spectral", sort_pop=False,
              sort_ind=False, order_cluster=None):
    cluster_cols = [c for c in df.columns if c.startswith("Cluster")]

    # determine ordering
    if order_cluster:
        cname = f"Cluster {order_cluster}"
        individual_order = df.sort_values(cname, ascending=False)["Individual"].tolist()
    else:
        if sort_pop:
            pop_means   = df.groupby("Population")[cluster_cols].mean()
            pop_dom     = pop_means.idxmax(axis=1)
            pops_sorted = sorted(pop_dom.index, key=lambda p: int(pop_dom[p].split()[1]))
        else:
            pops_sorted = df["Population"].unique().tolist()
        if sort_ind:
            doms = df.groupby("Population")[cluster_cols].mean().idxmax(axis=1)
            individual_order = []
            for pop in pops_sorted:
                sub = df[df.Population==pop].sort_values(doms[pop], ascending=False)
                individual_order += sub["Individual"].tolist()
        else:
            seen, individual_order = set(), []
            for i in df["Individual"]:
                if i not in seen:
                    seen.add(i); individual_order.append(i)

    # compute population ticks
    pop_map  = df[["Individual","Population"]].drop_duplicates().set_index("Individual")
    pops     = pop_map.loc[individual_order]["Population"]
    counts   = pops.value_counts().loc[pops.unique()]
    positions= counts.cumsum() - counts/2

    # melt for plotting
    id_vars = ["Individual","Population"]
    if "Prior Classification" in df.columns:
        id_vars.append("Prior Classification")
    df_long = df.melt(id_vars=id_vars, var_name="Cluster", value_name="Proportion")
    df_long["Individual"] = pd.Categorical(df_long["Individual"],
                                           categories=individual_order,
                                           ordered=True)

    # color sequence
    n = len(cluster_cols)
    if hasattr(px.colors.qualitative, palette):
        seq = getattr(px.colors.qualitative, palette)
    elif hasattr(px.colors.sequential, palette):
        seq = px.colors.sample_colorscale(getattr(px.colors.sequential, palette),
                   [i/(n-1) for i in range(n)])
    elif hasattr(px.colors.diverging, palette):
        seq = px.colors.sample_colorscale(getattr(px.colors.diverging, palette),
                   [i/(n-1) for i in range(n)])
    else:
        raise ValueError(f"palette '{palette}' not found")

    # build hover fields
    hover = ["Individual","Population","Cluster","Proportion"]
    if "Prior Classification" in df.columns:
        hover.append("Prior Classification")

    # generate bar plot
    fig = px.bar(
        df_long,
        x="Individual", y="Proportion", color="Cluster",
        category_orders={"Individual": individual_order},
        color_discrete_sequence=seq,
        hover_data=hover
    )
    fig.update_traces(marker_line_width=0)
    fig.update_layout(
        barmode="stack",
        yaxis=dict(title="Ancestry Proportion", range=[0,1], showgrid=False),
        xaxis=dict(
            tickmode="array",
            tickvals=positions.values,
            ticktext=positions.index,
            title="Population",
            showgrid=False
        ),
        margin=dict(t=60,b=60),
        title=dict(text="ADMIXTURE Ancestry Barplot", x=0.5),
        legend_title="Cluster",
        template="simple_white"
    )

    # save HTML
    html_body = fig.to_html(full_html=False, include_plotlyjs="cdn")
    if template:
        header = Path(template).read_text()
        Path(out_html).write_text(header + html_body)
    else:
        Path(out_html).write_text(html_body)

    print(f"Plot saved to {out_html}")

if __name__=="__main__":
    p = argparse.ArgumentParser(description="ADMIXTURE barplot")
    p.add_argument("--clumpp",     required=True, help="colon‐sep Q matrix")
    p.add_argument("--inds",       required=True, help="sample IDs")
    p.add_argument("--pops",       required=True, help="population IDs")
    p.add_argument("--candidates", help="TSV Individual<TAB>Prior Classification")
    p.add_argument("--out",        required=True, help="output HTML")
    p.add_argument("--template",   help="HTML header to prepend")
    p.add_argument("--palette",    default="Spectral", help="plotly.colors palette")
    p.add_argument("--sort_pop",   action="store_true")
    p.add_argument("--sort_ind",   action="store_true")
    p.add_argument("--order_cluster", type=int,
                   help="1-based cluster to sort by")
    args = p.parse_args()

    df = load_data(args.clumpp, args.inds, args.pops)
    if args.candidates:
        cand = pd.read_csv(args.candidates, sep="\t", header=None,
                           names=["Individual","Prior Classification"])
        df = df.merge(cand, on="Individual", how="left")
    make_plot(
        df, args.out,
        template=args.template,
        palette=args.palette,
        sort_pop=args.sort_pop,
        sort_ind=args.sort_ind,
        order_cluster=args.order_cluster
    )
