#!/usr/bin/env python3
import pandas as pd
import plotly.express as px
import argparse
from pathlib import Path

def read_nh_results(path):
    cats = ["P0","P1","F1","F2","Bx0","Bx1"]
    cols = ["Index","Individual"] + cats
    return pd.read_csv(path, sep=r'\s+', skiprows=1, header=None, names=cols)

def load_maps(nh_map, popmap, speciesmap):
    df_map = pd.read_csv(nh_map, sep="\t", header=0).rename(columns={"Sample":"Individual"})
    df_pop = pd.read_csv(popmap, sep="\t", header=None, names=["Individual","Population"])
    df_spc = pd.read_csv(speciesmap, sep="\t", header=None, names=["Individual","Species"])
    return df_map, df_pop, df_spc

def assemble(df_nh, df_map, df_pop, df_spc):
    df = df_nh.drop(columns="Individual").merge(df_map, on="Index")
    df = df.merge(df_pop, on="Individual", how="left")
    df = df.merge(df_spc, on="Individual", how="left")
    return df

def make_plot(df, out_html, template=None, palette="Viridis"):
    cats = ["P0","Bx0","F1","F2","Bx1","P1"]

    # category assignment and priority
    df["AssignedCategory"] = df[cats].idxmax(axis=1)

    # species ordering by mean P0
    species_order = df.groupby("Species")["P0"] \
                      .mean() \
                      .sort_values(ascending=False) \
                      .index.tolist()

    # individual ordering per species and category
    individual_order = []
    for sp in species_order:
        sub_sp = df[df["Species"] == sp]
        for cat in cats:
            sub_cat = sub_sp[sub_sp["AssignedCategory"] == cat]
            ascending = (cat == "P1")
            sub_sorted = sub_cat.sort_values(by=cat, ascending=ascending)
            individual_order += sub_sorted["Individual"].tolist()

    # species tick positions
    sp_map = df[["Individual","Species"]].drop_duplicates().set_index("Individual")
    species = sp_map.loc[individual_order]["Species"]
    counts  = species.value_counts().loc[species.unique()]
    ticks   = counts.cumsum() - counts/2

    # melt for plotting
    id_vars = ["Individual","Species","Population","AssignedCategory"]
    df_long = df.melt(id_vars=id_vars, value_vars=cats,
                      var_name="Category", value_name="Probability")
    df_long["Individual"] = pd.Categorical(df_long["Individual"],
                                           categories=individual_order,
                                           ordered=True)

    # color sequence
    n = len(cats)
    if hasattr(px.colors.qualitative, palette):
        color_seq = getattr(px.colors.qualitative, palette)
    elif hasattr(px.colors.sequential, palette):
        scale     = getattr(px.colors.sequential, palette)
        color_seq = px.colors.sample_colorscale(scale, [i/(n-1) for i in range(n)])
    elif hasattr(px.colors.diverging, palette):
        scale     = getattr(px.colors.diverging, palette)
        color_seq = px.colors.sample_colorscale(scale, [i/(n-1) for i in range(n)])
    else:
        raise ValueError(f"palette '{palette}' not found")

    # hover info
    hover = ["Individual","Species","Population","AssignedCategory","Category","Probability"]

    # bar plot
    fig = px.bar(df_long,
                 x="Individual", y="Probability", color="Category",
                 category_orders={"Individual": individual_order, "Category": cats},
                 color_discrete_sequence=color_seq,
                 hover_data=hover)
    fig.update_traces(marker_line_width=0)

    # layout
    fig.update_layout(
        barmode="stack",
        yaxis=dict(title="Posterior Probability", range=[0,1], showgrid=False),
        xaxis=dict(tickmode="array",
                   tickvals=ticks.values,
                   ticktext=ticks.index,
                   title="Species",
                   showgrid=False),
        margin=dict(t=60,b=60),
        title=dict(text="NewHybrids Classifications", x=0.5),
        legend_title="Category",
        template="simple_white"
    )

    html = fig.to_html(full_html=False, include_plotlyjs="cdn")
    if template:
        header = Path(template).read_text()
        Path(out_html).write_text(header + html)
    else:
        Path(out_html).write_text(html)

    print(f"✅ Plot saved to: {out_html}")

if __name__=="__main__":
    p = argparse.ArgumentParser(description="Plot NewHybrids classifications")
    p.add_argument("--result",     required=True, help="NH posterior file")
    p.add_argument("--result_map", required=True, help="Index→sample TSV")
    p.add_argument("--popmap",     required=True, help="Sample→population TSV")
    p.add_argument("--speciesmap", required=True, help="Sample→species TSV")
    p.add_argument("--template",   help="HTML header to prepend")
    p.add_argument("--palette",    default="Viridis", help="plotly.colors palette")
    p.add_argument("--out",        required=True, help="Output HTML file")

    args = p.parse_args()
    df_nh      = read_nh_results(args.result)
    df_map, df_pop, df_spc = load_maps(args.result_map, args.popmap, args.speciesmap)
    df         = assemble(df_nh, df_map, df_pop, df_spc)
    make_plot(df, args.out, template=args.template, palette=args.palette)
