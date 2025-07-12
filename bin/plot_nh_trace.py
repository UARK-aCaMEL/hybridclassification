#!/usr/bin/env python3
import pandas as pd
import plotly.express as px
import argparse
from pathlib import Path

def load_trace(trace_file):
    df = pd.read_csv(trace_file, sep="\t", header=0)
    rep_col = df.columns[0]
    df['Iteration'] = df[rep_col].str.replace("PI_TRACE:", "", regex=False).astype(int)
    return df.drop(columns=[rep_col]).set_index('Iteration')

def make_trace_plot(df, output_html, burnin, template_file=None):
    df_long = df.reset_index().melt(
        id_vars='Iteration', var_name='Category', value_name='Probability'
    )
    fig = px.line(
        df_long,
        x='Iteration', y='Probability', color='Category',
        title='NewHybrids PI Trace',
        template='simple_white'
    )
    fig.update_layout(
        xaxis_title='Iteration',
        yaxis_title='π Probability',
        legend_title='Category',
        margin=dict(t=60, b=60)
    )
    # explicit vertical line at burnin
    fig.add_shape(
        type="line",
        x0=burnin, x1=burnin,
        y0=0, y1=1,
        xref="x", yref="paper",
        line=dict(color="red", dash="dash", width=2),
        layer="above"
    )
    # label the burnin
    fig.add_annotation(
        x=burnin, y=1.0,
        xref="x", yref="paper",
        text=f"burnin={burnin}",
        showarrow=False,
        yanchor="bottom",
        font=dict(color="red", size=12)
    )

    html_body = fig.to_html(full_html=False, include_plotlyjs='cdn')
    if template_file:
        header = Path(template_file).read_text().rstrip() + "\n"
        Path(output_html).write_text(header + html_body)
    else:
        Path(output_html).write_text(html_body)
    print(f"✅ Trace plot saved to: {output_html}")

if __name__ == "__main__":
    p = argparse.ArgumentParser(
        description="Plot NewHybrids PI MCMC trace for each category."
    )
    p.add_argument("--trace",    required=True, help="PI_TRACE file")
    p.add_argument("--out",      required=True, help="Output HTML path")
    p.add_argument("--template", help="Optional HTML header file")
    p.add_argument("--burnin",   type=int, required=True, help="Burn-in iteration")
    args = p.parse_args()

    df = load_trace(args.trace)
    make_trace_plot(df, args.out, burnin=args.burnin, template_file=args.template)
