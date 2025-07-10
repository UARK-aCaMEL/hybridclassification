#!/usr/bin/env python3
import pandas as pd
import plotly.express as px
import argparse
from pathlib import Path

def load_trace(trace_file):
    """
    Read a PI_TRACE file like:

      PI_TRACE:Rep  Pure_0  Pure_1  F1   F2   0_BX  1_BX
      PI_TRACE:0    0.45    0.43    …
      PI_TRACE:1    0.44    0.40    …
      …

    Returns a DataFrame indexed by numeric iteration, with one column per category.
    """
    # whitespace‐separated, header=0
    df = pd.read_csv(trace_file, sep="\t", header=0)
    rep_col = df.columns[0]
    # extract the iteration number after the colon
    print(df)
    df['Iteration'] = df[rep_col].str.split(':', 1).str[1].astype(int)
    df = df.drop(columns=[rep_col]).set_index('Iteration')
    return df

def make_trace_plot(df, output_html, burnin, template_file=None):
    """
    Melt and plot all categories as lines, add a vertical red line at `burnin`,
    and write HTML (optionally with a prepended template header).
    """
    # long form for px.line
    df_long = df.reset_index().melt(
        id_vars='Iteration',
        var_name='Category',
        value_name='Probability'
    )

    fig = px.line(
        df_long,
        x='Iteration', y='Probability', color='Category',
        title='NewHybrids PI Trace',
    )
    fig.update_layout(
        xaxis_title='Iteration',
        yaxis_title='π Probability',
        template='simple_white',
        legend_title='Category',
        margin=dict(t=60, b=60)
    )
    # add burn‐in marker
    fig.add_vline(x=burnin, line_color='red', line_dash='dash')

    # generate HTML snippet
    html_body = fig.to_html(full_html=False, include_plotlyjs='cdn')

    # prepend template if requested
    if template_file:
        header = Path(template_file).read_text()
        Path(output_html).write_text(header + html_body)
    else:
        Path(output_html).write_text(html_body)

    print(f"✅ Trace plot saved to: {output_html}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Plot NewHybrids PI MCMC trace for each category."
    )
    parser.add_argument(
        "--trace", required=True,
        help="PI_TRACE file (whitespace‐separated; first column like PI_TRACE:<iter>)."
    )
    parser.add_argument(
        "--out", required=True,
        help="Output HTML file path."
    )
    parser.add_argument(
        "--template",
        help="Optional HTML header/template file to prepend."
    )
    parser.add_argument(
        "--burnin", type=int, required=True,
        help="Iteration at which burn‐in ends (vertical red line)."
    )
    args = parser.parse_args()

    df = load_trace(args.trace)
    make_trace_plot(df, args.out, burnin=args.burnin, template_file=args.template)
