#!/usr/bin/env python3
import argparse
import re
from pathlib import Path

import pandas as pd
import plotly.express as px


def load_trace(trace_file):
    cleaned_rows = []

    with open(trace_file, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue

            # Keep only actual trace rows
            if not line.startswith("PI_TRACE:"):
                continue

            # Extract iteration number from start of line
            m = re.match(r"^PI_TRACE:(\d+)", line)
            if not m:
                continue
            iteration = int(m.group(1))

            # Remove everything before/including PI_TRACE:<n>
            rest = re.sub(r"^PI_TRACE:\d+\s*", "", line)

            # Remove any embedded WARNING text and everything after it
            rest = re.sub(r"WARNING:.*$", "", rest).strip()

            # Split remaining fields on whitespace
            vals = re.split(r"\s+", rest)

            # Require at least one numeric value
            if not vals:
                continue

            cleaned_rows.append([iteration] + vals)

    if not cleaned_rows:
        raise ValueError(f"No valid PI_TRACE rows found in {trace_file}")

    # Determine expected number of value columns from the most common row length
    lengths = pd.Series([len(r) for r in cleaned_rows])
    expected_len = lengths.mode().iloc[0]

    # Keep only rows matching the common width
    cleaned_rows = [r for r in cleaned_rows if len(r) == expected_len]

    if not cleaned_rows:
        raise ValueError(f"No valid PI_TRACE rows with consistent width found in {trace_file}")

    ncat = expected_len - 1
    cols = ["Iteration"] + [f"Category_{i+1}" for i in range(ncat)]
    df = pd.DataFrame(cleaned_rows, columns=cols)

    # Convert numeric columns safely
    for c in cols[1:]:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    # Drop malformed rows
    df = df.dropna()

    # Sort and index
    df = df.sort_values("Iteration").drop_duplicates(subset=["Iteration"])
    return df.set_index("Iteration")


def make_trace_plot(df, output_html, burnin, template_file=None):
    df_long = df.reset_index().melt(
        id_vars="Iteration", var_name="Category", value_name="Probability"
    )

    fig = px.line(
        df_long,
        x="Iteration",
        y="Probability",
        color="Category",
        title="",
        template="simple_white",
    )

    fig.update_layout(
        xaxis_title="Iteration",
        yaxis_title="π",
        legend_title="Category",
        margin=dict(t=60, b=60),
    )

    fig.add_shape(
        type="line",
        x0=burnin,
        x1=burnin,
        y0=0,
        y1=1,
        xref="x",
        yref="paper",
        line=dict(color="red", dash="dash", width=2),
        layer="above",
    )

    fig.add_annotation(
        x=burnin,
        y=1.0,
        xref="x",
        yref="paper",
        text=f"burnin={burnin}",
        showarrow=False,
        yanchor="bottom",
        font=dict(color="red", size=12),
    )

    html_body = fig.to_html(full_html=False, include_plotlyjs="cdn")
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
    p.add_argument("--trace", required=True, help="PI_TRACE file")
    p.add_argument("--out", required=True, help="Output HTML path")
    p.add_argument("--template", help="Optional HTML header file")
    p.add_argument("--burnin", type=int, required=True, help="Burn-in iteration")
    args = p.parse_args()

    df = load_trace(args.trace)
    make_trace_plot(df, args.out, burnin=args.burnin, template_file=args.template)
