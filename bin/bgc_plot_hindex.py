#!/usr/bin/env python3

import re
import argparse
import pandas as pd
import plotly.graph_objs as go
import plotly.io as pio


def parse_template(template_file):
    """Extract metadata from the MultiQC-style HTML comment block."""
    with open(template_file, "r") as f:
        content = f.read()

    match = re.search(r"<!--(.*?)-->", content, re.DOTALL)
    if not match:
        raise ValueError("No metadata block found in template.")

    raw_block = match.group(1)
    meta = {}
    for line in raw_block.strip().splitlines():
        if ":" in line:
            key, value = line.strip().split(":", 1)
            meta[key.strip()] = value.strip().strip("\"'")
    return meta


def build_comment(meta_dict):
    lines = ["<!--"]
    for key, value in meta_dict.items():
        lines.append(f'{key}: "{value}"')
    lines.append("-->")
    return "\n".join(lines)


def read_order_file(order_file):
    """Read sample names from the input order file."""
    with open(order_file, "r") as f:
        samples = [line.strip() for line in f if line.strip()]

    if not samples:
        raise ValueError(f"No sample IDs found in order file: {order_file}")

    df = pd.DataFrame(
        {
            "input_order": range(1, len(samples) + 1),
            "sample": samples,
        }
    )
    return df


def read_hindex_file(hindex_file):
    """
    Read hybrid index table.

    Expected format:
        <row_id>  50%  2.5%  5%  95%  97.5%

    where the first column is the 1-based row index corresponding to
    the input sample order.
    """
    df = pd.read_csv(hindex_file, sep="\t", index_col=0)

    expected_cols = {"50%", "2.5%", "97.5%"}
    missing = expected_cols - set(df.columns)
    if missing:
        raise ValueError(
            f"Missing required columns in hindex file {hindex_file}: {sorted(missing)}"
        )

    df = df.reset_index()
    df = df.rename(columns={df.columns[0]: "input_order"})
    df["input_order"] = pd.to_numeric(df["input_order"], errors="raise")

    return df


def prepare_plot_df(order_file, hindex_file):
    order_df = read_order_file(order_file)
    hi_df = read_hindex_file(hindex_file)

    df = order_df.merge(hi_df, on="input_order", how="inner")

    if df.shape[0] != order_df.shape[0]:
        missing_n = order_df.shape[0] - df.shape[0]
        raise ValueError(
            f"Order file has {order_df.shape[0]} samples but only {df.shape[0]} "
            f"matched rows were found in the hindex file ({missing_n} missing)."
        )

    df["median_hi"] = pd.to_numeric(df["50%"], errors="raise")
    df["ci_lower"] = pd.to_numeric(df["2.5%"], errors="raise")
    df["ci_upper"] = pd.to_numeric(df["97.5%"], errors="raise")

    df["err_minus"] = df["median_hi"] - df["ci_lower"]
    df["err_plus"] = df["ci_upper"] - df["median_hi"]

    # Sort for plotting by ascending hybrid index
    df = df.sort_values(["median_hi", "sample"], ascending=[True, True]).reset_index(drop=True)

    return df


def generate_plot(order_file, hindex_file, output_file, header_comment):
    df = prepare_plot_df(order_file, hindex_file)

    # Reverse y so smallest hybrid index appears at the bottom and largest at the top
    # If you want smallest at the top, remove the [::-1] behavior and autorange reversal.
    plot_df = df.copy()

    fig = go.Figure()

    fig.add_trace(
        go.Scatter(
            x=plot_df["median_hi"],
            y=plot_df["sample"],
            mode="markers",
            name="Hybrid index",
            error_x=dict(
                type="data",
                symmetric=False,
                array=plot_df["err_plus"],
                arrayminus=plot_df["err_minus"],
                visible=True,
                thickness=1.2,
                width=0,
            ),
            marker=dict(size=7),
            customdata=plot_df[["input_order", "median_hi", "ci_lower", "ci_upper"]].to_numpy(),
            hovertemplate=(
                "<b>%{y}</b><br>"
                "Input order: %{customdata[0]}<br>"
                "Hybrid index (50%%): %{customdata[1]:.4f}<br>"
                "95%% CI: [%{customdata[2]:.4f}, %{customdata[3]:.4f}]<extra></extra>"
            ),
        )
    )

    fig.update_layout(
        title="Hybrid Index Estimates with 95% Credible Intervals",
        xaxis_title="Hybrid Index",
        yaxis_title="Sample",
        template="plotly_white",
        height=max(500, len(plot_df) * 18),
        hovermode="closest",
        showlegend=False,
    )

    fig.update_xaxes(range=[-0.02, 1.02])

    # For many samples, labels can get crowded. This still keeps them available in hover.
    if len(plot_df) > 75:
        fig.update_yaxes(showticklabels=False)
    else:
        fig.update_yaxes(automargin=True)

    html = pio.to_html(fig, full_html=False, include_plotlyjs="cdn")
    with open(output_file, "w") as f:
        f.write(header_comment + "\n" + html)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate hybrid index dot-and-whisker plot as HTML with embedded MultiQC metadata."
    )
    parser.add_argument(
        "--order",
        required=True,
        help="Input order file with one sample name per line",
    )
    parser.add_argument(
        "--hindex",
        required=True,
        help="Hybrid index results TSV file",
    )
    parser.add_argument(
        "--out",
        default="hindex_plot.html",
        help="Output HTML filename",
    )
    parser.add_argument(
        "--template",
        required=True,
        help="Path to HTML file containing MultiQC metadata comment block",
    )
    parser.add_argument("--id", help="Override for the 'id' field in the metadata")
    parser.add_argument("--title", help="Override for the 'title' field in the metadata")
    parser.add_argument("--section_name", help="Override for 'section_name'")
    parser.add_argument("--description", help="Override for 'description'")

    args = parser.parse_args()

    metadata = parse_template(args.template)
    if args.id:
        metadata["id"] = args.id
    if args.title:
        metadata["title"] = args.title
    if args.section_name:
        metadata["section_name"] = args.section_name
    if args.description:
        metadata["description"] = args.description

    header_comment = build_comment(metadata)
    generate_plot(args.order, args.hindex, args.out, header_comment)
