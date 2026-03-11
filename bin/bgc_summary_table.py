#!/usr/bin/env python3
"""
bgc_summary_table.py
--------------------
Parse Stan summary text files from BGC-HM and combine them into a single
summary table suitable for MultiQC custom content.

Parsing rule used here:
  - the Stan summary table starts after the first blank line in the file
  - the table ends at the second blank line after that
  - within that table:
      * lines beginning with leading whitespace are header lines
      * lines not beginning with leading whitespace are parameter rows
  - header lines may appear multiple times (e.g. mean...n_eff, then later Rhat)
  - parameter rows are merged across header blocks by parameter name

This handles cases where:
  - mean..97.5% are one block
  - n_eff is inline with the first block or on its own later block
  - Rhat is on its own later block
  - blocks are interleaved within the same summary table region

Usage:
    bgc_summary_table.py \
      --gencline CAMANO_CHRERY__gencline_hmc__summary.txt \
      --hindex   CAMANO_CHRERY__hi_hmc__summary.txt \
      --template multiqc_bgc_summary.html \
      --out      bgc_summary_mqc.json
"""

import argparse
import json
import re
from pathlib import Path

import pandas as pd


FINAL_ORDER = [
    "RowName",
    "Model",
    "Parameter",
    "n_eff",
    "Rhat",
    "mean",
    "se_mean",
    "sd",
    "2.5%",
    "25%",
    "50%",
    "75%",
    "97.5%",
]


def parse_html_header(path):
    meta = {}
    with open(path) as fh:
        for line in fh:
            txt = line.strip().lstrip("<!--").rstrip("-->").strip()
            m = re.match(r'([A-Za-z0-9_]+):\s*"(.*)"', txt)
            if m:
                meta[m.group(1)] = m.group(2)
    return meta


def is_float(x):
    try:
        float(x)
        return True
    except Exception:
        return False


def normalize_numeric_columns(df):
    for col in df.columns:
        if col not in {"RowName", "Model", "Parameter"}:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def extract_table_lines(lines):
    """
    Extract only the Stan summary table region:
      - starts after first blank line
      - ends at second blank line
    """
    blank_count = 0
    in_table = False
    table_lines = []

    for line in lines:
        if line.strip() == "":
            blank_count += 1
            if not in_table and blank_count == 1:
                in_table = True
                continue
            if in_table:
                break
        else:
            if in_table:
                table_lines.append(line.rstrip("\n"))

    return table_lines


def parse_stan_summary(path, model_label):
    lines = Path(path).read_text().splitlines()
    table_lines = extract_table_lines(lines)

    if not table_lines:
        raise ValueError(f"Could not extract Stan summary table from: {path}")

    rows = {}
    current_cols = None

    for line in table_lines:
        if not line.strip():
            continue

        # Header lines start with leading whitespace
        if line[:1].isspace():
            cols = line.split()
            if not cols:
                continue
            current_cols = cols
            continue

        # Data lines do not start with leading whitespace
        if current_cols is None:
            continue

        parts = line.split()
        if len(parts) < 1 + len(current_cols):
            continue

        param = parts[0]
        vals = parts[1:1 + len(current_cols)]

        if not all(is_float(v) for v in vals):
            continue

        if param not in rows:
            rows[param] = {"Parameter": param}

        for col, val in zip(current_cols, vals):
            rows[param][col] = float(val)

    df = pd.DataFrame(rows.values())

    if df.empty:
        raise ValueError(f"No Stan summary rows parsed from: {path}")

    for col in ["n_eff", "Rhat", "mean", "se_mean", "sd", "2.5%", "25%", "50%", "75%", "97.5%"]:
        if col not in df.columns:
            df[col] = pd.NA

    df["Model"] = model_label
    df["RowName"] = df["Model"] + "|" + df["Parameter"]

    df = df[[c for c in FINAL_ORDER if c in df.columns]]
    df = normalize_numeric_columns(df)

    return df


def build_headers(df):
    headers = {}
    placement = 1

    for col in df.columns:
        if col == "RowName":
            continue

        header = {"title": col, "placement": placement}

        if col == "n_eff":
            header["format"] = "{:,.0f}"
        elif col in {"Rhat", "mean", "se_mean", "sd", "2.5%", "25%", "50%", "75%", "97.5%"}:
            header["format"] = "{:.2f}"

        headers[col] = header
        placement += 1

    return headers


def write_mqc_json(df, metadata, output):
    data = {
        row["RowName"]: {
            k: (None if pd.isna(v) else v)
            for k, v in row.items()
            if k != "RowName"
        }
        for _, row in df.iterrows()
    }

    pconfig = {
        "id": metadata.get("id", "bgc_summary_table"),
        "title": metadata.get("section_name", metadata.get("title", "BGC-HM Stan Summary Table")),
    }

    out = {
        "data": data,
        "headers": build_headers(df),
        "pconfig": pconfig,
    }
    out.update({k: v for k, v in metadata.items() if k != "id"})

    with open(output, "w") as fh:
        json.dump(out, fh, indent=2)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--gencline", required=True, help="Stan summary text file for gencline")
    p.add_argument("--hindex", required=True, help="Stan summary text file for h-index")
    p.add_argument("--label_gencline", default="gencline", help="Label for gencline rows")
    p.add_argument("--label_hindex", default="hindex", help="Label for hindex rows")
    p.add_argument("--template", help="MultiQC HTML header template")
    p.add_argument("--out", required=True, help="Output path (.json or .tsv)")
    args = p.parse_args()

    df_gencline = parse_stan_summary(args.gencline, args.label_gencline)
    df_hindex = parse_stan_summary(args.hindex, args.label_hindex)

    summary = pd.concat([df_gencline, df_hindex], ignore_index=True)
    summary = summary[[c for c in FINAL_ORDER if c in summary.columns]]

    if args.template:
        meta = parse_html_header(args.template)
        write_mqc_json(summary, meta, args.out)
    else:
        summary.to_csv(args.out, sep="\t", index=False, float_format="%.2f")

    print(f"✅ Combined Stan summary table → {args.out}")


if __name__ == "__main__":
    main()
