#!/usr/bin/env python3
import pandas as pd
import json
import argparse
import re
from pathlib import Path

def read_nh_results(path):
    cats = ["P0","P1","F1","F2","Bx0","Bx1"]
    names = ["Index","Individual"] + cats
    return pd.read_csv(path, sep=r'\s+', skiprows=1, header=None, names=names)

def load_maps(nh_map, popmap, speciesmap):
    df_map = pd.read_csv(nh_map, sep="\t", header=0)\
               .rename(columns={"Sample":"Individual"})
    df_pop = pd.read_csv(popmap, sep="\t", header=None,
                         names=["Individual","Population"])
    df_spc = pd.read_csv(speciesmap, sep="\t", header=None,
                         names=["Individual","Species"])
    return df_map, df_pop, df_spc

def parse_html_header(path):
    meta = {}
    for line in open(path):
        txt = line.strip().lstrip("<!--").rstrip("-->").strip()
        m = re.match(r'([A-Za-z0-9_]+):\s*"(.*)"', txt)
        if m:
            meta[m.group(1)] = m.group(2)
    return meta

def write_mqc_json(df, metadata, output):
    data = {
        row['Group']: {k: v for k, v in row.items() if k != 'Group'}
        for _, row in df.iterrows()
    }
    pconfig = {
        'id':       metadata.get('id'),
        'ylab':     'Proportion',
        'xlab':     'Group',
        'xDecimals': False,
        'min':      0.0,
        'max':      1.0,
        'scale':    'YlGnBu'
    }
    out = {'data': data, 'pconfig': pconfig}
    out.update({k:v for k,v in metadata.items() if k!='id'})
    with open(output, 'w') as f:
        json.dump(out, f, indent=2)

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--result',      required=True, help='NH posterior file')
    p.add_argument('--result_map',  required=True, help='Index→sample TSV')
    p.add_argument('--popmap',      required=True, help='Sample→population TSV')
    p.add_argument('--speciesmap',  required=True, help='Sample→species TSV')
    p.add_argument('--template',    help='MultiQC HTML header')
    p.add_argument('--out',         required=True, help='Output path (tsv or JSON)')
    args = p.parse_args()

    # read and merge inputs
    df_nh      = read_nh_results(args.result)
    df_map, df_pop, df_spc = load_maps(args.result_map, args.popmap, args.speciesmap)
    df = (
        df_nh.drop(columns="Individual")
             .merge(df_map, on="Index")
             .merge(df_pop, on="Individual", how="left")
             .merge(df_spc, on="Individual", how="left")
    )

    cats    = ["P0","P1","F1","F2","Bx0","Bx1"]
    hybrids = ["F1","F2","Bx0","Bx1"]

    # assign each to highest-probability category
    df["AssignedCategory"] = df[cats].idxmax(axis=1)

    # species-level counts & proportions
    sp_counts = df.groupby("Species")["AssignedCategory"] \
                  .value_counts().unstack(fill_value=0)
    sp_n      = sp_counts.sum(axis=1)
    sp_prop   = sp_counts.div(sp_n, axis=0)
    # ensure all categories present
    for c in cats:
        if c not in sp_prop.columns:
            sp_prop[c] = 0.0
    sp_prop["Total_Hybrids"] = sp_prop[hybrids].sum(axis=1)
    sp_prop = sp_prop.reset_index().rename(columns={"Species":"Group"})
    sp_prop["N"] = sp_n.values.astype(int)
    sp_prop = sp_prop[["Group","N"] + cats + ["Total_Hybrids"]]

    # population-level counts & proportions
    pop_counts = df.groupby(["Species","Population"])["AssignedCategory"] \
                   .value_counts().unstack(fill_value=0)
    pop_n      = pop_counts.sum(axis=1)
    pop_prop   = pop_counts.div(pop_n, axis=0)
    for c in cats:
        if c not in pop_prop.columns:
            pop_prop[c] = 0.0
    pop_prop["Total_Hybrids"] = pop_prop[hybrids].sum(axis=1)
    pop_prop = pop_prop.reset_index()
    pop_prop["Group"] = pop_prop["Species"] + "|" + pop_prop["Population"]
    pop_n_df = pop_n.reset_index(name="N")
    pop_prop = pop_prop.merge(pop_n_df, on=["Species","Population"])
    pop_prop["N"] = pop_prop["N"].astype(int)
    pop_prop = pop_prop[["Group","N"] + cats + ["Total_Hybrids"]]

    # combine & format proportions to two decimals
    summary = pd.concat([sp_prop, pop_prop], ignore_index=True)
    summary[cats + ["Total_Hybrids"]] = summary[cats + ["Total_Hybrids"]].round(2)

    if args.template:
        meta = parse_html_header(args.template)
        write_mqc_json(summary, meta, args.out)
    else:
        summary.to_csv(args.out, sep="\t", index=False, float_format="%.2f")

if __name__ == '__main__':
    main()
