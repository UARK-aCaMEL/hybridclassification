#!/usr/bin/env python3
"""
nh_summary_table.py
------------------
Generate two MultiQC-compatible summary tables of NewHybrids assignments:
 1. Default summary of P0, P1, F1, F2, Bx0, Bx1, Unassigned, Total_Hybrids, and N
 2. Masked summary (optional) where samples listed in --mask are forced to Unassigned
    before computing proportions.

Usage:
    nh_summary_table.py \
      --result       nh_posteriors.txt \
      --result_map   index_map.tsv \
      --popmap       sample_pop.tsv \
      --speciesmap   sample_sp.tsv \
      --template     nh_summary_header.html \
      --out          nh_summary.json \
      [--threshold 0.5] \
      [--list       hybrid_list.txt] \
      [--mask       masked_ids.txt] \
      [--masked_template nh_masked_header.html] \
      [--out_mask   nh_summary_masked.json] \
      [--list_masked masked_hybrid_list.txt]
"""
import pandas as pd
import json
import argparse
import numpy as np
import re
from pathlib import Path


def read_nh_results(path):
    cats = ["P0","P1","F1","F2","Bx0","Bx1"]
    names = ["Index","Individual"] + cats
    return pd.read_csv(path, sep=r'\s+', skiprows=1, header=None, names=names)


def load_maps(nh_map, popmap, speciesmap):
    df_map = (pd.read_csv(nh_map, sep="\t", header=0)
                .rename(columns={"Sample":"Individual"}))
    df_pop = pd.read_csv(popmap, sep="\t", header=None, names=["Individual","Population"])
    df_spc = pd.read_csv(speciesmap, sep="\t", header=None, names=["Individual","Species"])
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
        row['Group']: {k: v for k, v in row.items() if k not in ('Group','N')}
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


def compute_summary(df, cats, hybrids):
    # species-level
    all_cats = cats + ['Unassigned']
    sp_counts = df.groupby('Species')['AssignedCategory'] \
                   .value_counts().unstack(fill_value=0)
    sp_n      = sp_counts.sum(axis=1)
    sp_prop   = sp_counts.div(sp_n, axis=0)
    for c in all_cats:
        if c not in sp_prop.columns:
            sp_prop[c] = 0.0
    sp_prop['Total_Hybrids'] = sp_prop[hybrids].sum(axis=1)
    sp_prop = sp_prop.reset_index().rename(columns={'Species':'Group'})
    sp_prop['N'] = sp_n.values.astype(int)
    sp_prop = sp_prop[['Group','N'] + cats + ['Unassigned','Total_Hybrids']]

    # population-level
    pop_counts = df.groupby(['Species','Population'])['AssignedCategory'] \
                    .value_counts().unstack(fill_value=0)
    pop_n      = pop_counts.sum(axis=1)
    pop_prop   = pop_counts.div(pop_n, axis=0)
    for c in all_cats:
        if c not in pop_prop.columns:
            pop_prop[c] = 0.0
    pop_prop['Total_Hybrids'] = pop_prop[hybrids].sum(axis=1)
    pop_prop = pop_prop.reset_index()
    pop_prop['Group'] = pop_prop['Species'] + '|' + pop_prop['Population']
    pop_n_df = pop_n.reset_index(name='N')
    pop_prop = pop_prop.merge(pop_n_df, on=['Species','Population'])
    pop_prop['N'] = pop_prop['N'].astype(int)
    pop_prop = pop_prop[['Group','N'] + cats + ['Unassigned','Total_Hybrids']]

    summary = pd.concat([sp_prop, pop_prop], ignore_index=True)
    summary[cats + ['Unassigned','Total_Hybrids']] = \
        summary[cats + ['Unassigned','Total_Hybrids']].round(2)
    return summary


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--result',           required=True, help='NH posterior file')
    p.add_argument('--result_map',       required=True, help='Index→sample TSV')
    p.add_argument('--popmap',           required=True, help='Sample→population TSV')
    p.add_argument('--speciesmap',       required=True, help='Sample→species TSV')
    p.add_argument('--template',         help='MultiQC HTML header for default table')
    p.add_argument('--out',              required=True, help='Output path (tsv or JSON) for default table')
    p.add_argument('--threshold',        type=float, default=0.0,
                   help='Minimum posterior probability; else Unassigned')
    p.add_argument('--list',             help='Optional output path for list of hybrids')
    p.add_argument('--mask',             help='File listing individuals to mask (one per line)')
    p.add_argument('--masked_template',  help='MultiQC HTML header for masked table')
    p.add_argument('--out_mask',         help='Output path (tsv or JSON) for masked table')
    p.add_argument('--list_masked',      help='Optional output path for list of hybrids after masking')
    args = p.parse_args()

    # read and merge inputs
    df_nh      = read_nh_results(args.result)
    df_map, df_pop, df_spc = load_maps(
        args.result_map, args.popmap, args.speciesmap
    )
    df = (
        df_nh.drop(columns='Individual')
             .merge(df_map, on='Index')
             .merge(df_pop, on='Individual', how='left')
             .merge(df_spc, on='Individual', how='left')
    )

    # categories
    cats    = ['P0','P1','F1','F2','Bx0','Bx1']
    hybrids = ['F1','F2','Bx0','Bx1']

    # assign categories
    df['MaxProb'] = df[cats].max(axis=1)
    df['AssignedCategory'] = df[cats].idxmax(axis=1)
    df.loc[df['MaxProb'] <= args.threshold, 'AssignedCategory'] = 'Unassigned'

    # write default hybrid list
    if args.list:
        is_hybrid = df['AssignedCategory'].isin(hybrids)
        hybrid_df = df[is_hybrid].copy()
        hybrid_df['Prob'] = hybrid_df.apply(
            lambda r: r[r['AssignedCategory']], axis=1
        )
        hybrid_df['Individual'].to_csv(args.list, index=False, header=False)
        print(f"✅ Hybrid list → {args.list}")

    # default summary
    summary = compute_summary(df, cats, hybrids)

    # write default table
    if args.template:
        meta = parse_html_header(args.template)
        write_mqc_json(summary, meta, args.out)
    else:
        summary.to_csv(args.out, sep='\t', index=False, float_format='%.2f')
    print(f"✅ Default summary table → {args.out}")

    # masked summary (optional)
    if args.mask:
        mask_ids = set(Path(args.mask).read_text().split())
        if mask_ids & set(df['Individual']):
            df_masked = df.copy()
            df_masked.loc[
                df_masked['Individual'].isin(mask_ids),
                'AssignedCategory'
            ] = 'Unassigned'
            summary_m = compute_summary(df_masked, cats, hybrids)

            # determine masked summary output path
            out_m = args.out_mask or \
                str(Path(args.out).with_name(
                    Path(args.out).stem + '.masked' + Path(args.out).suffix
                ))

            if args.masked_template:
                meta_m = parse_html_header(args.masked_template)
                write_mqc_json(summary_m, meta_m, out_m)
            else:
                summary_m.to_csv(out_m, sep='\t', index=False, float_format='%.2f')
            print(f"✅ Masked summary table → {out_m}")

            # write masked hybrid list
            if args.list_masked:
                is_hybrid_m = df_masked['AssignedCategory'].isin(hybrids)
                hybrid_m = df_masked[is_hybrid_m].copy()
                hybrid_m['Prob'] = hybrid_m.apply(
                    lambda r: r[r['AssignedCategory']], axis=1
                )
                hybrid_m['Individual'].to_csv(
                    args.list_masked, index=False, header=False
                )
                print(f"✅ Masked hybrid list → {args.list_masked}")

if __name__ == '__main__':
    main()
