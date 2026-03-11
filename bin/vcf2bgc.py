#!/usr/bin/env python3

from __future__ import annotations
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

BASES = ("A", "C", "G", "T")
BASE_TO_I = {b: i for i, b in enumerate(BASES)}


def load_popmap(path: Path) -> dict[str, str]:
    df = pd.read_csv(path, sep="\t", header=None, names=["sample", "pop"], dtype=str)
    return dict(zip(df["sample"], df["pop"]))


def recode(rec):
    v = np.full(len(rec.samples), np.nan, dtype=np.float32)
    for i, call in enumerate(rec.samples.values()):
        gt = call["GT"]
        if gt and None not in gt and min(gt) >= 0:
            v[i] = gt[0] + gt[1]
    return v


def vcf_matrix(vcf: Path):
    fh = pysam.VariantFile(str(vcf))
    samp, loci, g = list(fh.header.samples), [], []
    for r in fh.fetch():
        loci.append(r.id if r.id and r.id != "." else f"{r.contig}_{r.pos}")
        g.append(recode(r))
    if not g:
        sys.exit(f"{vcf} is empty")
    return samp, loci, np.vstack(g).T


def pseudo_glik_for_gt(gt, ref_base: str, alt_base: str, e: float) -> np.ndarray:
    out = np.full(4, 0.25, dtype=np.float32)

    # missing / invalid diploid
    if not gt or None in gt or len(gt) != 2 or min(gt) < 0:
        return out

    # only handle biallelic 0/1; anything else -> missing/uniform
    if gt[0] > 1 or gt[1] > 1:
        return out

    a0 = ref_base if gt[0] == 0 else alt_base
    a1 = ref_base if gt[1] == 0 else alt_base

    if a0 not in BASE_TO_I or a1 not in BASE_TO_I:
        return out

    if a0 == a1:
        out[:] = e / 3.0
        out[BASE_TO_I[a0]] = 1.0 - e
    else:
        out[:] = 0.0
        out[BASE_TO_I[a0]] = (1.0 - e) / 2.0
        out[BASE_TO_I[a1]] = (1.0 - e) / 2.0
        for b in BASES:
            if b != a0 and b != a1:
                out[BASE_TO_I[b]] = e / 2.0

    return out


def r_matrix_text(mat: np.ndarray) -> str:
    """
    Return an R matrix expression that prints like your example:
      structure(c(...), .Dim=c(nrow,ncol))
    IMPORTANT: R is column-major, so we flatten in Fortran order.
    """
    nrow, ncol = mat.shape
    flat = mat.reshape(-1, order="F")
    chunks = []
    for i in range(0, flat.size, 2000):
        chunk = ", ".join(f"{x:.6f}" for x in flat[i:i+2000])
        chunks.append(chunk)
    vec = ",\n  ".join(chunks)
    return f"structure(c(\n  {vec}\n), .Dim=c({nrow},{ncol}))"


def write_glik_r(path: Path, mats_by_base: dict[str, np.ndarray]):
    """
    Write an R-readable text file for dget(), containing:
      list( <matA>, <matC>, <matG>, <matT> )
    """
    mat_exprs = [r_matrix_text(mats_by_base[b]) for b in BASES]
    body = ",\n".join(mat_exprs)
    txt = f"list(\n{body}\n)\n"
    path.write_text(txt)


def write_sample_list(path: Path, samples: list[str]):
    """
    Write one sample name per line, in exactly the same order as rows
    in the corresponding Glik file.
    """
    with path.open("w") as fh:
        for s in samples:
            fh.write(f"{s}\n")


def write_locus_list(path: Path, loci: list[str]):
    """
    Write one locus ID per line, in exactly the same order as columns
    in the corresponding Glik files.
    Uses VCF ID when present; otherwise falls back to contig_pos.
    """
    with path.open("w") as fh:
        for loc in loci:
            fh.write(f"{loc}\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True, type=Path)
    ap.add_argument("--popmap", required=True, type=Path)
    ap.add_argument("--p0", required=True)
    ap.add_argument("--p1", required=True)
    ap.add_argument("--out_prefix", required=True)
    ap.add_argument("--af_diff_min", type=float, default=0.0)
    ap.add_argument("--error_rate", type=float, default=0.01)
    args = ap.parse_args()

    if not (0.0 <= args.error_rate < 0.5):
        sys.exit("error_rate must be between 0 and <0.5")

    Path(args.out_prefix).parent.mkdir(parents=True, exist_ok=True)

    popmap = load_popmap(args.popmap)
    if args.p0 == args.p1 or args.p0 not in popmap.values() or args.p1 not in popmap.values():
        sys.exit("Invalid parental labels")

    # --- locus filtering EXACTLY like triangle_metrics.py ---
    samp, loci, g_main = vcf_matrix(args.vcf)
    mask_p0 = np.array([popmap[s] == args.p0 for s in samp])
    mask_p1 = np.array([popmap[s] == args.p1 for s in samp])

    het = (g_main == 1).astype(float)
    hom2 = (g_main == 2).astype(float)
    af0 = (het[mask_p0].sum(0) + 2 * hom2[mask_p0].sum(0)) / (2 * np.isfinite(g_main[mask_p0]).sum(0))
    af1 = (het[mask_p1].sum(0) + 2 * hom2[mask_p1].sum(0)) / (2 * np.isfinite(g_main[mask_p1]).sum(0))

    diff = np.abs(af0 - af1)
    keep_user = diff >= args.af_diff_min

    kept_idx = np.where(keep_user)[0].astype(int)
    kept_loci = [loci[i] for i in kept_idx]

    # groups keep VCF header order
    group_samples = {
        "P0": [s for s in samp if popmap.get(s) == args.p0],
        "P1": [s for s in samp if popmap.get(s) == args.p1],
        "ADMIX": [s for s in samp if popmap.get(s) == "ADMIX"],
    }

    if len(group_samples["P0"]) == 0 or len(group_samples["P1"]) == 0:
        sys.exit("No samples found for P0 or P1 among VCF samples (check popmap labels/sample IDs).")

    # Pre-allocate: per group, 4 matrices (A,C,G,T), rows=inds, cols=kept loci
    mats = {}
    for gname, gsamps in group_samples.items():
        mats[gname] = {b: np.full((len(gsamps), len(kept_loci)), 0.25, dtype=np.float32) for b in BASES}

    # Fill by streaming VCF; only kept loci get filled; others ignored
    fh = pysam.VariantFile(str(args.vcf))
    kept_col = -1
    for j, rec in enumerate(fh.fetch()):
        if not keep_user[j]:
            continue
        kept_col += 1

        ref_base = rec.ref.upper()
        alt_base = rec.alts[0].upper() if rec.alts and len(rec.alts) >= 1 else None
        valid_site = (ref_base in BASE_TO_I) and (alt_base in BASE_TO_I)

        for gname, gsamps in group_samples.items():
            for row_i, sname in enumerate(gsamps):
                call = rec.samples[sname]
                gt = call.get("GT", None)
                if valid_site:
                    p = pseudo_glik_for_gt(gt, ref_base, alt_base, args.error_rate)
                else:
                    p = np.array([0.25, 0.25, 0.25, 0.25], dtype=np.float32)
                for bi, b in enumerate(BASES):
                    mats[gname][b][row_i, kept_col] = p[bi]

    prefix = Path(args.out_prefix)

    # Write the 3 Glik files
    write_glik_r(prefix.with_name(prefix.name + "_GlikP0.txt"), mats["P0"])
    write_glik_r(prefix.with_name(prefix.name + "_GlikP1.txt"), mats["P1"])
    write_glik_r(prefix.with_name(prefix.name + "_GlikADMIX.txt"), mats["ADMIX"])

    # Write matching sample-name files, one sample per line, row order
    write_sample_list(prefix.with_name(prefix.name + "_GlikP0_samples.txt"), group_samples["P0"])
    write_sample_list(prefix.with_name(prefix.name + "_GlikP1_samples.txt"), group_samples["P1"])
    write_sample_list(prefix.with_name(prefix.name + "_GlikADMIX_samples.txt"), group_samples["ADMIX"])

    # Write locus-order file, one locus per line, column order
    write_locus_list(prefix.with_name(prefix.name + "_locus_order.txt"), kept_loci)

    sys.stderr.write(
        f"[make_glik] loci_total={len(loci)} loci_kept={len(kept_loci)} af_diff_min={args.af_diff_min} error_rate={args.error_rate}\n"
        f"[make_glik] nP0={len(group_samples['P0'])} nP1={len(group_samples['P1'])} nADMIX={len(group_samples['ADMIX'])}\n"
        f"[make_glik] sample_lists_written="
        f"{prefix.name}_GlikP0_samples.txt,"
        f"{prefix.name}_GlikP1_samples.txt,"
        f"{prefix.name}_GlikADMIX_samples.txt\n"
        f"[make_glik] locus_order_written={prefix.name}_locus_order.txt\n"
    )


if __name__ == "__main__":
    main()
