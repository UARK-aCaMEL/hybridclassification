#!/usr/bin/env python3
"""
Spatial map of NewHybrids results
---------------------------------
* Folium pie‑charts at sampling sites
* Optional GeoJSON / shapefile overlays (with per‑layer style + z‑order)
* Legend for class colours and pie‑size scale
* Hover tooltip shows N + counts per class
* Writes an optional site‑summary TSV
"""

import geopandas as gpd
import pandas as pd
import folium
from branca.element import Element
import matplotlib.pyplot as plt
import plotly.express as px
import argparse, json, io, base64, re
from pathlib import Path


# ───────────────────────── helpers ─────────────────────────
def read_nh_results(path):
    cats = ["P0", "P1", "F1", "F2", "Bx0", "Bx1"]
    cols = ["Index", "Individual"] + cats
    return pd.read_csv(path, sep=r"\s+", skiprows=1, header=None, names=cols)


def load_maps(nh_map, popmap):
    df_map = (
        pd.read_csv(nh_map, sep="\t", header=0)
        .rename(columns={"Sample": "Individual"})
    )
    df_pop = pd.read_csv(popmap, sep="\t", header=None, names=["Individual", "Population"])
    return df_map, df_pop


def load_coords(path):
    return pd.read_csv(path, sep="\t", header=None, names=["ID", "Latitude", "Longitude"])


def merge_coords(df, coords):
    """
    Attach coordinates to each row, preferring the join that yields the
    greatest number of non‑missing Latitude/Longitude values.

    1. Try an Individual‑level join.
    2. Try a Population‑level join.
    3. Keep whichever join produces MORE successful coordinate matches.
    """
    by_ind = df.merge(coords, left_on="Individual",  right_on="ID", how="left")
    by_pop = df.merge(coords, left_on="Population",  right_on="ID", how="left")

    good_ind = by_ind["Latitude"].notna().sum()
    good_pop = by_pop["Latitude"].notna().sum()

    best = by_ind if good_ind >= good_pop else by_pop
    # Drop rows that still lack coordinates so they never reach plotting
    return best.dropna(subset=["Latitude", "Longitude"])


def aggregate_sites(df_nh, df_map, df_pop, coords, base_cats, th):
    df = (
        df_nh.drop(columns="Individual")
        .merge(df_map, on="Index")
        .merge(df_pop, on="Individual", how="left")
    )
    df = merge_coords(df, coords)
    df[base_cats] = df[base_cats].apply(pd.to_numeric, errors="coerce")

    def classify(r):
        ok = r[base_cats] >= th
        return base_cats[ok.argmax()] if ok.sum() == 1 else "Unclassified"

    df["Class"] = df.apply(classify, axis=1)

    counts = (
        df.groupby(["ID", "Latitude", "Longitude", "Class"], as_index=False)
        .size()
        .rename(columns={"size": "count"})
    )
    pivot = (
        counts.pivot_table(
            index=["ID", "Latitude", "Longitude"],
            columns="Class",
            values="count",
            fill_value=0,
        )
        .reset_index()
    )

    for c in base_cats + ["Unclassified"]:
        if c not in pivot:
            pivot[c] = 0
    pivot["total"] = pivot[base_cats + ["Unclassified"]].sum(axis=1)
    return pivot


def rgb_to_hex(col):
    m = re.match(r"rgb\((\d+),\s*(\d+),\s*(\d+)\)", col)
    return f"#{int(m[1]):02x}{int(m[2]):02x}{int(m[3]):02x}" if m else col


def palette(name, n):
    if hasattr(px.colors.qualitative, name):
        seq = getattr(px.colors.qualitative, name)
    elif hasattr(px.colors.sequential, name):
        seq = px.colors.sample_colorscale(
            getattr(px.colors.sequential, name),
            [i / (n - 1) for i in range(n)],
        )
    elif hasattr(px.colors.diverging, name):
        seq = px.colors.sample_colorscale(
            getattr(px.colors.diverging, name),
            [i / (n - 1) for i in range(n)],
        )
    else:
        raise ValueError(f"palette '{name}' not found")
    if len(seq) < n:
        seq = (seq * ((n // len(seq)) + 1))[:n]
    return [rgb_to_hex(c) for c in seq]


def pie_icon(vals, cols, size_px):
    fig, ax = plt.subplots(figsize=(1, 1), dpi=size_px)
    ax.pie(vals, startangle=90, colors=cols)
    ax.axis("equal")
    buf = io.BytesIO()
    plt.savefig(buf, format="png", transparent=True, bbox_inches="tight", pad_inches=0)
    plt.close(fig)
    return "data:image/png;base64," + base64.b64encode(buf.getvalue()).decode()


def load_overlays(json_path):
    if not json_path:
        return []
    data = json.loads(Path(json_path).read_text())
    if isinstance(data, dict):
        data = [data]
    return sorted(data, key=lambda d: int(d.get("z_order", 0)))


# ───────────────────────── main ─────────────────────────
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--result", required=True)
    ap.add_argument("--result_map", required=True)
    ap.add_argument("--popmap", required=True)
    ap.add_argument("--site_coords", required=True)
    ap.add_argument("--threshold", type=float, default=0.5)
    ap.add_argument("--geo_data_json", help="JSON list of overlay layers")
    ap.add_argument("--template", help="HTML header template to prepend")
    ap.add_argument("--out", required=True, help="Output HTML file")
    ap.add_argument("--palette", default="Spectral")
    ap.add_argument("--min_pie_px", type=int, default=20)
    ap.add_argument("--max_pie_px", type=int, default=80)
    ap.add_argument(
        "--table_out",
        help="Write site summary TSV (default: <out>.tsv)",
    )
    args = ap.parse_args()

    base_cats = ["P0", "Bx0", "F1", "F2", "Bx1", "P1"]
    cats = base_cats + ["Unclassified"]

    df_sites = aggregate_sites(
        read_nh_results(args.result),
        *load_maps(args.result_map, args.popmap),
        load_coords(args.site_coords),
        base_cats,
        args.threshold,
    )

    if df_sites.empty:
        raise ValueError("No valid site coordinates found")

    # ── optional table ‑‑──────────────────────────────────
    table_path = (
        Path(args.table_out)
        if args.table_out
        else Path(args.out).with_suffix(".tsv")
    )
    cols_order = ["ID", "Latitude", "Longitude", "total"] + cats
    df_sites[cols_order].to_csv(table_path, sep="\t", index=False)
    print(f"✅ Site summary TSV saved to: {table_path}")

    # ── Folium map ─────────────────────────────────────────
    m = folium.Map(
        location=[df_sites.Latitude.mean(), df_sites.Longitude.mean()],
        zoom_start=6,
        tiles="CartoDB Positron",
        control_scale=True,
        zoom_control=False,
        width="100%",
        height="100%",
    )

    m.get_root().header.add_child(
        Element(
            """
            <style>
              html,body,#map{height:100%!important;margin:0;}
              .leaflet-control-zoomslider,
              .leaflet-control-zoom{display:none!important;}
              .leaflet-control-scale{background:rgba(255,255,255,0.9);font-size:10px;}
              #map{border:2px solid #000;}
            </style>
            """
        )
    )
    m.get_root().html.add_child(
        Element(
            '<div style="position:absolute;top:8px;right:10px;z-index:999;'
            'font-size:1.2em;">↑ N</div>'
        )
    )

    # overlays in z‑order
    for layer in load_overlays(args.geo_data_json):
        style = layer.get("style", {})
        gdf = gpd.read_file(layer["path"]).to_crs(epsg=4326)
        folium.GeoJson(
            gdf.__geo_interface__,
            name=Path(layer["path"]).stem,
            style_function=lambda _f, s=style: s,
        ).add_to(m)

    # auto‑zoom
    m.fit_bounds(
        [
            [df_sites.Latitude.min(), df_sites.Longitude.min()],
            [df_sites.Latitude.max(), df_sites.Longitude.max()],
        ]
    )

    # markers
    colors = palette(args.palette, len(base_cats)) + ["#808080"]
    max_n = df_sites["total"].max()
    min_n = df_sites["total"].min()

    for _, row in df_sites.iterrows():
        scale = (row["total"] / max_n) ** 0.5
        px_size = int(
            args.min_pie_px + (args.max_pie_px - args.min_pie_px) * scale
        )
        tooltip = (
            f"{row.ID} (N={int(row['total'])}): "
            + ", ".join(f"{c}={int(row[c])}" for c in cats)
        )
        folium.Marker(
            [row.Latitude, row.Longitude],
            icon=folium.CustomIcon(
                pie_icon([row[c] for c in cats], colors, px_size),
                icon_size=(px_size // 2, px_size // 2),
            ),
            tooltip=tooltip,
        ).add_to(m)

    # legend
    def marker_dim(n):
        return int(
            args.min_pie_px
            + (args.max_pie_px - args.min_pie_px) * ((n / max_n) ** 0.5)
        ) // 2

    size_refs = [
        (min_n, f"{min_n} sample" if min_n == 1 else f"{min_n} samples"),
        (max_n, f"{max_n} samples"),
    ]
    size_items = ""
    for n, lbl in size_refs:
        marker_d   = marker_dim(n)        # visible pie diameter
        svg_dim    = marker_d + 6         # add margin for the circle
        radius     = marker_d // 2
        size_items += (
            f'<div style="display:flex;align-items:center;margin-bottom:2px;">'
            f'  <svg width="{svg_dim}" height="{svg_dim}" style="margin-right:6px;">'
            f'    <circle cx="{svg_dim//2}" cy="{svg_dim//2}" r="{radius}" fill="#999"/></svg>{lbl}'
            f'</div>'
        )
    class_items = "".join(
        f'<div style="margin-bottom:2px;">'
        f'<i style="background:{col};display:inline-block;width:12px;'
        f'height:12px;margin-right:6px;"></i>{cat}</div>'
        for cat, col in zip(cats, colors)
    )
    legend_html = (
        '<div style="border:1px solid #ccc;padding:10px;background:#f8f8f8;'
        'font-size:0.9em;">'
        "<b>Class</b>"
        + class_items
        + '<hr style="margin:6px 0;"><b>Pie size</b>'
        + size_items
        + "</div>"
    )

    # wrap map + legend
    snippet = m._repr_html_()
    wrapped = f"""
    <div style="border:1px solid #ddd;border-radius:4px;padding:10px;margin-bottom:1em;
                display:flex;flex-wrap:wrap;">
        <div style="flex:1 1 500px;min-width:400px;">{snippet}</div>
        <div style="flex:0 0 250px;margin-left:20px;">{legend_html}</div>
    </div>
    """

    header = Path(args.template).read_text() if args.template else ""
    Path(args.out).write_text(header + wrapped)
    print(f"✅ Spatial map saved to: {args.out}")


if __name__ == "__main__":
    main()
