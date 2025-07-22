import glob, os, sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# ── 1.  auto‑detect the metrics + popmap files ────────────────────
try:
    hindex_path = glob.glob("*_hindex.tsv")[0]             # e.g. CAMANO_CAMOLI_hindex.tsv
except IndexError:
    sys.exit("❌  No *_hindex.tsv file found in this directory.")

prefix = hindex_path.rsplit("_hindex.tsv", 1)[0]
popmap_path = f"{prefix}_classification_popmap.tsv"

if not os.path.exists(popmap_path):
    sys.exit(f"❌  Expected popmap file '{popmap_path}' not found.")

print(f"→ using prefix '{prefix}'")

# ── 2.  load data ────────────────────────────────────────────────
df_metrics = pd.read_csv(hindex_path, sep="\t")       # Sample, HybridIndex, Heterozygosity, PercMissing
df_popmap  = pd.read_csv(popmap_path, sep="\t")       # Sample, Group

data = df_metrics.merge(df_popmap, on="Sample", how="left")

# ── 3.  assign colours per group ─────────────────────────────────
groups = data["Group"].fillna("Unknown").unique()
palette = plt.cm.get_cmap("tab10", len(groups))
colour  = {g: palette(i) for i, g in enumerate(groups)}

# ── 4.  draw triangle plot ───────────────────────────────────────
fig, ax = plt.subplots(figsize=(5, 4))

# outline
ax.plot([0, 0.5, 1], [0, 1, 0], color="black", lw=1)        # straight edges
x = np.linspace(0, 1, 200)
ax.plot(x, 2 * x * (1 - x), color="black", lw=1)            # curved lower bound

# points
for g in groups:
    sel = data["Group"] == g
    ax.scatter(
        data.loc[sel, "HybridIndex"],
        data.loc[sel, "Heterozygosity"],
        s=30, c=[colour[g]], label=g, edgecolors="k", linewidths=0.3
    )

ax.set_xlim(-0.05, 1.05)
ax.set_ylim(-0.05, 1.05)
ax.set_xlabel("Hybrid Index")
ax.set_ylabel("Inter‑class Heterozygosity")
ax.legend(title="Group", fontsize=8, markerscale=1.2, frameon=False)

plt.tight_layout()
plt.show()

