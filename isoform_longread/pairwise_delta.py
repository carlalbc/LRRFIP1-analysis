"""
Main figure: Pairwise isoform delta comparisons (long-read FLNC)
Only generational comparisons:
  - Gen III: Son (affected) vs Son (WT)
  - Gen II: Father (affected) vs Mother (WT)
Journal specs: 600 dpi, 8pt font, max 18 cm width, TIFF
"""
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Helvetica", "Arial", "DejaVu Sans"]
matplotlib.rcParams["font.size"] = 8

OUT_DIR = os.path.dirname(os.path.abspath(__file__)) if "__file__" in dir() else "."

# FLNC read counts: [Son_WT, Mother_WT, Son_aff, Father_aff]
raw_counts = {
    "ENST00000289175.10": [346, 213, 189, 236],
    "ENST00000863709.1":  [111, 112,  92, 165],
    "ENST00000244815.9":  [159, 119,  87,  82],
    "ENST00000478958.5":  [ 29,  17,  31,  16],
    "ENST00000863712.1":  [ 11,   7,   6,  12],
    "ENST00000863720.1":  [  6,   2,   9,   9],
    "ENST00000863726.1":  [  8,  10,   7,   1],
    "ENST00000970885.1":  [  3,   3,   6,   4],
    "ENST00000970916.1":  [  2,   5,   3,   4],
    "ENST00000933070.1":  [  2,   3,   0,   1],
    "ENST00000392000.4":  [  0,   3,   0,   0],
    "ENST00000474195.1":  [  0,   3,   0,   1],
    "ENST00000970891.1":  [  0,   3,   0,   1],
    "ENST00000863727.1":  [  2,   0,   1,   0],
    "ENST00000863729.1":  [  0,   2,   2,   0],
    "ENST00000970900.1":  [  0,   1,   1,   1],
    "ENST00000970883.1":  [  2,   0,   0,   0],
    "ENST00000970923.1":  [  0,   0,   2,   0],
    "ENST00000863713.1":  [  0,   0,   2,   0],
    "ENST00000970914.1":  [  0,   0,   1,   0],
    "ENST00000933072.1":  [  1,   0,   0,   1],
    "ENST00000863711.1":  [  0,   0,   0,   1],
    "ENST00000933071.1":  [  0,   1,   0,   0],
    "ENST00000863737.1":  [  0,   1,   0,   0],
    "ENST00000465870.5":  [  0,   0,   0,   1],
}

iso_map = {
    "ENST00000308482.14": "Isoform 1",
    "ENST00000863709.1":  "Isoform 2",
    "ENST00000392000.4":  "Isoform 3",
    "ENST00000244815.9":  "Isoform 4",
    "ENST00000289175.10": "Isoform 5",
    "ENST00000478958.5":  "Isoform 6",
}

iso_colors = {
    "Isoform 1": "#2d004b",
    "Isoform 2": "#5e4fa2",
    "Isoform 3": "#9e6ebd",
    "Isoform 4": "#d94f8e",
    "Isoform 5": "#f47d4d",
    "Isoform 6": "#fdbe85",
    "Other":     "#d9d9d9",
}

# Aggregate
grouped = {}
for enst, cts in raw_counts.items():
    name = iso_map.get(enst, "Other")
    if name not in grouped:
        grouped[name] = [0, 0, 0, 0]
    for i in range(4):
        grouped[name][i] += cts[i]

# Order by isoform name (1-6, then Other); skip absent ones
display_order = ["Isoform 1", "Isoform 2", "Isoform 3", "Isoform 4", "Isoform 5", "Isoform 6", "Other"]
display_order = [n for n in display_order if n in grouped and sum(grouped[n]) > 0]

totals = [0] * 4
for cts in raw_counts.values():
    for i in range(4):
        totals[i] += cts[i]

tx_pcts = []
for name in display_order:
    cts = grouped[name]
    pcts = [cts[i] / totals[i] * 100 if totals[i] > 0 else 0 for i in range(4)]
    tx_pcts.append(pcts)

# Comparisons: generational only
# Indices: 0=Son_WT, 1=Mother_WT, 2=Son_aff, 3=Father_aff
# Gen II first, then Gen III
comparisons = [
    ("Gen. II\nFather (affected) vs Mother (WT)", 3, 1),
    ("Gen. III\nSon (affected) vs Son (WT)", 2, 0),
]

fig_w = 14 / 2.54
fig_h = 8 / 2.54

fig, axes = plt.subplots(1, 2, figsize=(fig_w, fig_h), sharey=True)
fig.subplots_adjust(left=0.18, right=0.95, top=0.85, bottom=0.12, wspace=0.35)

n_tx = len(display_order)
y_pos = np.arange(n_tx) * 1.0

for ci, (title, idx_a, idx_b) in enumerate(comparisons):
    ax = axes[ci]
    ax.set_title(title, fontsize=8, fontweight="bold", pad=6, color="black")

    deltas = [tx_pcts[j][idx_a] - tx_pcts[j][idx_b] for j in range(n_tx)]
    bar_colors = [iso_colors.get(display_order[j], "#999") for j in range(n_tx)]

    ax.barh(y_pos, deltas, 0.6, color=bar_colors, edgecolor="white", linewidth=0.3)

    for j in range(n_tx):
        if abs(deltas[j]) > 0.3:
            ha = "left" if deltas[j] > 0 else "right"
            off = 0.3 if deltas[j] > 0 else -0.3
            ax.text(deltas[j] + off, y_pos[j],
                    f"{deltas[j]:+.1f}", va="center", ha=ha, fontsize=7, color="#333")

    ax.axvline(0, color="#999", lw=0.5, zorder=0)
    ax.set_yticks(y_pos)
    if ci == 0:
        ax.set_yticklabels(display_order, fontsize=8, fontweight="bold")
    ax.invert_yaxis()
    ax.set_xlabel("Δ (pp)", fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    max_d = max(abs(d) for d in deltas) + 3
    ax.set_xlim(-max_d, max_d)

# Panel labels
axes[0].text(-0.18, 1.12, "A", transform=axes[0].transAxes, fontsize=12, fontweight="bold", va="top")
axes[1].text(-0.18, 1.12, "B", transform=axes[1].transAxes, fontsize=12, fontweight="bold", va="top")

for ext in ["png", "pdf", "svg"]:
    fig.savefig(os.path.join(OUT_DIR, f"main_figure_pairwise_delta.{ext}"),
                dpi=600, bbox_inches="tight", facecolor="white")
    print(f"Saved: main_figure_pairwise_delta.{ext}")
fig.savefig(os.path.join(OUT_DIR, "main_figure_pairwise_delta.tiff"),
            dpi=600, bbox_inches="tight", facecolor="white",
            pil_kwargs={"compression": "tiff_lzw"})
print("Saved: main_figure_pairwise_delta.tiff")

plt.close()

# Print summary
for title, idx_a, idx_b in comparisons:
    print(f"\n{title.replace(chr(10), ' ')}:")
    for j, name in enumerate(display_order):
        d = tx_pcts[j][idx_a] - tx_pcts[j][idx_b]
        print(f"  {name}: {d:+.1f} pp ({tx_pcts[j][idx_a]:.1f}% vs {tx_pcts[j][idx_b]:.1f}%)")
