"""
Main figure: Combined pairwise isoform delta comparisons.
Left two panels: short-read (Salmon), right two panels: long-read (PacBio FLNC).
Only generational comparisons:
  - Gen II: Father (affected) vs Mother (WT)
  - Gen III: Son (affected) vs Son (WT)
Journal specs: 600 dpi, 8pt font, max 18 cm width, TIFF
"""
import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Helvetica", "Arial", "DejaVu Sans"]
matplotlib.rcParams["font.size"] = 8

OUT_DIR = os.path.dirname(os.path.abspath(__file__)) if "__file__" in dir() else "."

ENST_COLS = {
    "ENST00000308482.14": "Isoform 1",
    "ENST00000863709.1":  "Isoform 2",
    "ENST00000392000.4":  "Isoform 3",
    "ENST00000244815.9":  "Isoform 4",
    "ENST00000289175.10": "Isoform 5",
    "ENST00000478958.5":  "Isoform 6",
}
ISO_ORDER = ["Isoform 1", "Isoform 2", "Isoform 3", "Isoform 4", "Isoform 5", "Isoform 6"]

iso_colors = {
    "Isoform 1": "#2d004b",
    "Isoform 2": "#5e4fa2",
    "Isoform 3": "#9e6ebd",
    "Isoform 4": "#d94f8e",
    "Isoform 5": "#f47d4d",
    "Isoform 6": "#fdbe85",
    "Other":     "#d9d9d9",
}

INDIV_ORDER = ["II-8.5", "II-8", "III-13", "III-12"]

# ══════════════════════════════════════════════════
# SHORT-READ data (Salmon TPM)
# ══════════════════════════════════════════════════
tpm = pd.read_csv("./data/salmon_quant_fibroblast/lrrfip1_iso6_summary/LRRFIP1_isoforms_TPM.tsv", sep="\t")
tpm = tpm.set_index("sample")
tpm.columns = [ENST_COLS[c] for c in tpm.columns]

sample_to_individual = {}
for s in tpm.index:
    if "33016" in s: sample_to_individual[s] = "III-13"
    elif "33017" in s: sample_to_individual[s] = "II-8.5"
    elif "33018" in s: sample_to_individual[s] = "III-12"
    elif "33019" in s: sample_to_individual[s] = "II-8"

tpm["individual"] = tpm.index.map(sample_to_individual)
tpm["total"] = tpm[ISO_ORDER].sum(axis=1)
for iso in ISO_ORDER:
    tpm[f"{iso}_frac"] = tpm[iso] / tpm["total"]

sr_fracs = {}
for ind in INDIV_ORDER:
    sub = tpm[tpm["individual"] == ind]
    sr_fracs[ind] = {iso: sub[f"{iso}_frac"].mean() * 100 for iso in ISO_ORDER}

# ══════════════════════════════════════════════════
# LONG-READ data (PacBio FLNC counts via IsoQuant)
# ══════════════════════════════════════════════════
# Indices: Son_WT=III-13, Mother_WT=II-8.5, Son_aff=III-12, Father_aff=II-8
lr_raw = {
    "ENST00000289175.10": {"III-13": 346, "II-8.5": 213, "III-12": 189, "II-8": 236},
    "ENST00000863709.1":  {"III-13": 111, "II-8.5": 112, "III-12":  92, "II-8": 165},
    "ENST00000244815.9":  {"III-13": 159, "II-8.5": 119, "III-12":  87, "II-8":  82},
    "ENST00000478958.5":  {"III-13":  29, "II-8.5":  17, "III-12":  31, "II-8":  16},
    "ENST00000863712.1":  {"III-13":  11, "II-8.5":   7, "III-12":   6, "II-8":  12},
    "ENST00000863720.1":  {"III-13":   6, "II-8.5":   2, "III-12":   9, "II-8":   9},
    "ENST00000863726.1":  {"III-13":   8, "II-8.5":  10, "III-12":   7, "II-8":   1},
    "ENST00000970885.1":  {"III-13":   3, "II-8.5":   3, "III-12":   6, "II-8":   4},
    "ENST00000970916.1":  {"III-13":   2, "II-8.5":   5, "III-12":   3, "II-8":   4},
}
# Add remaining minor isoforms
for enst, counts in [
    ("ENST00000933070.1",  {"III-13": 2, "II-8.5": 3, "III-12": 0, "II-8": 1}),
    ("ENST00000392000.4",  {"III-13": 0, "II-8.5": 3, "III-12": 0, "II-8": 0}),
    ("ENST00000474195.1",  {"III-13": 0, "II-8.5": 3, "III-12": 0, "II-8": 1}),
    ("ENST00000970891.1",  {"III-13": 0, "II-8.5": 3, "III-12": 0, "II-8": 1}),
    ("ENST00000863727.1",  {"III-13": 2, "II-8.5": 0, "III-12": 1, "II-8": 0}),
    ("ENST00000863729.1",  {"III-13": 0, "II-8.5": 2, "III-12": 2, "II-8": 0}),
    ("ENST00000970900.1",  {"III-13": 0, "II-8.5": 1, "III-12": 1, "II-8": 1}),
    ("ENST00000970883.1",  {"III-13": 2, "II-8.5": 0, "III-12": 0, "II-8": 0}),
    ("ENST00000970923.1",  {"III-13": 0, "II-8.5": 0, "III-12": 2, "II-8": 0}),
    ("ENST00000863713.1",  {"III-13": 0, "II-8.5": 0, "III-12": 2, "II-8": 0}),
    ("ENST00000970914.1",  {"III-13": 0, "II-8.5": 0, "III-12": 1, "II-8": 0}),
    ("ENST00000933072.1",  {"III-13": 1, "II-8.5": 0, "III-12": 0, "II-8": 1}),
    ("ENST00000863711.1",  {"III-13": 0, "II-8.5": 0, "III-12": 0, "II-8": 1}),
    ("ENST00000933071.1",  {"III-13": 0, "II-8.5": 1, "III-12": 0, "II-8": 0}),
    ("ENST00000863737.1",  {"III-13": 0, "II-8.5": 1, "III-12": 0, "II-8": 0}),
    ("ENST00000465870.5",  {"III-13": 0, "II-8.5": 0, "III-12": 0, "II-8": 1}),
]:
    lr_raw[enst] = counts

iso_map = {
    "ENST00000308482.14": "Isoform 1",
    "ENST00000863709.1":  "Isoform 2",
    "ENST00000392000.4":  "Isoform 3",
    "ENST00000244815.9":  "Isoform 4",
    "ENST00000289175.10": "Isoform 5",
    "ENST00000478958.5":  "Isoform 6",
}

lr_grouped = {}
for enst, counts in lr_raw.items():
    name = iso_map.get(enst, "Other")
    if name not in lr_grouped:
        lr_grouped[name] = {ind: 0 for ind in INDIV_ORDER}
    for ind in INDIV_ORDER:
        lr_grouped[name][ind] += counts.get(ind, 0)

lr_display = [n for n in ISO_ORDER + ["Other"] if n in lr_grouped and any(lr_grouped[n][ind] > 0 for ind in INDIV_ORDER)]
lr_totals = {ind: sum(lr_grouped[iso][ind] for iso in lr_display) for ind in INDIV_ORDER}
lr_fracs = {}
for ind in INDIV_ORDER:
    lr_fracs[ind] = {iso: lr_grouped[iso][ind] / lr_totals[ind] * 100 if lr_totals[ind] > 0 else 0 for iso in lr_display}

# ══════════════════════════════════════════════════
# Combined figure: 4 panels (2 short-read, 2 long-read)
# ══════════════════════════════════════════════════
comparisons = [
    ("Gen. II\nFather (aff) vs Mother (WT)", "II-8", "II-8.5"),
    ("Gen. III\nSon (aff) vs Son (WT)", "III-12", "III-13"),
]

fig_w = 18 / 2.54
fig_h = 9 / 2.54

fig, axes = plt.subplots(1, 4, figsize=(fig_w, fig_h))
fig.subplots_adjust(left=0.12, right=0.97, top=0.78, bottom=0.14, wspace=0.45)

# Short-read panels (0, 1)
n_sr = len(ISO_ORDER)
y_sr = np.arange(n_sr) * 1.0

for ci, (title, aff_id, wt_id) in enumerate(comparisons):
    ax = axes[ci]
    ax.set_title(title, fontsize=7, fontweight="bold", pad=6)
    deltas = [sr_fracs[aff_id][iso] - sr_fracs[wt_id][iso] for iso in ISO_ORDER]
    bar_colors = [iso_colors[iso] for iso in ISO_ORDER]
    ax.barh(y_sr, deltas, 0.6, color=bar_colors, edgecolor="white", linewidth=0.3)
    for j in range(n_sr):
        if abs(deltas[j]) > 0.3:
            ha = "left" if deltas[j] > 0 else "right"
            off = 0.3 if deltas[j] > 0 else -0.3
            ax.text(deltas[j] + off, y_sr[j], f"{deltas[j]:+.1f}", va="center", ha=ha, fontsize=6, color="#333")
    ax.axvline(0, color="#999", lw=0.5, zorder=0)
    ax.set_yticks(y_sr)
    if ci == 0:
        ax.set_yticklabels(ISO_ORDER, fontsize=7, fontweight="bold")
    else:
        ax.set_yticklabels([])
    ax.invert_yaxis()
    ax.set_xlabel("Δ (pp)", fontsize=7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    max_d = max(abs(d) for d in deltas) + 3
    ax.set_xlim(-max_d, max_d)

# Long-read panels (2, 3)
n_lr = len(lr_display)
y_lr = np.arange(n_lr) * 1.0

for ci, (title, aff_id, wt_id) in enumerate(comparisons):
    ax = axes[ci + 2]
    ax.set_title(title, fontsize=7, fontweight="bold", pad=6)
    deltas = [lr_fracs[aff_id][iso] - lr_fracs[wt_id][iso] for iso in lr_display]
    bar_colors = [iso_colors.get(iso, "#d9d9d9") for iso in lr_display]
    ax.barh(y_lr, deltas, 0.6, color=bar_colors, edgecolor="white", linewidth=0.3)
    for j in range(n_lr):
        if abs(deltas[j]) > 0.3:
            ha = "left" if deltas[j] > 0 else "right"
            off = 0.3 if deltas[j] > 0 else -0.3
            ax.text(deltas[j] + off, y_lr[j], f"{deltas[j]:+.1f}", va="center", ha=ha, fontsize=6, color="#333")
    ax.axvline(0, color="#999", lw=0.5, zorder=0)
    ax.set_yticks(y_lr)
    if ci == 0:
        ax.set_yticklabels(lr_display, fontsize=7, fontweight="bold")
    else:
        ax.set_yticklabels([])
    ax.invert_yaxis()
    ax.set_xlabel("Δ (pp)", fontsize=7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    max_d = max(abs(d) for d in deltas) + 3
    ax.set_xlim(-max_d, max_d)

# Section headers
fig.text(0.30, 0.95, "Short-read (Salmon)", ha="center", fontsize=9, fontweight="bold")
fig.text(0.75, 0.95, "Long-read (PacBio FLNC)", ha="center", fontsize=9, fontweight="bold")

for ext in ["png", "pdf", "svg"]:
    fig.savefig(os.path.join(OUT_DIR, f"main_figure_pairwise_delta_combined.{ext}"),
                dpi=600, bbox_inches="tight", facecolor="white")
    print(f"Saved: main_figure_pairwise_delta_combined.{ext}")
fig.savefig(os.path.join(OUT_DIR, "main_figure_pairwise_delta_combined.tiff"),
            dpi=600, bbox_inches="tight", facecolor="white",
            pil_kwargs={"compression": "tiff_lzw"})
print("Saved: main_figure_pairwise_delta_combined.tiff")

plt.close()
