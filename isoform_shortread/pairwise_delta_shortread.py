"""
Supplementary figure: Pairwise isoform delta comparisons (short-read, Salmon)
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
}

# Individual mapping: sample ID substring → pedigree ID
INDIV_ORDER = ["II-8.5", "II-8", "III-13", "III-12"]
INDIV_GENO = {"III-13": "WT", "II-8.5": "WT", "III-12": "Affected", "II-8": "Affected"}

# Load short-read data (Salmon TPM)
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

# Mean isoform fraction per individual (mean of clone means)
indiv_fracs = {}
for ind in INDIV_ORDER:
    sub = tpm[tpm["individual"] == ind]
    indiv_fracs[ind] = {iso: sub[f"{iso}_frac"].mean() * 100 for iso in ISO_ORDER}

# Comparisons: generational only
# Gen II: Father (II-8, affected) vs Mother (II-8.5, WT)
# Gen III: Son (III-12, affected) vs Son (III-13, WT)
comparisons = [
    ("Gen. II\nFather (affected) vs Mother (WT)", "II-8", "II-8.5"),
    ("Gen. III\nSon (affected) vs Son (WT)", "III-12", "III-13"),
]

fig_w = 14 / 2.54
fig_h = 8 / 2.54

fig, axes = plt.subplots(1, 2, figsize=(fig_w, fig_h), sharey=True)
fig.subplots_adjust(left=0.18, right=0.95, top=0.85, bottom=0.12, wspace=0.35)

n_tx = len(ISO_ORDER)
y_pos = np.arange(n_tx) * 1.0

for ci, (title, aff_id, wt_id) in enumerate(comparisons):
    ax = axes[ci]
    ax.set_title(title, fontsize=8, fontweight="bold", pad=6, color="black")

    deltas = [indiv_fracs[aff_id][iso] - indiv_fracs[wt_id][iso] for iso in ISO_ORDER]
    bar_colors = [iso_colors[iso] for iso in ISO_ORDER]

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
        ax.set_yticklabels(ISO_ORDER, fontsize=8, fontweight="bold")
    ax.invert_yaxis()
    ax.set_xlabel("Δ (pp)", fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    max_d = max(abs(d) for d in deltas) + 3
    ax.set_xlim(-max_d, max_d)

axes[0].text(-0.18, 1.12, "A", transform=axes[0].transAxes, fontsize=12, fontweight="bold", va="top")
axes[1].text(-0.18, 1.12, "B", transform=axes[1].transAxes, fontsize=12, fontweight="bold", va="top")

for ext in ["png", "pdf", "svg"]:
    fig.savefig(os.path.join(OUT_DIR, f"supp_figure_pairwise_delta_shortread.{ext}"),
                dpi=600, bbox_inches="tight", facecolor="white")
    print(f"Saved: supp_figure_pairwise_delta_shortread.{ext}")
fig.savefig(os.path.join(OUT_DIR, "supp_figure_pairwise_delta_shortread.tiff"),
            dpi=600, bbox_inches="tight", facecolor="white",
            pil_kwargs={"compression": "tiff_lzw"})
print("Saved: supp_figure_pairwise_delta_shortread.tiff")

plt.close()

for title, aff_id, wt_id in comparisons:
    print(f"\n{title.replace(chr(10), ' ')}:")
    for iso in ISO_ORDER:
        d = indiv_fracs[aff_id][iso] - indiv_fracs[wt_id][iso]
        print(f"  {iso}: {d:+.1f} pp ({indiv_fracs[aff_id][iso]:.1f}% vs {indiv_fracs[wt_id][iso]:.1f}%)")
