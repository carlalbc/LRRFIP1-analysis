"""
Main figure: 4-panel fibroblast LRRFIP1 isoform composite
Top row: Salmon short-read (A: stacked bars, B: Aitchison PCA)
Bottom row: IsoQuant long-read (C: stacked bars, D: Aitchison PCA)
Bars ordered by generation: Gen II (Mom WT, Dad affected) | Gen III (Son WT, Son affected)
Journal specs: 600 dpi, 8pt font, max 18 cm width, max 17.25 cm height, TIFF
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from sklearn.decomposition import PCA
import os, warnings
warnings.filterwarnings("ignore")

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Helvetica", "Arial", "DejaVu Sans"]
plt.rcParams["font.size"] = 8

OUTDIR = "./data/salmon_quant_iPSC"

ISO_ORDER = ["Isoform 1", "Isoform 2", "Isoform 3", "Isoform 4", "Isoform 5", "Isoform 6"]
ENST_COLS = {
    "ENST00000308482.14": "Isoform 1",
    "ENST00000863709.1":  "Isoform 2",
    "ENST00000392000.4":  "Isoform 3",
    "ENST00000244815.9":  "Isoform 4",
    "ENST00000289175.10": "Isoform 5",
    "ENST00000478958.5":  "Isoform 6",
}
ISO_COLORS = {
    "Isoform 1": "#2d004b",
    "Isoform 2": "#5e4fa2",
    "Isoform 3": "#9e6ebd",
    "Isoform 4": "#d94f8e",
    "Isoform 5": "#f47d4d",
    "Isoform 6": "#fdbe85",
    "Other":     "#d9d9d9",
}

# Generation order: Gen II first, Gen III second; WT then Affected within each
INDIV_ORDER = ["II-8.5", "II-8", "III-13", "III-12"]
INDIV_GENO = {"III-13": "WT", "II-8.5": "WT", "III-12": "Affected", "II-8": "Affected"}
INDIV_GENERATION = {"II-8.5": "Gen. II", "II-8": "Gen. II", "III-13": "Gen. III", "III-12": "Gen. III"}

# Dot colors: match COL6A3 style (solid grey/blue, high opacity)
DOT_ALPHA = 0.9
DOT_PALETTE = {"WT": "#888888", "Affected": "#2874A6"}

# ══════════════════════════════════════════════════
# Load SHORT-READ data (Salmon)
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
tpm["genotype"] = tpm["individual"].map(INDIV_GENO)
tpm["total"] = tpm[ISO_ORDER].sum(axis=1)
for iso in ISO_ORDER:
    tpm[f"{iso}_frac"] = tpm[iso] / tpm["total"]

sr_indiv_fracs = {}
for ind in INDIV_ORDER:
    sub = tpm[tpm["individual"] == ind]
    sr_indiv_fracs[ind] = {iso: sub[f"{iso}_frac"].mean() for iso in ISO_ORDER}

# ══════════════════════════════════════════════════
# Load LONG-READ data (IsoQuant PacBio FLNC)
# ══════════════════════════════════════════════════
lr_counts = {
    "Isoform 1": {"III-13": 0,   "II-8.5": 0,   "III-12": 0,   "II-8": 0},
    "Isoform 2": {"III-13": 111, "II-8.5": 112,  "III-12": 92,  "II-8": 165},
    "Isoform 3": {"III-13": 0,   "II-8.5": 3,    "III-12": 0,   "II-8": 0},
    "Isoform 4": {"III-13": 159, "II-8.5": 119,  "III-12": 87,  "II-8": 82},
    "Isoform 5": {"III-13": 346, "II-8.5": 213,  "III-12": 189, "II-8": 236},
    "Isoform 6": {"III-13": 29,  "II-8.5": 17,   "III-12": 31,  "II-8": 16},
    "Other":     {"III-13": 37,  "II-8.5": 44,   "III-12": 40,  "II-8": 37},
}

all_lr_isos = ISO_ORDER + ["Other"]
lr_totals = {ind: sum(lr_counts[iso][ind] for iso in all_lr_isos) for ind in INDIV_ORDER}
lr_fracs = {}
for ind in INDIV_ORDER:
    lr_fracs[ind] = {iso: lr_counts[iso][ind] / lr_totals[ind] if lr_totals[ind] > 0 else 0 for iso in all_lr_isos}

# ══════════════════════════════════════════════════
# FIGURE: 2×2 layout
# ══════════════════════════════════════════════════
fig_w = 18 / 2.54   # 18 cm -> inches
fig_h = 16 / 2.54   # 16 cm -> inches

fig = plt.figure(figsize=(fig_w, fig_h))
gs = fig.add_gridspec(2, 2, width_ratios=[1.2, 1],
                      left=0.10, right=0.95, top=0.90, bottom=0.10,
                      wspace=0.45, hspace=0.55)

def plot_stacked_bars(ax, indiv_fracs, isos, panel_label, ylabel, show_legend=False):
    x = np.arange(len(INDIV_ORDER))
    bar_w = 0.65
    bottoms = np.zeros(len(INDIV_ORDER))
    for iso in isos:
        vals = [indiv_fracs[ind].get(iso, 0) for ind in INDIV_ORDER]
        ax.bar(x, vals, bar_w, bottom=bottoms, color=ISO_COLORS[iso],
               edgecolor="white", linewidth=0.3, label=iso)
        bottoms += vals

    # Generation separator
    ax.axvline(1.5, color="black", linewidth=1, linestyle="-")

    ax.set_xticks(x)
    ax.set_xticklabels(INDIV_ORDER, fontsize=7)
    ax.set_ylabel(ylabel, fontsize=8)
    ax.set_ylim(0, 1.0)
    ax.set_xlim(-0.5, len(INDIV_ORDER) - 0.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    # Generation headers
    ax.text(0.5, 1.05, "Gen. II", ha="center", va="bottom", fontsize=8,
            fontweight="bold", transform=ax.get_xaxis_transform())
    ax.text(2.5, 1.05, "Gen. III", ha="center", va="bottom", fontsize=8,
            fontweight="bold", transform=ax.get_xaxis_transform())

    ax.text(-0.15, 1.08, panel_label, transform=ax.transAxes, fontsize=12, fontweight="bold", va="top")

    if show_legend:
        handles = [mpatches.Patch(facecolor=ISO_COLORS[iso], label=iso) for iso in reversed(isos)]
        ncol = 4 if "Other" in isos else 3
        ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.18),
                  fontsize=6, frameon=False, ncol=ncol,
                  handlelength=1.0, handleheight=0.7, columnspacing=0.8)


def plot_aitchison_pca(ax, indiv_fracs, pca_isos, panel_label):
    indiv_comp = np.array([[indiv_fracs[ind].get(iso, 0) for iso in pca_isos] for ind in INDIV_ORDER])
    pseudocount = 1e-6
    X = indiv_comp + pseudocount
    X = X / X.sum(axis=1, keepdims=True)
    gm = np.exp(np.mean(np.log(X), axis=1, keepdims=True))
    clr = np.log(X / gm)

    pca = PCA(n_components=2)
    scores = pca.fit_transform(clr)
    var_exp = pca.explained_variance_ratio_ * 100

    for i, ind in enumerate(INDIV_ORDER):
        c = DOT_PALETTE[INDIV_GENO[ind]]
        ax.scatter(scores[i, 0], scores[i, 1], c=c, s=70, zorder=5,
                   edgecolor="white", linewidth=0.5, alpha=DOT_ALPHA)
        ax.text(scores[i, 0] + 0.008, scores[i, 1] + 0.008, ind,
                fontsize=6.5, color=c, ha="left", va="bottom")

    ax.axhline(0, color="#cccccc", linewidth=0.5, zorder=0)
    ax.axvline(0, color="#cccccc", linewidth=0.5, zorder=0)
    ax.set_xlabel(f"PC1 ({var_exp[0]:.1f}% var)", fontsize=8)
    ax.set_ylabel(f"PC2 ({var_exp[1]:.1f}% var)", fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    legend_handles = [
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor="#888888",
                   markersize=7, alpha=DOT_ALPHA, label="WT"),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor="#2874A6",
                   markersize=7, alpha=DOT_ALPHA, label="Affected"),
    ]
    ax.legend(handles=legend_handles, fontsize=6.5, loc="best", frameon=True, framealpha=0.9)
    ax.text(-0.15, 1.08, panel_label, transform=ax.transAxes, fontsize=12, fontweight="bold", va="top")

    # Compute Aitchison distance
    wt_idx = [i for i, ind in enumerate(INDIV_ORDER) if INDIV_GENO[ind] == "WT"]
    aff_idx = [i for i, ind in enumerate(INDIV_ORDER) if INDIV_GENO[ind] == "Affected"]
    wt_centroid = clr[wt_idx].mean(axis=0)
    aff_centroid = clr[aff_idx].mean(axis=0)
    ait_dist = np.sqrt(np.sum((wt_centroid - aff_centroid)**2))
    return ait_dist, var_exp

# Panel A: Short-read stacked bars
ax_a = fig.add_subplot(gs[0, 0])
plot_stacked_bars(ax_a, sr_indiv_fracs, ISO_ORDER, "A", "Mean isoform fraction", show_legend=False)

# Panel B: Short-read Aitchison PCA
ax_b = fig.add_subplot(gs[0, 1])
sr_pca_isos = ["Isoform 2", "Isoform 4", "Isoform 5", "Isoform 6"]
sr_ait, sr_var = plot_aitchison_pca(ax_b, sr_indiv_fracs, sr_pca_isos, "B")

# Panel C: Long-read stacked bars
ax_c = fig.add_subplot(gs[1, 0])
plot_stacked_bars(ax_c, lr_fracs, ISO_ORDER + ["Other"], "C",
                  "Isoform fraction (FLNC reads)", show_legend=True)

# Panel D: Long-read Aitchison PCA
ax_d = fig.add_subplot(gs[1, 1])
lr_pca_isos = ["Isoform 2", "Isoform 4", "Isoform 5", "Isoform 6"]
lr_ait, lr_var = plot_aitchison_pca(ax_d, lr_fracs, lr_pca_isos, "D")

# Title removed — info goes in figure caption

# Save in multiple formats
for ext in ["png", "pdf", "svg"]:
    fig.savefig(os.path.join(OUTDIR, f"main_figure_fibroblast_composite.{ext}"),
                dpi=600, bbox_inches="tight", facecolor="white")
    print(f"Saved: main_figure_fibroblast_composite.{ext}")
# TIFF separately with pil_kwargs
fig.savefig(os.path.join(OUTDIR, "main_figure_fibroblast_composite.tiff"),
            dpi=600, bbox_inches="tight", facecolor="white",
            pil_kwargs={"compression": "tiff_lzw"})
print("Saved: main_figure_fibroblast_composite.tiff")

plt.close()

# Print stats
print(f"\n=== Short-read (Salmon) ===")
print(f"Aitchison distance: {sr_ait:.4f}")
print(f"PC1: {sr_var[0]:.1f}%, PC2: {sr_var[1]:.1f}%")

print(f"\n=== Long-read (IsoQuant) ===")
print(f"Aitchison distance: {lr_ait:.4f}")
print(f"PC1: {lr_var[0]:.1f}%, PC2: {lr_var[1]:.1f}%")

for label, fracs_dict in [("Short-read", sr_indiv_fracs), ("Long-read", lr_fracs)]:
    print(f"\n=== {label} individual fractions ===")
    for ind in INDIV_ORDER:
        print(f"\n{ind} ({INDIV_GENO[ind]}, {INDIV_GENERATION[ind]}):")
        for iso in ISO_ORDER:
            print(f"  {iso}: {fracs_dict[ind].get(iso, 0)*100:.1f}%")
        if "Other" in fracs_dict[ind]:
            print(f"  Other: {fracs_dict[ind]['Other']*100:.1f}%")
