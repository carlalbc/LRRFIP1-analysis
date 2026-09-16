"""
Supplementary figure: iPSC LRRFIP1 isoform composite
Panel A: stacked bars, Panel B: Aitchison PCA (all 6 isoforms)
Journal specs: 600 dpi, 8pt font, max 18 cm width, TIFF
Dot colors matched to main figure (grey WT / blue Hetero, high opacity)
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from sklearn.decomposition import PCA
import os, glob, warnings
warnings.filterwarnings("ignore")

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Helvetica", "Arial", "DejaVu Sans"]
plt.rcParams["font.size"] = 8

IPSC_ROOT = "./data/salmon_quant_iPSC"
OUTDIR = os.path.join(IPSC_ROOT, "Final_results_30.06.26/figures_and_tables/Isoform_analyses/supplementary/figures")

ENST_TO_ISO = {
    "ENST00000308482.14": "Isoform 1",
    "ENST00000863709.1":  "Isoform 2",
    "ENST00000392000.4":  "Isoform 3",
    "ENST00000244815.9":  "Isoform 4",
    "ENST00000289175.10": "Isoform 5",
    "ENST00000478958.5":  "Isoform 6",
}
ISO_ORDER = ["Isoform 1", "Isoform 2", "Isoform 3", "Isoform 4", "Isoform 5", "Isoform 6"]
keep = list(ENST_TO_ISO.keys())

ISO_COLORS = {
    "Isoform 1": "#2d004b",
    "Isoform 2": "#5e4fa2",
    "Isoform 3": "#9e6ebd",
    "Isoform 4": "#d94f8e",
    "Isoform 5": "#f47d4d",
    "Isoform 6": "#fdbe85",
}

DOT_ALPHA = 0.9
DOT_PALETTE = {"WT": "#888888", "Hetero": "#2874A6"}

# ── Load data ──
rows = []
for qf in sorted(glob.glob(os.path.join(IPSC_ROOT, "*/quant.sf"))):
    sample = os.path.basename(os.path.dirname(qf))
    df = pd.read_csv(qf, sep="\t")
    sub = df[df["Name"].isin(keep)].set_index("Name")
    for enst in keep:
        if enst in sub.index:
            rows.append({"sample": sample, "Transcript_ID": enst,
                         "Isoform": ENST_TO_ISO[enst], "TPM": sub.loc[enst, "TPM"]})

data = pd.DataFrame(rows)
data["genotype"] = data["sample"].apply(lambda s: "Hetero" if s.startswith("hetero") else "WT")
data["clone"] = data["sample"].apply(lambda s: s.split("_Repl")[0])

sample_totals = data.groupby("sample")["TPM"].sum()
data["LRRFIP1_total"] = data["sample"].map(sample_totals)
data["fraction"] = np.where(data["LRRFIP1_total"] > 0, data["TPM"] / data["LRRFIP1_total"], 0)

# Exclude replicate 2 samples (batch QC)
good_mask = ~data["sample"].str.contains("Repl2")
n_excluded = data[~good_mask]["sample"].nunique()
print(f"Excluding {n_excluded} replicate 2 samples (batch QC)")
data_good = data[good_mask].copy()

clone_fracs = data_good.groupby(["clone", "genotype", "Isoform"]).agg(
    mean_frac=("fraction", "mean"), sd_frac=("fraction", "std"),
).reset_index()
clone_fracs["sd_frac"] = clone_fracs["sd_frac"].fillna(0)

# ── Figure ──
fig_w = 18 / 2.54
fig_h = 10 / 2.54

fig = plt.figure(figsize=(fig_w, fig_h))
gs = fig.add_gridspec(1, 2, width_ratios=[1.2, 1], left=0.08, right=0.95,
                      top=0.82, bottom=0.18, wspace=0.45)

# ── Panel A: Stacked bars ──
ax_a = fig.add_subplot(gs[0, 0])

wt_clones = sorted([c for c in clone_fracs["clone"].unique() if c.startswith("wt")])
het_clones = sorted([c for c in clone_fracs["clone"].unique() if c.startswith("hetero")])
clone_order = wt_clones + het_clones
clone_labels = {"wt_1": "WT-1", "wt_2": "WT-2", "wt_3": "WT-3",
                "hetero_1": "Het-1", "hetero_2": "Het-2", "hetero_3": "Het-3"}

x = np.arange(len(clone_order))
bar_w = 0.65

bottoms = np.zeros(len(clone_order))
for iso in ISO_ORDER:
    vals = []
    for c in clone_order:
        row = clone_fracs[(clone_fracs["clone"] == c) & (clone_fracs["Isoform"] == iso)]
        vals.append(row["mean_frac"].values[0] if len(row) > 0 else 0)
    ax_a.bar(x, vals, bar_w, bottom=bottoms, color=ISO_COLORS[iso],
             edgecolor="white", linewidth=0.3, label=iso)
    bottoms += vals

sep_x = len(wt_clones) - 0.5
ax_a.axvline(sep_x, color="black", linewidth=1, linestyle="-")

ax_a.set_xticks(x)
ax_a.set_xticklabels([clone_labels.get(c, c) for c in clone_order], fontsize=7)
ax_a.set_ylabel("Mean isoform fraction", fontsize=8)
ax_a.set_ylim(0, 1.0)
ax_a.set_xlim(-0.5, len(clone_order) - 0.5)
ax_a.spines["top"].set_visible(False)
ax_a.spines["right"].set_visible(False)

ax_a.text(np.mean(range(len(wt_clones))), 1.02, "$\\it{LRRFIP1}$\nWT",
          ha="center", va="bottom", fontsize=7, transform=ax_a.get_xaxis_transform())
ax_a.text(len(wt_clones) + np.mean(range(len(het_clones))), 1.02,
          "$\\it{LRRFIP1}$\nc.96+34737C>T",
          ha="center", va="bottom", fontsize=7, color="#2874A6",
          transform=ax_a.get_xaxis_transform())

ax_a.text(-0.12, 1.08, "A", transform=ax_a.transAxes, fontsize=12, fontweight="bold", va="top")

handles = [mpatches.Patch(facecolor=ISO_COLORS[iso], label=iso) for iso in reversed(ISO_ORDER)]
ax_a.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.12),
            fontsize=6, frameon=False, ncol=3,
            handlelength=1.0, handleheight=0.7, columnspacing=1.0)

# Annotate n= per clone (all should be n=2 now)
n_good_per_clone = data_good.groupby("clone")["sample"].nunique()
for i, c in enumerate(clone_order):
    n = n_good_per_clone.get(c, 0)
    ax_a.text(i, -0.10, f"n={n}", ha="center", va="top", fontsize=6,
              color="#555555", transform=ax_a.get_xaxis_transform())

# ── Panel B: Aitchison PCA (all 6 isoforms) ──
ax_b = fig.add_subplot(gs[0, 1])

comp_data = []
clone_meta = []
for c in clone_order:
    row = {}
    for iso in ISO_ORDER:
        match = clone_fracs[(clone_fracs["clone"] == c) & (clone_fracs["Isoform"] == iso)]
        row[iso] = match["mean_frac"].values[0] if len(match) > 0 else 0
    comp_data.append(row)
    geno = "WT" if c.startswith("wt") else "Hetero"
    clone_meta.append({"clone": c, "genotype": geno, "label": clone_labels.get(c, c)})

comp_df = pd.DataFrame(comp_data, index=clone_order)
meta_df = pd.DataFrame(clone_meta)

PCA_ISOS = ISO_ORDER
comp_df_pca = comp_df[PCA_ISOS]
pseudocount = 1e-6
X = comp_df_pca.values + pseudocount
X = X / X.sum(axis=1, keepdims=True)
gm = np.exp(np.mean(np.log(X), axis=1, keepdims=True))
clr = np.log(X / gm)

pca = PCA(n_components=2)
scores = pca.fit_transform(clr)
var_exp = pca.explained_variance_ratio_ * 100

for i, row in meta_df.iterrows():
    c = DOT_PALETTE[row["genotype"]]
    ax_b.scatter(scores[i, 0], scores[i, 1], c=c, s=70, zorder=5,
                 edgecolor="white", linewidth=0.5, alpha=DOT_ALPHA)
    ax_b.text(scores[i, 0] + 0.008, scores[i, 1] + 0.008, row["label"],
              fontsize=6.5, color=c, ha="left", va="bottom")

ax_b.axhline(0, color="#cccccc", linewidth=0.5, zorder=0)
ax_b.axvline(0, color="#cccccc", linewidth=0.5, zorder=0)
ax_b.set_xlabel(f"PC1 ({var_exp[0]:.1f}% var)", fontsize=8)
ax_b.set_ylabel(f"PC2 ({var_exp[1]:.1f}% var)", fontsize=8)
ax_b.spines["top"].set_visible(False)
ax_b.spines["right"].set_visible(False)

legend_handles = [
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor="#888888",
               markersize=7, alpha=DOT_ALPHA, label="WT"),
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor="#2874A6",
               markersize=7, alpha=DOT_ALPHA, label="Hetero"),
]
ax_b.legend(handles=legend_handles, fontsize=6.5, loc="best", frameon=True, framealpha=0.9)

ax_b.text(-0.12, 1.08, "B", transform=ax_b.transAxes, fontsize=12, fontweight="bold", va="top")

for ext in ["png", "pdf", "svg"]:
    fig.savefig(os.path.join(OUTDIR, f"supp_figure_iPSC_composite.{ext}"),
                dpi=600, bbox_inches="tight", facecolor="white")
    print(f"Saved: supp_figure_iPSC_composite.{ext}")
fig.savefig(os.path.join(OUTDIR, "supp_figure_iPSC_composite.tiff"),
            dpi=600, bbox_inches="tight", facecolor="white",
            pil_kwargs={"compression": "tiff_lzw"})
print("Saved: supp_figure_iPSC_composite.tiff")

plt.close()

# Stats
print(f"\nAitchison distance (WT vs Hetero, all 6 isoforms):")
wt_cidx = [i for i, r in meta_df.iterrows() if r["genotype"] == "WT"]
het_cidx = [i for i, r in meta_df.iterrows() if r["genotype"] == "Hetero"]
wt_centroid = clr[wt_cidx].mean(axis=0)
het_centroid = clr[het_cidx].mean(axis=0)
ait_dist = np.sqrt(np.sum((wt_centroid - het_centroid)**2))
print(f"  {ait_dist:.4f}")
print(f"PC1: {var_exp[0]:.1f}%, PC2: {var_exp[1]:.1f}%")
print(f"\nSamples used: {data_good['sample'].nunique()}")
print(f"Samples excluded: {n_excluded}")
