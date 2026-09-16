"""
IsoQuant (PacBio FLNC) LRRFIP1 isoform composite figure
Matching the Salmon fibroblast composite format: Panel A stacked bars + Panel B Aitchison PCA
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from sklearn.decomposition import PCA
import os, warnings
warnings.filterwarnings("ignore")

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
plt.rcParams["font.size"] = 9

OUT_DIR = os.path.dirname(os.path.abspath(__file__)) if "__file__" in dir() else "."

# ── Raw FLNC IsoQuant read counts ──
counts = {
    "Isoform 1": {"III-13": 0,   "II-8.5": 0,   "III-12": 0,   "II-8": 0},
    "Isoform 2": {"III-13": 111, "II-8.5": 112,  "III-12": 92,  "II-8": 165},
    "Isoform 3": {"III-13": 0,   "II-8.5": 3,    "III-12": 0,   "II-8": 0},
    "Isoform 4": {"III-13": 159, "II-8.5": 119,  "III-12": 87,  "II-8": 82},
    "Isoform 5": {"III-13": 346, "II-8.5": 213,  "III-12": 189, "II-8": 236},
    "Isoform 6": {"III-13": 29,  "II-8.5": 17,   "III-12": 31,  "II-8": 16},
    "Other":     {"III-13": 37,  "II-8.5": 44,   "III-12": 40,  "II-8": 37},
}

ISO_ORDER = ["Isoform 1", "Isoform 2", "Isoform 3", "Isoform 4", "Isoform 5", "Isoform 6"]
ISO_COLORS = {
    "Isoform 1": "#2d004b",
    "Isoform 2": "#5e4fa2",
    "Isoform 3": "#9e6ebd",
    "Isoform 4": "#d94f8e",
    "Isoform 5": "#f47d4d",
    "Isoform 6": "#fdbe85",
    "Other":     "#d9d9d9",
}

indiv_order = ["III-13", "II-8.5", "III-12", "II-8"]
indiv_geno = {"III-13": "WT", "II-8.5": "WT", "III-12": "Affected", "II-8": "Affected"}
indiv_labels = {
    "III-13": "III-13\nhealthy", "II-8.5": "II-8.5\nhealthy",
    "III-12": "III-12\naffected", "II-8": "II-8\naffected",
}

# Compute fractions (including Other in denominator)
all_isos = ISO_ORDER + ["Other"]
totals = {ind: sum(counts[iso][ind] for iso in all_isos) for ind in indiv_order}
fracs = {}
for ind in indiv_order:
    fracs[ind] = {}
    for iso in all_isos:
        fracs[ind][iso] = counts[iso][ind] / totals[ind] if totals[ind] > 0 else 0

# ── Figure ──
fig = plt.figure(figsize=(20/2.54, 13/2.54))
gs = fig.add_gridspec(1, 2, width_ratios=[1.2, 1], left=0.08, right=0.92,
                       top=0.82, bottom=0.18, wspace=0.45)

# ── Panel A: Stacked bars ──
ax_a = fig.add_subplot(gs[0, 0])

x = np.arange(len(indiv_order))
bar_w = 0.65

bottoms = np.zeros(len(indiv_order))
for iso in ISO_ORDER + ["Other"]:
    vals = [fracs[ind][iso] for ind in indiv_order]
    ax_a.bar(x, vals, bar_w, bottom=bottoms, color=ISO_COLORS[iso],
             edgecolor="white", linewidth=0.3, label=iso)
    bottoms += vals

# Separator
ax_a.axvline(1.5, color="black", linewidth=1, linestyle="-")

ax_a.set_xticks(x)
ax_a.set_xticklabels([indiv_labels[ind] for ind in indiv_order], fontsize=7.5)
ax_a.set_ylabel("Isoform fraction (FLNC reads)", fontsize=9)
ax_a.set_ylim(0, 1.0)
ax_a.set_xlim(-0.5, len(indiv_order) - 0.5)
ax_a.spines["top"].set_visible(False)
ax_a.spines["right"].set_visible(False)

# Group headers
ax_a.text(0.5, 1.05, "$\it{LRRFIP1}$\nWT",
          ha="center", va="bottom", fontsize=8, transform=ax_a.get_xaxis_transform())
ax_a.text(2.5, 1.05, "$\it{LRRFIP1}$\nc.96+34737C>T",
          ha="center", va="bottom", fontsize=8, color="#2874A6",
          transform=ax_a.get_xaxis_transform())

ax_a.text(-0.15, 1.05, "A", transform=ax_a.transAxes, fontsize=14, fontweight="bold", va="top")

# Legend below — include Other
legend_isos = ISO_ORDER + ["Other"]
handles = [mpatches.Patch(facecolor=ISO_COLORS[iso], label=iso) for iso in reversed(legend_isos)]
ax_a.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.15),
            fontsize=6, frameon=False, ncol=4,
            handlelength=1.0, handleheight=0.7, columnspacing=0.8)

# ── Panel B: Aitchison PCA ──
ax_b = fig.add_subplot(gs[0, 1])

# CLR on 4 isoforms (Iso2/4/5/6) — excluding Iso1 (absent), Iso3 (negligible), Other
PCA_ISOS = ["Isoform 2", "Isoform 4", "Isoform 5", "Isoform 6"]
indiv_comp = np.array([[fracs[ind][iso] for iso in PCA_ISOS] for ind in indiv_order])
pseudocount = 1e-6
X = indiv_comp + pseudocount
X = X / X.sum(axis=1, keepdims=True)
gm = np.exp(np.mean(np.log(X), axis=1, keepdims=True))
clr_indiv = np.log(X / gm)

pca = PCA(n_components=2)
scores = pca.fit_transform(clr_indiv)
var_exp = pca.explained_variance_ratio_ * 100

palette = {"WT": "#888888", "Affected": "#2874A6"}
for i, ind in enumerate(indiv_order):
    c = palette[indiv_geno[ind]]
    ax_b.scatter(scores[i, 0], scores[i, 1], c=c, s=80, zorder=5,
                 edgecolor="white", linewidth=0.5, alpha=0.65)
    ax_b.text(scores[i, 0] + 0.008, scores[i, 1] + 0.008, ind,
              fontsize=7, color=c, ha="left", va="bottom")

ax_b.axhline(0, color="#cccccc", linewidth=0.5, zorder=0)
ax_b.axvline(0, color="#cccccc", linewidth=0.5, zorder=0)
ax_b.set_xlabel(f"PC1 ({var_exp[0]:.1f}% var)", fontsize=9)
ax_b.set_ylabel(f"PC2 ({var_exp[1]:.1f}% var)", fontsize=9)
ax_b.spines["top"].set_visible(False)
ax_b.spines["right"].set_visible(False)

legend_handles = [
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor="#888888", markersize=8, alpha=0.65, label="WT"),
    plt.Line2D([0], [0], marker='o', color='w', markerfacecolor="#2874A6", markersize=8, alpha=0.65, label="Affected"),
]
ax_b.legend(handles=legend_handles, fontsize=7, loc="best", frameon=True, framealpha=0.9)

ax_b.text(-0.15, 1.05, "B", transform=ax_b.transAxes, fontsize=14, fontweight="bold", va="top")

fig.suptitle("Fibroblast $\it{LRRFIP1}$ isoform composition — IsoQuant (PacBio FLNC)",
             fontsize=10, fontweight="bold", y=0.97)

for ext in ["png", "pdf", "svg"]:
    fig.savefig(os.path.join(OUT_DIR, f"isoquant_LRRFIP1_isoform_composite.{ext}"),
                dpi=300, bbox_inches="tight", facecolor="white")
    print(f"Saved: isoquant_LRRFIP1_isoform_composite.{ext}")

plt.close()

# ── Stats ──
print("\n=== Individual fractions (6 isoforms only, excl Other) ===")
for ind in indiv_order:
    print(f"\n{ind} ({indiv_geno[ind]}):")
    iso_total = sum(counts[iso][ind] for iso in ISO_ORDER)
    for iso in ISO_ORDER:
        f = counts[iso][ind] / iso_total * 100 if iso_total > 0 else 0
        print(f"  {iso}: {f:.1f}% ({counts[iso][ind]} reads)")

# Aitchison distance
wt_idx = [i for i, ind in enumerate(indiv_order) if indiv_geno[ind] == "WT"]
aff_idx = [i for i, ind in enumerate(indiv_order) if indiv_geno[ind] == "Affected"]
wt_centroid = clr_indiv[wt_idx].mean(axis=0)
aff_centroid = clr_indiv[aff_idx].mean(axis=0)
ait_dist = np.sqrt(np.sum((wt_centroid - aff_centroid)**2))
print(f"\nAitchison distance (WT vs Affected, 4 isoforms): {ait_dist:.4f}")
print(f"  (Salmon fibroblast Aitchison: 0.36 with 4 isoforms, 0.78 all 6)")
print(f"PC1: {var_exp[0]:.1f}%, PC2: {var_exp[1]:.1f}%")
