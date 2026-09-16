"""
FLNC IsoQuant: Isoform distribution — WT vs Affected
Panel A: Individual samples (4 bars)
Panel B: Pooled WT vs Pooled Affected
"""
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
matplotlib.rcParams["font.size"] = 9

# ── Data (FLNC IsoQuant counts, mapped to our isoform names) ──
isoforms = {
    "Isoform5\n(ENST289175)\nGENCODE Primary": {"Son_WT": 346, "Mother_WT": 213, "Son_aff": 189, "Father_aff": 236},
    "Isoform2\n(ENST863709)\nAPPRIS P1":  {"Son_WT": 111, "Mother_WT": 112, "Son_aff": 92,  "Father_aff": 165},
    "Isoform4\n(ENST244815)\nGencode Primary": {"Son_WT": 159, "Mother_WT": 119, "Son_aff": 87,  "Father_aff": 82},
    "Isoform6\n(ENST478958)\nret. intron": {"Son_WT": 29,  "Mother_WT": 17,  "Son_aff": 31,  "Father_aff": 16},
    "Other\n(17 minor tx)":                {"Son_WT": 37,  "Mother_WT": 44,  "Son_aff": 40,  "Father_aff": 37},
}

iso_names = list(isoforms.keys())
samples_individual = ["Son_WT", "Mother_WT", "Son_aff", "Father_aff"]
sample_labels = ["Son (WT)", "Mother (WT)", "Son (aff)", "Father (aff)"]
sample_colors = ["#5B9BD5", "#2E75B6", "#E07B54", "#C0392B"]

# Compute totals per sample
totals = {}
for s in samples_individual:
    totals[s] = sum(isoforms[iso][s] for iso in iso_names)

# Compute proportions
props = {}
for iso in iso_names:
    props[iso] = {}
    for s in samples_individual:
        props[iso][s] = isoforms[iso][s] / totals[s] * 100

# Pooled
pooled_wt = {}
pooled_aff = {}
for iso in iso_names:
    pooled_wt[iso] = isoforms[iso]["Son_WT"] + isoforms[iso]["Mother_WT"]
    pooled_aff[iso] = isoforms[iso]["Son_aff"] + isoforms[iso]["Father_aff"]

total_wt = sum(pooled_wt.values())
total_aff = sum(pooled_aff.values())

pooled_wt_pct = {iso: pooled_wt[iso] / total_wt * 100 for iso in iso_names}
pooled_aff_pct = {iso: pooled_aff[iso] / total_aff * 100 for iso in iso_names}

# ── Isoform colors ──
iso_colors = ["#2E75B6", "#E07B54", "#7FB685", "#9B72AA", "#BBBBBB"]

# ── Figure ──
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(24/2.54, 12/2.54),
                                 gridspec_kw={"width_ratios": [2.2, 1]})
fig.subplots_adjust(wspace=0.35)

# ── Panel A: Individual samples — grouped bars ──
x = np.arange(len(iso_names))
bar_w = 0.18
offsets = [-1.5, -0.5, 0.5, 1.5]

for si, (s, label, col) in enumerate(zip(samples_individual, sample_labels, sample_colors)):
    vals = [props[iso][s] for iso in iso_names]
    bars = ax1.bar(x + offsets[si] * bar_w, vals, bar_w * 0.9,
                   color=col, label=label, edgecolor="white", linewidth=0.5)
    # Value labels on top
    for bar, v in zip(bars, vals):
        if v > 3:
            ax1.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.8,
                     f"{v:.0f}%", ha="center", va="bottom", fontsize=6, color=col)

ax1.set_xticks(x)
ax1.set_xticklabels([n.split("\n")[0] for n in iso_names], fontsize=8.5, fontweight="bold")

# Add subtitle lines under isoform names
for i, iso in enumerate(iso_names):
    parts = iso.split("\n")
    subtitle = "\n".join(parts[1:]) if len(parts) > 1 else ""
    if subtitle:
        ax1.text(i, -4.5, subtitle, ha="center", va="top", fontsize=6,
                 color="#666", style="italic")

ax1.set_ylabel("% of assigned LRRFIP1 reads", fontsize=9)
ax1.set_ylim(0, 62)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)
ax1.legend(fontsize=7.5, loc="upper right", framealpha=0.9, edgecolor="#ddd")
ax1.set_title("A) Individual samples", fontsize=10, fontweight="bold", loc="left", pad=8)

# ── Panel B: Pooled WT vs Affected — stacked horizontal bars ──
groups = ["WT\n(pooled)", "Affected\n(pooled)"]
bottoms_wt = 0
bottoms_aff = 0

for i, iso in enumerate(iso_names):
    short = iso.split("\n")[0]
    # WT bar
    w = pooled_wt_pct[iso]
    ax2.barh(1, w, left=bottoms_wt, color=iso_colors[i], edgecolor="white",
             linewidth=0.5, height=0.5)
    if w > 4:
        ax2.text(bottoms_wt + w/2, 1, f"{w:.0f}%", ha="center", va="center",
                 fontsize=7, fontweight="bold", color="white")
    bottoms_wt += w

    # Affected bar
    a = pooled_aff_pct[iso]
    ax2.barh(0, a, left=bottoms_aff, color=iso_colors[i], edgecolor="white",
             linewidth=0.5, height=0.5, label=short if i == 0 else "")
    if a > 4:
        ax2.text(bottoms_aff + a/2, 0, f"{a:.0f}%", ha="center", va="center",
                 fontsize=7, fontweight="bold", color="white")
    bottoms_aff += a

ax2.set_yticks([0, 1])
ax2.set_yticklabels(["Affected\n(pooled)", "WT\n(pooled)"], fontsize=9)
ax2.set_xlim(0, 108)
ax2.set_xlabel("% of assigned LRRFIP1 reads", fontsize=8)
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)
ax2.set_title("B) Pooled comparison", fontsize=10, fontweight="bold", loc="left", pad=8)

# Stacked bar legend
from matplotlib.patches import Patch
legend_patches = [Patch(facecolor=iso_colors[i], label=iso_names[i].split("\n")[0])
                  for i in range(len(iso_names))]
ax2.legend(handles=legend_patches, fontsize=7, loc="upper center",
           bbox_to_anchor=(0.5, -0.18), ncol=3,
           framealpha=0.9, edgecolor="#ddd")

# ── Annotation: highlight Isoform2 shift ──
# Arrow between WT and Affected for Isoform2 in the pooled panel
wt_iso2_start = pooled_wt_pct[iso_names[0]]  # after Isoform5
wt_iso2_mid = wt_iso2_start + pooled_wt_pct[iso_names[1]] / 2
aff_iso2_start = pooled_aff_pct[iso_names[0]]
aff_iso2_mid = aff_iso2_start + pooled_aff_pct[iso_names[1]] / 2

# Add n= labels
ax2.text(100, 1, f"n={total_wt:.0f}", ha="left", va="center", fontsize=7, color="#666",
         style="italic")
ax2.text(100, 0, f"n={total_aff:.0f}", ha="left", va="center", fontsize=7, color="#666",
         style="italic")

import os
RESULTS_DIR = os.path.dirname(os.path.abspath(__file__)) if "__file__" in dir() else "."

for ext in ("png", "pdf", "svg"):
    out = os.path.join(RESULTS_DIR, f"flnc_isoform_distribution.{ext}")
    fig.savefig(out, dpi=300, bbox_inches="tight")
    print(f"Saved: {out}")

plt.close()
print("Done")
