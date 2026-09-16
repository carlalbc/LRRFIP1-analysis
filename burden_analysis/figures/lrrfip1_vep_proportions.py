"""
LRRFIP1 VEP proportional composition — Panel B only (Sleeve Bite palette).
Coding+regulatory variants excluding intronic, stacked by consequence severity.
JCI specs: 600 dpi, 8pt Arial/Helvetica, TIFF preferred.
"""
import os, numpy as np, matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
matplotlib.rcParams["font.size"] = 8

RESULTS_DIR = os.path.dirname(os.path.abspath(__file__)) if "__file__" in dir() else "."

# --- TWEAKABLE PARAMS ---
FIG_W, FIG_H = 9.0 / 2.54, 12.0 / 2.54   # figure size (cm -> inches)
FS = 8                                       # main font size
FS_SM = 7                                    # small font size
BAR_W = 0.55                                 # bar width
PCT_THRESH = 3                               # min % to show label
# -------------------------

# Data: [LoF, Missense, Synonymous, Splice, UTR, In-frame indel]
cu_focus_orig = [1, 1, 1, 0, 5, 0]
sh_focus_orig = [0, 1, 3, 1, 2, 0]
ko_focus_orig = [47, 230, 118, 106, 458, 5]

# Reorder for stacking (bottom -> top): LoF, Splice, Missense, In-frame indel, UTR, Synonymous
stack_order = [0, 3, 1, 5, 4, 2]
cu_focus = [cu_focus_orig[i] for i in stack_order]
sh_focus = [sh_focus_orig[i] for i in stack_order]
ko_focus = [ko_focus_orig[i] for i in stack_order]

cu_total = sum(cu_focus)
sh_total = sum(sh_focus)
ko_total = sum(ko_focus)

cu_pct = [v / cu_total * 100 if cu_total > 0 else 0 for v in cu_focus]
sh_pct = [v / sh_total * 100 if sh_total > 0 else 0 for v in sh_focus]
ko_pct = [v / ko_total * 100 if ko_total > 0 else 0 for v in ko_focus]

# Sleeve Bite palette (color-hex.com): #175282, #436892, #7582a4, #a79cb7, #d3b3c7
# + 1 interpolated (#8e8fae) for 6 categories. Darkest = most severe at bottom.
colors_stack = [
    "#175282",   # LoF — deepest blue
    "#436892",   # Splice region — dark blue
    "#7582a4",   # Missense — medium blue-purple
    "#8e8fae",   # In-frame indel — interpolated
    "#a79cb7",   # UTR — light mauve
    "#d3b3c7",   # Synonymous — lightest pink
]

labels_stack = ["LoF", "Splice region", "Missense", "In-frame indel", "UTR", "Synonymous"]

groups = [cu_pct, sh_pct, ko_pct]
group_labels = [
    f"Case-only\n(n={cu_total})",
    f"Shared\n(n={sh_total})",
    f"Control-only\n(n={ko_total})",
]

# Plot
fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
x2 = np.arange(3)

for i, (label, color) in enumerate(zip(labels_stack, colors_stack)):
    vals = [g[i] for g in groups]
    ax.bar(x2, vals, BAR_W, bottom=[sum(groups[j][:i]) for j in range(3)],
           color=color, edgecolor="white", linewidth=0.5, label=label)

ax.set_xticks(x2)
ax.set_xticklabels(group_labels, fontsize=FS)
ax.set_ylabel("Proportion of coding + regulatory variants (%)", fontsize=FS)
ax.set_ylim(0, 105)
ax.tick_params(axis="x", pad=2, labelsize=FS)
ax.tick_params(axis="y", pad=2, labelsize=FS_SM)

# Legend: inverted so top-to-bottom matches visual stack
handles, lbls = ax.get_legend_handles_labels()
ax.legend(handles[::-1], lbls[::-1], fontsize=FS_SM, loc="upper right",
          bbox_to_anchor=(1.0, 1.0), frameon=True, framealpha=0.9,
          edgecolor="#ccc", borderpad=0.3, handlelength=1.0, handleheight=0.7)

# Percentage labels — show if >= PCT_THRESH
for j in range(3):
    bottom = 0
    for i in range(len(labels_stack)):
        val = groups[j][i]
        if val >= PCT_THRESH:
            txt_color = "white" if i <= 2 else "#333"
            ax.text(x2[j], bottom + val / 2, f"{val:.0f}%", ha="center", va="center",
                    fontsize=FS, fontweight="bold", color=txt_color)
        bottom += val

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()

for ext in ("png", "pdf", "svg", "tiff"):
    out = os.path.join(RESULTS_DIR, f"lrrfip1_vep_proportions.{ext}")
    if ext == "tiff":
        fig.savefig(out, dpi=600, bbox_inches="tight",
                    pil_kwargs={"compression": "tiff_lzw"})
    else:
        fig.savefig(out, dpi=600, bbox_inches="tight")
    print(f"Saved: {out}")
plt.close()
print("Done")
