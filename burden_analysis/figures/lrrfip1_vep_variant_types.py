"""
LRRFIP1 VEP variant types by carrier status — Sleeve Bite palette.
SNV / Deletion / Insertion proportional stacked bars.
JCI specs: 600 dpi, 8pt Arial/Helvetica, TIFF preferred.
"""
import os, numpy as np, matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
matplotlib.rcParams["font.size"] = 8

RESULTS_DIR = os.path.dirname(os.path.abspath(__file__)) if "__file__" in dir() else "."

# --- TWEAKABLE PARAMS ---
FIG_W, FIG_H = 9.0 / 2.54, 8.0 / 2.54   # figure size (cm -> inches)
FS = 8                                     # main font size
FS_SM = 7                                  # small font size
BAR_W = 0.55                               # bar width
PCT_THRESH = 3                             # min % to show label
# -------------------------

# Data: Case-only, Shared, Control-only
cu_counts = [131, 8, 2]     # SNV, Deletion, Insertion
sh_counts = [3615, 260, 117]
ko_counts = [39493, 3060, 1207]

cu_total = sum(cu_counts)
sh_total = sum(sh_counts)
ko_total = sum(ko_counts)

cu_pct = [v / cu_total * 100 for v in cu_counts]
sh_pct = [v / sh_total * 100 for v in sh_counts]
ko_pct = [v / ko_total * 100 for v in ko_counts]

# Sleeve Bite palette (3 colors): start from dark blue, not deepest
colors = ["#436892", "#7582a4", "#d3b3c7"]
labels = ["SNV", "Deletion", "Insertion"]

groups = [cu_pct, sh_pct, ko_pct]
group_labels = [
    f"Case-only\n(n={cu_total:,})",
    f"Shared\n(n={sh_total:,})",
    f"Control-only\n(n={ko_total:,})",
]

# Plot
fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
x = np.arange(3)

for i, (label, color) in enumerate(zip(labels, colors)):
    vals = [g[i] for g in groups]
    ax.bar(x, vals, BAR_W, bottom=[sum(groups[j][:i]) for j in range(3)],
           color=color, edgecolor="white", linewidth=0.5, label=label)

ax.set_xticks(x)
ax.set_xticklabels(group_labels, fontsize=FS)
ax.set_ylabel("Proportion of variants (%)", fontsize=FS)
ax.set_ylim(0, 105)
ax.tick_params(axis="x", pad=2, labelsize=FS)
ax.tick_params(axis="y", pad=2, labelsize=FS_SM)

# Legend
handles, lbls = ax.get_legend_handles_labels()
ax.legend(handles[::-1], lbls[::-1], fontsize=FS_SM, loc="upper right",
          bbox_to_anchor=(1.0, 1.0), frameon=True, framealpha=0.9,
          edgecolor="#ccc", borderpad=0.3, handlelength=1.0, handleheight=0.7)

# Percentage labels
for j in range(3):
    bottom = 0
    for i in range(len(labels)):
        val = groups[j][i]
        if val >= PCT_THRESH:
            txt_color = "white" if i == 0 else "#333"
            ax.text(x[j], bottom + val / 2, f"{val:.0f}%", ha="center", va="center",
                    fontsize=FS, fontweight="bold", color=txt_color)
        bottom += val

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()

for ext in ("png", "pdf", "svg", "tiff"):
    out = os.path.join(RESULTS_DIR, f"lrrfip1_vep_variant_types.{ext}")
    if ext == "tiff":
        fig.savefig(out, dpi=600, bbox_inches="tight",
                    pil_kwargs={"compression": "tiff_lzw"})
    else:
        fig.savefig(out, dpi=600, bbox_inches="tight")
    print(f"Saved: {out}")
plt.close()
print("Done")
