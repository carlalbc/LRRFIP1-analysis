"""
LRRFIP1-only burden forest plot — presentation panel.
Direction labels (▲/▼) only for significant masks (p<0.05).
"""
import os, numpy as np, matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
matplotlib.rcParams["font.size"] = 10

RESULTS_DIR = os.path.dirname(os.path.abspath(__file__)) if "__file__" in dir() else "."

# --- TWEAKABLE PARAMS ---
FIG_W, FIG_H = 10, 3.8          # figure size in inches
CI_CLIP = 10.0                    # max x-axis (CIs beyond this get arrow)
FS_LABEL = 10                     # y-axis mask label font size
FS_OR_SIG = 11                    # significant OR annotation font size
FS_OR_NS = 10                     # non-significant OR annotation font size
FS_DIR = 11                       # PROTECTIVE / RISK-INCREASING font size
DOT_SIZE = 40                     # point estimate marker size
LINE_W = 1.5                      # CI line width
# -------------------------

MASK_INFO = {
    "M1_cds": "CDS",
    "M2_coding": "CDS + UTR",
    "M3_regulatory": "UTR + promoter",
    "M4_intronic": "Intronic",
    "M5_genic": "CDS + UTR + intronic",
    "M6_all": "All regions",
}
MASK_ORDER = list(MASK_INFO.keys())
SKIP = {"LRRFIP1_whole", "RBM44_whole"}

# Y-axis labels (en-dash style for presentation)
MASK_LABEL = {
    "M1_cds": "LRRFIP1 – CDS",
    "M2_coding": "LRRFIP1 – CDS + UTR",
    "M3_regulatory": "LRRFIP1 – UTR + promoter",
    "M4_intronic": "LRRFIP1 – intronic",
    "M5_genic": "LRRFIP1 – exonic + UTR + intronic",
    "M6_all": "LRRFIP1 – all regions",
}

def parse(fp, tier):
    rows = []
    with open(fp) as f:
        for line in f:
            if line.startswith("#") or line.startswith("CHROM"):
                continue
            p = line.strip().split()
            if len(p) < 12 or p[7] != "ADD":
                continue
            rid = p[2]
            gene = mask = rest = None
            for m in MASK_INFO:
                i = rid.find("." + m)
                if i >= 0:
                    gene, mask, rest = rid[:i], m, rid[i + 1:]
                    break
            if gene is None or gene in SKIP:
                continue
            bn = rest[len(mask):]
            bn = bn[1:] if bn.startswith(".") else bn
            try:
                b, s, l = float(p[8]), float(p[9]), float(p[11])
            except:
                continue
            rows.append({
                "tier": tier, "gene": gene, "mask": mask, "bin": bn,
                "beta": b, "OR": np.exp(b),
                "CI_lo": np.exp(b - 1.96 * s),
                "CI_hi": np.exp(b + 1.96 * s),
                "pvalue": 10 ** (-l),
            })
    return rows

# Load data — LRRFIP1 only, primary tier, singleton bin
data = []
fp = os.path.join(RESULTS_DIR, "scoliosis_skato_v3_primary_scoliosis.regenie")
if os.path.isfile(fp):
    data.extend(parse(fp, "primary"))

fdata = [r for r in data
         if r["bin"] == "singleton" and r["tier"] == "primary"
         and r["gene"] == "LRRFIP1"]

rows = []
for m in MASK_ORDER:
    hit = [r for r in fdata if r["mask"] == m]
    if hit:
        rows.append(hit[0])
print(f"LRRFIP1 rows: {len(rows)}")

# Plot
fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
y_pos = list(range(len(rows)))

for i, r in enumerate(rows):
    yi = y_pos[i]
    orv = r["OR"]
    clo = max(r["CI_lo"], 1e-3)
    chi = r["CI_hi"]
    unr = chi > CI_CLIP
    chip = min(chi, CI_CLIP)
    col = "#c0392b" if r["beta"] > 0 else "#1a6faf"

    # CI line
    ax.plot([clo, chip], [yi, yi], color=col, lw=LINE_W,
            alpha=0.35 if unr else 0.85, zorder=4,
            ls="--" if unr else "-")
    if unr:
        ax.annotate("", xy=(CI_CLIP, yi), xytext=(CI_CLIP - 0.15, yi),
                    arrowprops=dict(arrowstyle="->", color=col, lw=1.0), zorder=4)

    # Point estimate
    ax.scatter(min(orv, CI_CLIP), yi, color=col, s=DOT_SIZE, marker="o",
               zorder=5, edgecolors="white", linewidths=0.5)

    # OR annotation — direction label only for significant
    sig = r["pvalue"] < 0.05
    if sig:
        ds = "▲" if r["beta"] > 0 else "▼"
        pval = r["pvalue"]
        ps = f"p = {pval:.3f}" if pval >= 0.001 else f"p = {pval:.1e}"
        txt = f"OR = {orv:.2f}  {ds}  {ps}"
        tcol = "#1a6faf" if r["beta"] < 0 else "#cc0000"
        fsz = FS_OR_SIG
    else:
        if unr:
            txt = f"OR = {orv:.2f}  ▲  (wide CI)"
        else:
            txt = f"OR = {orv:.2f}"
        tcol = "#888"
        fsz = FS_OR_NS

    tw = "bold" if sig else "normal"
    ax.text(CI_CLIP + 0.3, yi, txt, va="center", ha="left", fontsize=fsz,
            color=tcol, fontweight=tw, clip_on=False, zorder=5)

# Null line
ax.axvline(1.0, color="#999", lw=0.8, ls="--", zorder=1)

# Shaded direction regions
ax.axvspan(1.0, CI_CLIP + 0.3, alpha=0.04, color="#c0392b", zorder=0)
ax.axvspan(0, 1.0, alpha=0.04, color="#1a6faf", zorder=0)

# Y-axis labels
ylbl = [MASK_LABEL[r["mask"]] for r in rows]
ax.set_yticks(y_pos)
ax.set_yticklabels(ylbl, fontsize=FS_LABEL)
ax.invert_yaxis()
ax.set_xlim(0, CI_CLIP + 0.3)
ax.set_ylim(len(rows) - 0.5, -1.5)
ax.tick_params(axis="x", labelsize=9)
ax.tick_params(axis="y", pad=4)

# Direction labels at top
ax.text(CI_CLIP * 0.75, -0.9, "RISK-INCREASING →", color="#c0392b",
        ha="center", fontsize=FS_DIR, fontweight="bold", style="italic")
ax.text(0.5, -0.9, "← PROTECTIVE", color="#1a6faf",
        ha="center", fontsize=FS_DIR, fontweight="bold", style="italic")

# Legend
r_marker = Line2D([0], [0], marker="^", color="w", markerfacecolor="#c0392b",
                   markersize=10, label="Risk-increasing (OR > 1)")
b_marker = Line2D([0], [0], marker="v", color="w", markerfacecolor="#1a6faf",
                   markersize=10, label="Protective (OR < 1)")
ax.legend(handles=[r_marker, b_marker], loc="lower center", ncol=2, fontsize=9,
          frameon=True, framealpha=0.9, edgecolor="#ccc",
          bbox_to_anchor=(0.45, -0.25), borderpad=0.4)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()

for ext in ("png", "pdf", "svg"):
    out = os.path.join(RESULTS_DIR, f"lrrfip1_burden_forest.{ext}")
    fig.savefig(out, dpi=600, bbox_inches="tight")
    print(f"Saved: {out}")
plt.close()
print("Done")
