"""
Forest plot: Gene names as left-side labels with vertical grouping lines.
Mask names only on y-axis (no gene prefix). Clean layout with correct
mask terminology (CDS, CDS + UTR, etc.).
"""
import os, numpy as np, matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
matplotlib.rcParams["font.size"] = 8

RESULTS_DIR = os.path.dirname(os.path.abspath(__file__)) if "__file__" in dir() else "."
FIG_W, FIG_H = 23.0 / 2.54, 17.25 / 2.54

MASK_INFO = {
    "M1_cds": "Coding",
    "M2_coding": "Coding + UTR",
    "M3_regulatory": "UTR+promoter",
    "M4_intronic": "Intronic",
    "M5_genic": "Coding+UTR+intronic",
    "M6_all": "All regions",
}
MASK_ORDER = list(MASK_INFO.keys())
GENE_ORDER = ["COL6A3", "MLPH", "PRLH", "RAB17",
              "LRRFIP1", "LRRFIP1_RBM44_overlap", "RBM44"]
SKIP = {"LRRFIP1_whole", "RBM44_whole"}
GENE_DISPLAY = {
    "COL6A3": "COL6A3", "MLPH": "MLPH", "PRLH": "PRLH", "RAB17": "RAB17",
    "LRRFIP1": "LRRFIP1",
    "LRRFIP1_RBM44_overlap": "LRRFIP1/RBM44\noverlap",
    "RBM44": "RBM44",
}
CI_CLIP = 10.0
BONF = 0.05 / 42
FS = 8; FS_SM = 7

# Load from pre-computed TSV (directionality_burden.tsv)
TSV = os.path.join(RESULTS_DIR, "directionality_burden.tsv")
import csv
data = []
with open(TSV) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        data.append({
            "tier": row["tier"],
            "gene": row["gene"],
            "mask": row["mask"],
            "bin": row["bin"],
            "beta": float(row["beta"]),
            "OR": float(row["OR"]),
            "CI_lo": float(row["CI_lo"]),
            "CI_hi": float(row["CI_hi"]),
            "pvalue": float(row["pvalue"]),
        })

fdata = [r for r in data
         if r["bin"] == "singleton" and r["tier"] == "primary"
         and r["gene"] in GENE_ORDER]

rows = []
for g in GENE_ORDER:
    for m in MASK_ORDER:
        hit = [r for r in fdata if r["gene"] == g and r["mask"] == m]
        if hit:
            rows.append(hit[0])
print(f"Loaded {len(rows)} rows ({len(set(r['gene'] for r in rows))} genes x {len(MASK_ORDER)} masks)")

# Build y-positions with gene-group gaps
y_pos = []
gene_spans = {}
GAP = 0.8
y = 0
prev_gene = None
for i, r in enumerate(rows):
    if prev_gene is not None and r["gene"] != prev_gene:
        y += GAP
    y_pos.append(y)
    gene_spans.setdefault(r["gene"], [y, y])
    gene_spans[r["gene"]][1] = y
    prev_gene = r["gene"]
    y += 1
y_max = y_pos[-1]

# Plot — wider left margin for gene labels
fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
fig.subplots_adjust(left=0.36)

# Alternating gene-group backgrounds
gseen = []
for r in rows:
    if r["gene"] not in gseen:
        gseen.append(r["gene"])
gbg = ["#f0f4fa", "#ffffff"]
for gi, g in enumerate(gseen):
    ymin, ymax_g = gene_spans[g]
    ax.axhspan(ymin - 0.5, ymax_g + 0.5, color=gbg[gi % 2], zorder=0)

# CI lines + point estimates + OR annotations
for i, r in enumerate(rows):
    yi = y_pos[i]
    orv = r["OR"]
    clo = max(r["CI_lo"], 1e-3)
    chi = r["CI_hi"]
    unr = chi > CI_CLIP
    chip = min(chi, CI_CLIP)
    col = "#c0392b" if r["beta"] > 0 else "#1a6faf"

    ax.plot([clo, chip], [yi, yi], color=col, lw=0.8,
            alpha=0.35 if unr else 0.8, zorder=4,
            ls="--" if unr else "-")
    if unr:
        ax.annotate("", xy=(CI_CLIP, yi), xytext=(CI_CLIP - 0.2, yi),
                    arrowprops=dict(arrowstyle="->", color=col, lw=0.5), zorder=4)

    mk = "D" if r["pvalue"] < BONF else "o"
    sz = 16 if r["pvalue"] < BONF else 8
    ax.scatter(min(orv, CI_CLIP), yi, color=col, s=sz, marker=mk, zorder=5,
               edgecolors="white" if r["pvalue"] < BONF else "none",
               linewidths=0.3, alpha=0.4 if unr else 1.0)

    # OR annotation
    sig = r["pvalue"] < 0.05
    if sig:
        ds = "▲ risk" if r["beta"] > 0 else "▼ prot."
        ps = f"  p={r['pvalue']:.2e}"
        txt = f"OR={orv:.2f}  {ds}{ps}"
        tcol = "#1a6faf" if r["beta"] < 0 else "#cc0000"
    else:
        txt = f"OR={orv:.2f}"
        tcol = "#888"

    tw = "bold" if sig else "normal"
    ax.text(1.02, yi, txt, va="center", ha="left", fontsize=FS_SM,
            color=tcol, fontweight=tw, clip_on=False, zorder=5,
            transform=ax.get_yaxis_transform())

# Null line
ax.axvline(1.0, color="#999", lw=0.7, ls="--", zorder=1)

# Y-axis: mask names only (no gene prefix)
ylbl = [MASK_INFO[r["mask"]] for r in rows]
ax.set_yticks(y_pos)
ax.set_yticklabels(ylbl, fontsize=FS_SM)
ax.invert_yaxis()
ax.set_xlim(0, CI_CLIP + 0.5)
ax.set_ylim(y_max + 0.5, -3.2)
ax.tick_params(axis="x", labelsize=FS_SM)
ax.tick_params(axis="y", pad=10, length=0)

# Gene labels with vertical bracket lines
for g in gseen:
    ymin_g, ymax_g = gene_spans[g]
    ymid = (ymin_g + ymax_g) / 2.0

    # Gene name — placed well left of y-axis
    label = GENE_DISPLAY.get(g, g)
    label_fs = 7.5 if g == "LRRFIP1_RBM44_overlap" else FS
    ax.text(-0.36, ymid, label, va="center", ha="center", fontsize=label_fs,
            fontweight="bold", transform=ax.get_yaxis_transform(),
            clip_on=False)

    # Vertical bracket line
    bracket_x = -0.24
    ax.plot([bracket_x, bracket_x],
            [ymin_g - 0.3, ymax_g + 0.3],
            transform=ax.get_yaxis_transform(),
            color="#999999", lw=0.8, solid_capstyle="round",
            clip_on=False, zorder=10)
    # Small horizontal ticks pointing RIGHT (toward labels), short
    tick_len = 0.006
    ax.plot([bracket_x, bracket_x + tick_len],
            [ymin_g - 0.3, ymin_g - 0.3],
            transform=ax.get_yaxis_transform(),
            color="#999999", lw=0.8, clip_on=False, zorder=10)
    ax.plot([bracket_x, bracket_x + tick_len],
            [ymax_g + 0.3, ymax_g + 0.3],
            transform=ax.get_yaxis_transform(),
            color="#999999", lw=0.8, clip_on=False, zorder=10)

ax.set_xlabel("Odds Ratio (ADD burden, MAC ≤ 5 singleton)  [CI capped at 10; → = wider]", fontsize=FS)

# Shaded direction regions
xmn, xmx = ax.get_xlim()
ax.axvspan(1.0, xmx, alpha=0.03, color="#c0392b", zorder=0)
ax.axvspan(xmn, 1.0, alpha=0.03, color="#1a6faf", zorder=0)

# Direction labels
ax.text(xmx * 0.99, -1, "RISK-INCREASING →", color="#c0392b",
        ha="right", fontsize=9.5, fontweight="bold", style="italic")
ax.text(xmn + 0.01, -1, "← PROTECTIVE", color="#1a6faf",
        ha="left", fontsize=9.5, fontweight="bold", style="italic")

# Legend
rp = mpatches.Patch(color="#c0392b", label="OR > 1 (risk-increasing)")
bp = mpatches.Patch(color="#1a6faf", label="OR < 1 (protective)")
ax.legend(handles=[rp, bp], loc="upper center", ncol=2, fontsize=FS_SM,
          frameon=True, framealpha=0.90, edgecolor="#ccc",
          bbox_to_anchor=(0.50, 1.038), borderpad=0.2,
          handlelength=1, handleheight=0.1)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()
fig.subplots_adjust(left=0.36)

for ext in ("png", "pdf", "svg", "tiff"):
    out = os.path.join(RESULTS_DIR, f"directionality_burden.{ext}")
    if ext == "tiff":
        fig.savefig(out, dpi=600, bbox_inches="tight",
                    pil_kwargs={"compression": "tiff_lzw"})
    else:
        fig.savefig(out, dpi=600, bbox_inches="tight")
    print(f"Saved: {out}")
plt.close()
print("Done")
