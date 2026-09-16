"""
LRRFIP1 burden signal — lollipop plot on gene schematic.
Shows all 14 variants in the CDS+UTR mask (OR=2.66, p=0.015).
Filled circles = case-only, open circles = shared with controls.
Nature-style, compressed introns, light background.

"""
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import seaborn as sns

sns.set_theme(style="white", font_scale=1.1)

RESULTS_DIR = os.path.dirname(os.path.abspath(__file__)) if "__file__" in dir() else "."

# ── Variant data (all 14 in M2_coding mask) ──────────────────────────────────
variants = [
    (237_627_735, "Missense",   "p.Arg31Trp",        "Exon 1",  "Shared"),
    (237_739_252, "Synonymous", "p.His634=",          "Exon 10", "Shared"),
    (237_751_259, "Synonymous", "p.Leu913=",          "Exon 14", "Shared"),
    (237_756_156, "LoF",        "p.Ala368ProfsTer2",  "Exon 16", "Case-only"),
    (237_756_183, "Missense",   "p.Leu376Pro",        "Exon 16", "Case-only"),
    (237_758_797, "Synonymous", "p.Val1351=",         "Exon 20", "Shared"),
    (237_769_948, "Synonymous", "",                    "Exon 24", "Case-only"),
    (237_779_668, "UTR",        "",                    "3'UTR",   "Case-only"),
    (237_780_461, "UTR",        "",                    "3'UTR",   "Case-only"),
    (237_780_597, "UTR",        "",                    "3'UTR",   "Shared"),
    (237_780_871, "UTR",        "",                    "3'UTR",   "Case-only"),
    (237_781_023, "UTR",        "",                    "3'UTR",   "Case-only"),
    (237_781_081, "UTR",        "",                    "3'UTR",   "Case-only"),
    (237_781_086, "UTR",        "",                    "3'UTR",   "Shared"),
]

CONS_COLORS = {
    "LoF":        "#011f4b",
    "Missense":   "#03396c",
    "Synonymous": "#6497b1",
    "UTR":        "#4c516d",
}

# ── Gene structure (key exons only) ──────────────────────────────────────────
exon_regions = [
    (237_627_600, 237_627_900, "1"),
    (237_738_900, 237_739_400, "10"),
    (237_751_000, 237_751_500, "14"),
    (237_755_900, 237_756_300, "16"),
    (237_758_500, 237_759_000, "20"),
    (237_769_700, 237_770_200, "24"),
]
utr_region = (237_779_400, 237_781_200, "3'UTR")

EXON_WIDTH = 1.2
INTRON_GAP = 0.6
UTR_WIDTH = 3.5

regions = []
x_cursor = 0.0
all_blocks = exon_regions + [utr_region]

for i, block in enumerate(all_blocks):
    gstart, gend, label = block
    w = UTR_WIDTH if label == "3'UTR" else EXON_WIDTH
    regions.append((gstart, gend, x_cursor, x_cursor + w, label))
    x_cursor += w
    if i < len(all_blocks) - 1:
        x_cursor += INTRON_GAP

total_width = x_cursor

def genomic_to_display(pos):
    for gstart, gend, dstart, dend, label in regions:
        if gstart <= pos <= gend:
            frac = (pos - gstart) / (gend - gstart) if gend > gstart else 0.5
            return dstart + frac * (dend - dstart)
    best_dist = 1e12
    best_x = 0
    for gstart, gend, dstart, dend, label in regions:
        dist = min(abs(pos - gstart), abs(pos - gend))
        if dist < best_dist:
            best_dist = dist
            best_x = (dstart + dend) / 2
    return best_x

# ── Figure ────────────────────────────────────────────────────────────────────
BG_COLOR = "#f7f9fc"
fig, ax = plt.subplots(figsize=(15, 7.5))
fig.patch.set_facecolor(BG_COLOR)
ax.set_facecolor(BG_COLOR)

GENE_Y = 0.0
STEM_BASE = 0.45
gene_bar_h = 0.28
intron_y = GENE_Y

# Intron line — plain
ax.plot([regions[0][2] - 0.15, regions[-1][3] + 0.15], [intron_y, intron_y],
        color="#b3cde0", lw=3, solid_capstyle="round", zorder=1)

# Exon blocks
for gstart, gend, dstart, dend, label in regions:
    if label == "3'UTR":
        rect = FancyBboxPatch((dstart, intron_y - gene_bar_h * 0.35),
                              dend - dstart, gene_bar_h * 0.7,
                              boxstyle="round,pad=0.04",
                              facecolor="#d4e6f1", edgecolor="#6497b1", lw=1.5, zorder=2)
        ax.add_patch(rect)
        ax.text((dstart + dend) / 2, intron_y - gene_bar_h * 1.5,
                "3' UTR", ha="center", va="top", fontsize=12, color="#4c516d",
                fontweight="bold")
    else:
        rect = FancyBboxPatch((dstart, intron_y - gene_bar_h / 2),
                              dend - dstart, gene_bar_h,
                              boxstyle="round,pad=0.04",
                              facecolor="#005b96", edgecolor="#03396c", lw=1.5, zorder=2)
        ax.add_patch(rect)
        ax.text((dstart + dend) / 2, intron_y - gene_bar_h * 1.5,
                f"Ex. {label}", ha="center", va="top", fontsize=12, color="#444")

# ── Positions and heights ────────────────────────────────────────────────────
x_positions = [genomic_to_display(v[0]) for v in variants]

CDS_H = 2.4
EX16_H = 3.8
UTR_LEVELS = [2.4, 3.1, 3.8]

heights = []
utr_count = 0
for i, (pos, cons, hgvsp, exon, status) in enumerate(variants):
    if exon == "Exon 16" and cons in ("LoF", "Missense"):
        heights.append(EX16_H)
    elif cons == "UTR":
        heights.append(UTR_LEVELS[utr_count % 3])
        utr_count += 1
    else:
        heights.append(CDS_H)

# ── Draw lollipops ────────────────────────────────────────────────────────────
for i, (pos, cons, hgvsp, exon, status) in enumerate(variants):
    x = x_positions[i]
    y_top = heights[i]
    color = CONS_COLORS[cons]
    is_case = status == "Case-only"

    ax.plot([x, x], [STEM_BASE, y_top], color=color, lw=2.0, alpha=0.55, zorder=3,
            solid_capstyle="round")

    if is_case:
        ax.scatter(x, y_top, s=160, color=color, edgecolors=color,
                   linewidths=2.0, zorder=5)
    else:
        ax.scatter(x, y_top, s=160, facecolors=BG_COLOR, edgecolors=color,
                   linewidths=2.0, zorder=5)

# ── Annotations — text only, no arrows ────────────────────────────────────────
lof_idx = 3
ax.text(x_positions[lof_idx] - 1.8, EX16_H + 1.3,
        "p.Ala368ProfsTer2\nframeshift · absent from gnomAD",
        fontsize=12, color="#011f4b", fontweight="bold", ha="center", va="bottom")

mis_idx = 4
ax.text(x_positions[mis_idx] + 2.2, EX16_H + 1.1,
        "p.Leu376Pro\nSIFT: deleterious\nPolyPhen: prob. damaging",
        fontsize=12, color="#03396c", fontweight="bold", ha="center", va="bottom")

# Bracket for exon 16 pair
brace_y = EX16_H + 0.45
ax.plot([x_positions[lof_idx], x_positions[mis_idx]], [brace_y, brace_y],
        color="#03396c", lw=2.2, zorder=4, solid_capstyle="round")
tick_h = 0.06
for xi in [x_positions[lof_idx], x_positions[mis_idx]]:
    ax.plot([xi, xi], [brace_y - tick_h, brace_y + tick_h], color="#03396c", lw=2.2, zorder=4)
ax.text((x_positions[lof_idx] + x_positions[mis_idx]) / 2, brace_y + 0.18,
        "27 bp apart · PF09738 coiled-coil domain", ha="center", va="bottom",
        fontsize=12, color="#03396c", style="italic")

# UTR cluster bracket
utr_indices = [i for i, v in enumerate(variants) if v[1] == "UTR"]
utr_x_min = min(x_positions[i] for i in utr_indices) - 0.12
utr_x_max = max(x_positions[i] for i in utr_indices) + 0.12
utr_brace_y = max(heights[i] for i in utr_indices) + 0.5
ax.plot([utr_x_min, utr_x_max], [utr_brace_y, utr_brace_y],
        color="#4c516d", lw=1.8, zorder=4, solid_capstyle="round")
for xi in [utr_x_min, utr_x_max]:
    ax.plot([xi, xi], [utr_brace_y - tick_h, utr_brace_y + tick_h], color="#4c516d", lw=1.8, zorder=4)
ax.text((utr_x_min + utr_x_max) / 2, utr_brace_y + 0.15,
        "7 variants in 1.4 kb  (5 case-only · 2 shared)", ha="center", va="bottom",
        fontsize=12, color="#4c516d", style="italic")

# Shared missense exon 1
ax.text(x_positions[0], heights[0] + 0.5, "p.Arg31Trp",
        fontsize=12, color="#03396c", ha="center", va="bottom")

# ── Title ─────────────────────────────────────────────────────────────────────
ax.set_title("LRRFIP1 — Variants in the CDS + UTR burden mask\n"
             "OR = 2.66  (95% CI 1.3–5.3)  ·  p = 0.015  ·  14 variants: 8 case-only, 6 shared",
             fontsize=14, fontweight="bold", pad=12, linespacing=1.6)

# ── Legend ────────────────────────────────────────────────────────────────────
legend_elements = [
    plt.scatter([], [], s=80, color="#333", edgecolors="white", linewidths=0.8,
                label="Case-only (filled)"),
    plt.scatter([], [], s=80, facecolors="white", edgecolors="#333", linewidths=1.8,
                label="Shared with controls (open)"),
    mpatches.Patch(color="#011f4b", label="Loss-of-function"),
    mpatches.Patch(color="#03396c", label="Missense"),
    mpatches.Patch(color="#6497b1", label="Synonymous"),
    mpatches.Patch(color="#4c516d", label="3' UTR"),
]
leg = ax.legend(handles=legend_elements, loc="upper center", fontsize=12,
          frameon=True, framealpha=0.95, edgecolor="#ccc", ncol=3,
          bbox_to_anchor=(0.5, 1.0))
leg.get_frame().set_linewidth(0.8)

# ── Axes cleanup ──────────────────────────────────────────────────────────────
ax.set_xlim(-1.0, total_width + 1.5)
ax.set_ylim(-1.2, EX16_H + 3.0)
ax.set_yticks([])
ax.set_xticks([])
for spine in ax.spines.values():
    spine.set_visible(False)

# Gene label — close to gene bar
ax.text(-0.3, GENE_Y, "LRRFIP1", ha="right", va="center",
        fontsize=14, fontweight="bold", color="#011f4b", style="italic")

# Chromosome coords — black
ax.text((regions[0][2] + regions[-1][3]) / 2, -1.0,
        "chr2: 237,627,735 – 237,781,086  (153 kb)",
        ha="center", fontsize=12, color="black")

plt.tight_layout()

# Uncomment to save:
# for ext in ("png", "pdf", "svg"):
#                 dpi=200, bbox_inches="tight", facecolor=BG_COLOR)

plt.show()
