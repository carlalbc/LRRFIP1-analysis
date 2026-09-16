"""
LRRFIP1: Unbiased IsoSeq isoform distribution analysis.
Maps all 102 PB isoforms to Gencode v49 reference transcripts and shows
per-sample distribution — the long-read equivalent of the Illumina figure.

Key evidence:
1. Isoform 1 (MANE Select, ENST00000308482.14): ZERO IsoSeq matches
2. Isoform 2 (APPRIS P1, ENST00000863709.1): dominant transcript
3. Isoform 6 (retained intron): unexpectedly strong support

Reads directly from pigeon classification + isoform summary files.
Magma palette. JCI specs: 600 dpi, 8pt Arial/Helvetica.
"""
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import warnings

matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
matplotlib.rcParams["font.size"] = 8

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data") if "__file__" in dir() else "./data"
RESULTS_DIR = os.path.dirname(os.path.abspath(__file__)) if "__file__" in dir() else "."
os.makedirs(RESULTS_DIR, exist_ok=True)

# --- TWEAKABLE PARAMS ---
FIG_W, FIG_H = 22.0 / 2.54, 18.0 / 2.54
FS = 8
FS_SM = 7
FS_LG = 9
FS_XS = 6

# --- MAGMA PALETTE ---
MAG1 = "#2d1160"   # deep indigo
MAG2 = "#6e1d81"   # purple
MAG3 = "#b63679"   # magenta
MAG4 = "#f1605d"   # coral
MAG5 = "#feaf77"   # peach

# ============================================================
# 1. READ v49 CLASSIFICATION → PB-to-transcript mapping
# ============================================================
v49_class_file = os.path.join(
    BASE, "pigeon_v49_rerun", "outputs",
    "collapse_isoforms_classification.filtered_lite_classification.txt")

pb_to_v49 = {}  # PB_ID → (category, transcript)
with open(v49_class_file) as f:
    header = f.readline()  # skip header
    for line in f:
        cols = line.strip().split("\t")
        pb_id = cols[0]
        if not pb_id.startswith("PB.4795."):
            continue
        cat = cols[5]       # structural_category
        tx = cols[7]        # associated_transcript
        pb_to_v49[pb_id] = (cat, tx)

print(f"\n  {'─'*56}")
print(f"  │  LRRFIP1 IsoSeq Isoform Distribution — Magma palette")
print(f"  {'─'*56}\n")
print(f"  ▸ v49 classification loaded: {len(pb_to_v49)} LRRFIP1 isoforms")

# ============================================================
# 2. READ ISOFORM SUMMARY → per-sample FL counts
# ============================================================
summary_file = os.path.join(BASE, "LRRFIP1_isoform_summary.tsv")

# Columns: isoform, length, exons, cat, subcat, tx, FL_BS1..BS4, FL_total, ...
# BS1=III-13 healthy son, BS2=II-8.5 healthy mother,
# BS3=III-12 affected son, BS4=II-8 affected father
sample_names = ["III-13\nhealthy", "II-8.5\nhealthy",
                "III-12\naffected", "II-8\naffected"]
sample_labels = ["Son", "Mother", "Son", "Father"]
sample_status = ["WT", "WT", "c.96+34737C>T", "c.96+34737C>T"]

pb_fl = {}  # PB_ID → [BS1, BS2, BS3, BS4]
with open(summary_file) as f:
    for line in f:
        cols = line.strip().split("\t")
        if cols[0] == "isoform" or not cols[0].startswith("PB.4795."):
            continue
        pb_id = cols[0]
        fl = [int(cols[6]), int(cols[7]), int(cols[8]), int(cols[9])]
        pb_fl[pb_id] = fl

print(f"  ▸ Isoform summary loaded: {len(pb_fl)} isoforms with FL counts")

# ============================================================
# 3. DEFINE ISOFORM GROUPS
# ============================================================
# The 6 selected reference isoforms
ref_isoforms = {
    "ENST00000308482.14": ("Isoform 1\n(MANE Select)", "#175282"),
    "ENST00000863709.1":  ("Isoform 2\n(APPRIS P1)", "#2a6a9e"),
    "ENST00000392000.4":  ("Isoform 3", "#436892"),
    "ENST00000244815.9":  ("Isoform 4", "#7582a4"),
    "ENST00000289175.10": ("Isoform 5", "#a79cb7"),
    "ENST00000478958.5":  ("Isoform 6\n(ret. intron)", "#d3b3c7"),
}

# Group all PB isoforms
groups = {
    "Iso2 (APPRIS P1)":       {"color": MAG1, "fl": [0, 0, 0, 0], "count": 0,
                                 "enst": "ENST00000863709.1"},
    "Iso5":                    {"color": MAG4, "fl": [0, 0, 0, 0], "count": 0,
                                 "enst": "ENST00000289175.10"},
    "Iso6 (ret. intron)":      {"color": MAG5, "fl": [0, 0, 0, 0], "count": 0,
                                 "enst": "ENST00000478958.5"},
    "Other annotated":         {"color": MAG3, "fl": [0, 0, 0, 0], "count": 0,
                                 "enst": None},
    "Novel (NIC/NNC)":         {"color": MAG2, "fl": [0, 0, 0, 0], "count": 0,
                                 "enst": None},
}

# Classify each PB isoform
other_enst_set = set()
for pb_id, fl in pb_fl.items():
    if pb_id not in pb_to_v49:
        continue
    cat, tx = pb_to_v49[pb_id]

    if tx == "ENST00000863709.1":
        grp = "Iso2 (APPRIS P1)"
    elif tx == "ENST00000289175.10":
        grp = "Iso5"
    elif tx == "ENST00000478958.5":
        grp = "Iso6 (ret. intron)"
    elif tx == "novel":
        grp = "Novel (NIC/NNC)"
    else:
        grp = "Other annotated"
        other_enst_set.add(tx)

    groups[grp]["fl"] = [groups[grp]["fl"][i] + fl[i] for i in range(4)]
    groups[grp]["count"] += 1

print(f"  ▸ Other annotated refs: {len(other_enst_set)} unique transcripts")

# Print group summary
total_fl = [0, 0, 0, 0]
for gdata in groups.values():
    total_fl = [total_fl[i] + gdata["fl"][i] for i in range(4)]
grand = sum(total_fl)

print(f"\n  ┌{'─'*75}┐")
print(f"  │  {'Group':<23} {'#PB':>4}  {'III-13':>7} {'II-8.5':>7} {'III-12':>7} {'II-8':>7} {'Total':>6} {'%':>6}  │")
print(f"  │  {'':23} {'':>4}  {'healthy':>7} {'healthy':>7} {'affect':>7} {'affect':>7}               │")
print(f"  ├{'─'*75}┤")
for gname, gdata in groups.items():
    t = sum(gdata["fl"])
    pct = t / grand * 100 if grand > 0 else 0
    print(f"  │  {gname:<23} {gdata['count']:>4}  {gdata['fl'][0]:>7} {gdata['fl'][1]:>7} "
          f"{gdata['fl'][2]:>7} {gdata['fl'][3]:>7} {t:>6} {pct:>5.1f}%  │")
print(f"  ├{'─'*75}┤")
print(f"  │  {'TOTAL':<23} {'':>4}  {total_fl[0]:>7} {total_fl[1]:>7} "
      f"{total_fl[2]:>7} {total_fl[3]:>7} {grand:>6}        │")
print(f"  └{'─'*75}┘")

# ============================================================
# 4. FIGURE
# ============================================================
fig = plt.figure(figsize=(FIG_W, FIG_H))
gs = fig.add_gridspec(2, 2, width_ratios=[1.2, 1], height_ratios=[1.2, 1],
                       hspace=0.45, wspace=0.35)

# ---- PANEL A: IsoSeq isoform fraction per sample (stacked bars) ----
ax1 = fig.add_subplot(gs[0, 0])
ax1.set_title("A  IsoSeq isoform distribution per sample",
              fontsize=FS_LG, fontweight="bold", loc="left", pad=8)

group_order = ["Iso2 (APPRIS P1)", "Novel (NIC/NNC)", "Other annotated",
               "Iso5", "Iso6 (ret. intron)"]

x = np.arange(4)
bar_w = 0.55

# Compute fractions
fractions = {}
for gname in group_order:
    gfl = groups[gname]["fl"]
    fracs = [gfl[i] / total_fl[i] if total_fl[i] > 0 else 0 for i in range(4)]
    fractions[gname] = fracs

bottoms = np.zeros(4)
for gname in group_order:
    fracs = fractions[gname]
    color = groups[gname]["color"]
    ax1.bar(x, fracs, bar_w, bottom=bottoms, color=color,
            edgecolor="white", linewidth=0.5, label=gname)

    # Labels inside bars (if tall enough)
    for i in range(4):
        if fracs[i] >= 0.06:
            pct = f"{fracs[i]*100:.0f}%"
            txt_color = "#333" if gname in ("Iso6 (ret. intron)",) else "white"
            ax1.text(i, bottoms[i] + fracs[i] / 2, pct,
                     ha="center", va="center", fontsize=FS_XS, fontweight="bold",
                     color=txt_color)
    bottoms += np.array(fracs)

ax1.set_xticks(x)
ax1.set_xticklabels(sample_labels, fontsize=FS_SM, ha="center")
ax1.set_ylabel("IsoSeq isoform fraction", fontsize=FS)
ax1.set_ylim(0, 1.05)
ax1.tick_params(axis="y", labelsize=FS_SM)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)

# Legend — figure-level, centered horizontally below titles
handles = [mpatches.Patch(facecolor=groups[g]["color"], edgecolor="white",
                           label=g) for g in reversed(group_order)]
fig.legend(handles=handles, fontsize=FS_XS, loc="upper center",
           bbox_to_anchor=(0.5, 0.97), ncol=len(group_order),
           frameon=True, framealpha=0.9, edgecolor="#ccc",
           borderpad=0.3, handlelength=1.0, handleheight=0.7,
           columnspacing=1.0)

# Total FL counts above bars
for i in range(4):
    ax1.text(i, 1.02, f"n={total_fl[i]}", ha="center", va="bottom",
             fontsize=FS_XS, color="#666")

# ---- PANEL B: Absolute FL counts per group per sample ----
ax2 = fig.add_subplot(gs[0, 1])
ax2.set_title("B  Absolute FL read counts",
              fontsize=FS_LG, fontweight="bold", loc="left", pad=8)

bottoms2 = np.zeros(4)
for gname in group_order:
    gfl = np.array(groups[gname]["fl"], dtype=float)
    color = groups[gname]["color"]
    ax2.bar(x, gfl, bar_w, bottom=bottoms2, color=color,
            edgecolor="white", linewidth=0.5)

    for i in range(4):
        if gfl[i] >= 15:
            txt_color = "#333" if gname in ("Iso6 (ret. intron)",) else "white"
            ax2.text(i, bottoms2[i] + gfl[i] / 2, str(int(gfl[i])),
                     ha="center", va="center", fontsize=FS_XS, fontweight="bold",
                     color=txt_color)
    bottoms2 += gfl

ax2.set_xticks(x)
ax2.set_xticklabels([sample_labels[i] for i in range(4)], fontsize=FS_SM)
ax2.set_ylabel("Full-length reads", fontsize=FS)
ax2.tick_params(axis="y", labelsize=FS_SM)
ax2.spines["top"].set_visible(False)
ax2.spines["right"].set_visible(False)

# ---- PANEL C: Isoform 1 absence evidence ----
ax3 = fig.add_subplot(gs[1, 0])
ax3.set_xlim(0, 10)
ax3.set_ylim(0, 6)
ax3.axis("off")
ax3.set_title("C  Isoform 1 (MANE Select) is absent from IsoSeq",
              fontsize=FS_LG, fontweight="bold", loc="left", pad=8)

evidence = [
    ("PIGEON v49 filtered classification", "0 PB isoforms matched", "#c44"),
    ("PIGEON v49 unfiltered classification", "0 PB isoforms matched", "#c44"),
    ("PIGEON v39 filtered classification", "0 PB isoforms matched", "#c44"),
    ("Salmon full transcriptome (Illumina)", "0 TPM — no unique k-mers", "#c44"),
    ("Salmon restricted to 6 isoforms", "detectable (forced assignment)", "#e8a735"),
]

for i, (source, result, col) in enumerate(evidence):
    y = 5.0 - i * 0.85
    # Bullet with X or checkmark (ASCII-safe for Arial)
    ax3.plot(0.3, y, "o", color=col, markersize=6, zorder=3)
    symbol = "X" if "0" in result else "+"
    ax3.text(0.3, y, symbol, ha="center", va="center", fontsize=5,
             color="white", fontweight="bold", zorder=4,
             fontfamily="sans-serif")
    # Source
    ax3.text(0.7, y + 0.1, source, fontsize=FS_SM, va="center",
             fontweight="bold", color="#333")
    # Result
    ax3.text(0.7, y - 0.2, result, fontsize=FS_XS, va="center",
             color=col, style="italic")

# Conclusion box
ax3.text(0.3, 0.6,
         "ENST00000308482.14 (NM_001137550, 640aa) is not detected by\n"
         "long-read IsoSeq in any sample or annotation version.\n"
         "Illumina short-read detection was an artifact of restricted quantification.",
         fontsize=FS_SM, va="top", color="#333",
         bbox=dict(boxstyle="round,pad=0.4", facecolor="#fff3f3",
                   edgecolor="#c44", alpha=0.8))

# ---- PANEL D: Summary of well-supported isoforms ----
ax4 = fig.add_subplot(gs[1, 1])
ax4.set_xlim(0, 10)
ax4.set_ylim(0, 6)
ax4.axis("off")
ax4.set_title("D  IsoSeq-validated LRRFIP1 isoforms",
              fontsize=FS_LG, fontweight="bold", loc="left", pad=8)

validated = [
    ("Isoform 2", "ENST00000863709.1", "APPRIS P1, 404aa",
     "7 FSM + 1 ISM", "462 FL", MAG1),
    ("Isoform 6", "ENST00000478958.5", "Retained intron",
     "4 FSM + 3 ISM", "141 FL", MAG5),
    ("Isoform 5", "ENST00000289175.10", "752aa",
     "3 FSM", "11 FL", MAG4),
]

not_validated = [
    ("Isoform 1", "ENST00000308482.14", "MANE Select, 640aa", "0 matches"),
    ("Isoform 3", "ENST00000392000.4", "808aa", "3 filtered (IntraPriming)"),
    ("Isoform 4", "ENST00000244815.9", "784aa", "3 filtered (IntraPriming)"),
]

y = 5.2
ax4.text(0.2, y, "Validated (genuine IsoSeq support):", fontsize=FS_SM,
         fontweight="bold", color=MAG1, va="center")
y -= 0.15
for name, enst, desc, matches, fl, color in validated:
    y -= 0.6
    ax4.add_patch(plt.Rectangle((0.2, y - 0.2), 0.15, 0.4,
                                 facecolor=color, edgecolor="none", zorder=2))
    ax4.text(0.5, y + 0.08, f"{name}  ({enst})",
             fontsize=FS_SM, va="center", fontweight="bold")
    ax4.text(0.5, y - 0.18, f"{desc} — {matches} — {fl}",
             fontsize=FS_XS, va="center", color="#555")

y -= 0.7
ax4.text(0.2, y, "Not validated:", fontsize=FS_SM,
         fontweight="bold", color="#c44", va="center")
y -= 0.15
for name, enst, desc, status in not_validated:
    y -= 0.5
    ax4.text(0.5, y + 0.05, f"{name}  ({enst})",
             fontsize=FS_SM, va="center", color="#999")
    ax4.text(0.5, y - 0.18, f"{desc} — {status}",
             fontsize=FS_XS, va="center", color="#999", style="italic")

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    plt.tight_layout()

for ext in ("png", "pdf", "svg", "tiff"):
    out = os.path.join(RESULTS_DIR, f"lrrfip1_isoseq_isoform_distribution.{ext}")
    if ext == "tiff":
        fig.savefig(out, dpi=600, bbox_inches="tight",
                    pil_kwargs={"compression": "tiff_lzw"})
    else:
        fig.savefig(out, dpi=600, bbox_inches="tight")
    print(f"Saved: {out}")
plt.close()

# ============================================================
# 5. EVIDENCE SUMMARY
# ============================================================
print(f"\n  ┌{'─'*60}┐")
print(f"  │  ISOFORM 1 (MANE SELECT) — ABSENCE EVIDENCE              │")
print(f"  │  ENST00000308482.14 / NM_001137550 / 640aa                │")
print(f"  ├{'─'*60}┤")
print(f"  │  ✗  IsoSeq pigeon v49 filtered    →  0 matches            │")
print(f"  │  ✗  IsoSeq pigeon v49 unfiltered  →  0 matches            │")
print(f"  │  ✗  IsoSeq pigeon v39 filtered    →  0 matches            │")
print(f"  │  ✗  Salmon full transcriptome     →  0 TPM (no k-mers)    │")
print(f"  │  ⚠  Salmon restricted (6 only)    →  forced assignment    │")
print(f"  ├{'─'*60}┤")
print(f"  │  VERDICT: Not expressed in fibroblasts.                    │")
print(f"  │  Illumina detection = restricted quantification artifact.  │")
print(f"  └{'─'*60}┘")

print(f"\n  ▸ Figures saved to: {RESULTS_DIR}/")
print(f"    lrrfip1_isoseq_isoform_distribution.{{png,pdf,svg,tiff}}")
print(f"\n  Done.\n")
