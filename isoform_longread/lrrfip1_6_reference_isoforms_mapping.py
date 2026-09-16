"""
LRRFIP1: Map PB isoforms to the 6 selected Gencode v49 reference transcripts.
Produces a summary table showing IsoSeq support for each reference isoform.
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

RESULTS_DIR = os.path.dirname(os.path.abspath(__file__)) if "__file__" in dir() else "."
os.makedirs(RESULTS_DIR, exist_ok=True)

# --- TWEAKABLE PARAMS ---
FIG_W, FIG_H = 18.0 / 2.54, 16.0 / 2.54
FS = 8
FS_SM = 7
FS_LG = 9
FS_XS = 6

# --- MAGMA PALETTE (6 stops, dark to light) ---
MAG = ["#140e36", "#3b0f70", "#721f81", "#b63679", "#f1605d", "#feaf77"]
MAG_FSM = "#2d1160"   # FSM bar color
MAG_ISM = "#b63679"   # ISM bar color
# -------------------------

# ============================================================
# 6 SELECTED GENCODE v49 REFERENCE TRANSCRIPTS
# ============================================================
ref_isoforms = [
    {"name": "Isoform 1", "enst": "ENST00000308482.14", "refseq": "NM_001137550",
     "bp": 4092, "protein": "640aa", "biotype": "Protein coding",
     "flags": "MANE Select, Ensembl Canonical",
     "color": MAG[0]},
    {"name": "Isoform 2", "enst": "ENST00000863709.1", "refseq": "NM_001137551",
     "bp": 3721, "protein": "404aa", "biotype": "Protein coding",
     "flags": "APPRIS P1",
     "color": MAG[1]},
    {"name": "Isoform 3", "enst": "ENST00000392000.4", "refseq": "NM_001137552",
     "bp": 3365, "protein": "808aa", "biotype": "Protein coding",
     "flags": "TSL:1",
     "color": MAG[2]},
    {"name": "Isoform 4", "enst": "ENST00000244815.9", "refseq": "NM_004735",
     "bp": 4370, "protein": "784aa", "biotype": "Protein coding",
     "flags": "GENCODE Primary, TSL:1",
     "color": MAG[3]},
    {"name": "Isoform 5", "enst": "ENST00000289175.10", "refseq": "NM_001137553",
     "bp": 3709, "protein": "752aa", "biotype": "Protein coding",
     "flags": "GENCODE Primary, TSL:1",
     "color": MAG[4]},
    {"name": "Isoform 6", "enst": "ENST00000478958.5", "refseq": "—",
     "bp": 2869, "protein": "No protein", "biotype": "Retained intron",
     "flags": "TSL:1",
     "color": MAG[5]},
]

# ============================================================
# PB ISOFORM MATCHES (from pigeon v49 filtered classification)
# ============================================================
# FSM = full-splice_match, ISM = incomplete-splice_match

pb_matches = {
    "ENST00000308482.14": {
        "fsm": [],
        "ism": [],
        "filtered": [],  # none at all
        "note": "No PB isoforms matched — MANE Select canonical not expressed in IsoSeq",
    },
    "ENST00000863709.1": {
        "fsm": [
            ("PB.4795.1",    2, "alternative_3end"),
            ("PB.4795.4",  418, "alternative_3end"),     # DOMINANT
            ("PB.4795.30",   2, "alternative_3end"),
            ("PB.4795.32",  12, "alternative_3end"),
            ("PB.4795.35",   2, "reference_match"),
            ("PB.4795.38",   2, "alternative_3end"),
            ("PB.4795.58",   8, "alternative_3end"),
        ],
        "ism": [
            ("PB.4795.176", 16, "3prime_fragment"),
        ],
        "filtered": [],
        "note": "7 FSM + 1 ISM — includes dominant PB.4795.4 (FL=418)",
    },
    "ENST00000392000.4": {
        "fsm": [],
        "ism": [],
        "filtered": [
            ("PB.4795.91", "FSM", "IntraPriming"),
            ("PB.4795.159", "FSM", "IntraPriming"),
            ("PB.4795.205", "ISM", "IntraPriming"),
        ],
        "note": "3 matches filtered (all IntraPriming artifacts)",
    },
    "ENST00000244815.9": {
        "fsm": [],
        "ism": [],
        "filtered": [
            ("PB.4795.110", "FSM", "IntraPriming"),
            ("PB.4795.160", "FSM", "IntraPriming"),
            ("PB.4795.124", "ISM", "IntraPriming"),
            # Plus ~15 mono-exon ISM artifacts
        ],
        "note": "All matches filtered (IntraPriming) — many mono-exon artifacts",
    },
    "ENST00000289175.10": {
        "fsm": [
            ("PB.4795.94",   5, "reference_match"),
            ("PB.4795.108",  2, "alternative_3end"),
            ("PB.4795.154",  4, "alternative_5end"),
        ],
        "ism": [],
        "filtered": [],
        "note": "3 FSM — consistent across v39 and v49",
    },
    "ENST00000478958.5": {
        "fsm": [
            ("PB.4795.284", 73, "alternative_5end"),
            ("PB.4795.285",  2, "alternative_3end5end"),
            ("PB.4795.288",  4, "alternative_3end5end"),
            ("PB.4795.290",  2, "alternative_5end"),
        ],
        "ism": [
            ("PB.4795.293", 37, "3prime_fragment"),
            ("PB.4795.295",  2, "3prime_fragment"),
            ("PB.4795.298", 21, "3prime_fragment"),
        ],
        "filtered": [],
        "note": "4 FSM + 3 ISM — retained intron isoform with strong IsoSeq support",
    },
}

# ============================================================
# FIGURE
# ============================================================

fig, axes = plt.subplots(2, 1, figsize=(FIG_W, FIG_H),
                          gridspec_kw={"height_ratios": [1, 1.2]})

# ---- PANEL A: IsoSeq support per reference isoform (bar chart) ----
ax1 = axes[0]
ax1.set_title("A  IsoSeq support for 6 LRRFIP1 reference transcripts (Gencode v49)",
              fontsize=FS_LG, fontweight="bold", loc="left", pad=8)

x = np.arange(6)
bar_w = 0.35

fsm_counts = []
ism_counts = []
total_fl = []
for ref in ref_isoforms:
    enst = ref["enst"]
    matches = pb_matches[enst]
    n_fsm = len(matches["fsm"])
    n_ism = len(matches["ism"])
    fl_sum = sum(fl for _, fl, _ in matches["fsm"]) + sum(fl for _, fl, _ in matches["ism"])
    fsm_counts.append(n_fsm)
    ism_counts.append(n_ism)
    total_fl.append(fl_sum)

# Stacked bars: FSM + ISM
ax1.bar(x, fsm_counts, bar_w, color=MAG_FSM, edgecolor="white",
        linewidth=0.5, label="Full-splice match (FSM)")
ax1.bar(x, ism_counts, bar_w, bottom=fsm_counts, color=MAG_ISM,
        edgecolor="white", linewidth=0.5, label="Incomplete-splice match (ISM)")

# FL count annotation on top
for i in range(6):
    total = fsm_counts[i] + ism_counts[i]
    if total > 0:
        ax1.text(i, total + 0.3, f"FL={total_fl[i]}",
                 ha="center", va="bottom", fontsize=FS_SM, fontweight="bold",
                 color=MAG_FSM)
    else:
        ax1.text(i, 0.15, "0", ha="center", va="bottom", fontsize=FS_SM,
                 color="#999")

# Labels inside bars
for i in range(6):
    if fsm_counts[i] > 0:
        ax1.text(i, fsm_counts[i] / 2, str(fsm_counts[i]),
                 ha="center", va="center", fontsize=FS, fontweight="bold",
                 color="white")
    if ism_counts[i] > 0:
        ax1.text(i, fsm_counts[i] + ism_counts[i] / 2, str(ism_counts[i]),
                 ha="center", va="center", fontsize=FS, fontweight="bold",
                 color="#333")

ax1.set_xticks(x)
xlabels = []
for ref in ref_isoforms:
    xlabels.append(f"{ref['name']}\n{ref['enst'].split('.')[0]}")
ax1.set_xticklabels(xlabels, fontsize=FS_XS, ha="center")
ax1.set_ylabel("PB isoforms matched", fontsize=FS)
ax1.set_ylim(0, 12)
ax1.tick_params(axis="y", labelsize=FS_SM)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)

ax1.legend(fontsize=FS_SM, loc="upper left", frameon=True,
           framealpha=0.9, edgecolor="#ccc")

# Mark filtered isoforms
for i, ref in enumerate(ref_isoforms):
    enst = ref["enst"]
    n_filt = len(pb_matches[enst]["filtered"])
    if n_filt > 0 and fsm_counts[i] == 0 and ism_counts[i] == 0:
        ax1.text(i, 0.5, f"({n_filt} filtered:\nIntraPriming)",
                 ha="center", va="bottom", fontsize=FS_XS, color="#999",
                 style="italic")

# ---- PANEL B: Detailed mapping table ----
ax2 = axes[1]
ax2.set_xlim(0, 10)
ax2.set_ylim(-1.5, 7.5)
ax2.axis("off")
ax2.set_title("B  Reference transcript details and IsoSeq classification",
              fontsize=FS_LG, fontweight="bold", loc="left", pad=8)

# Table header
header_y = 7.0
headers = ["Ref. isoform", "ENST ID", "bp", "Biotype", "FSM", "ISM",
           "Total FL", "Key PB isoform"]
header_x = [0.0, 1.4, 3.5, 4.2, 5.3, 6.0, 6.7, 7.6]
for hx, htext in zip(header_x, headers):
    ax2.text(hx, header_y, htext, fontsize=FS_XS, fontweight="bold",
             color="#333", va="center")
ax2.plot([0, 10], [header_y - 0.3, header_y - 0.3], color="#999", linewidth=0.5)

for i, ref in enumerate(ref_isoforms):
    y = header_y - 1.0 - i * 1.15
    enst = ref["enst"]
    matches = pb_matches[enst]

    # Alternating bg
    if i % 2 == 0:
        ax2.fill_between([0, 10], y - 0.5, y + 0.5,
                         color="#f5f5f5", alpha=0.5, zorder=0)

    # Colored indicator
    ax2.add_patch(plt.Rectangle((0, y - 0.35), 0.12, 0.7,
                                 facecolor=ref["color"], edgecolor="none", zorder=2))

    # Ref isoform name + flags
    ax2.text(header_x[0] + 0.2, y + 0.12, ref["name"],
             fontsize=FS_SM, va="center", fontweight="bold")
    ax2.text(header_x[0] + 0.2, y - 0.18, ref["flags"],
             fontsize=FS_XS - 0.5, va="center", color="#888", style="italic")

    # ENST ID
    ax2.text(header_x[1], y, ref["enst"], fontsize=FS_XS, va="center",
             fontfamily="monospace")

    # bp
    ax2.text(header_x[2], y, str(ref["bp"]), fontsize=FS_SM, va="center")

    # Biotype
    bt_color = "#333" if ref["biotype"] == "Protein coding" else "#c44"
    ax2.text(header_x[3], y, ref["biotype"], fontsize=FS_XS, va="center",
             color=bt_color)

    # FSM count
    n_fsm = len(matches["fsm"])
    fsm_col = MAG_FSM if n_fsm > 0 else "#ccc"
    ax2.text(header_x[4], y, str(n_fsm), fontsize=FS, va="center",
             fontweight="bold", color=fsm_col, ha="center")

    # ISM count
    n_ism = len(matches["ism"])
    ism_col = MAG_ISM if n_ism > 0 else "#ccc"
    ax2.text(header_x[5], y, str(n_ism), fontsize=FS, va="center",
             fontweight="bold", color=ism_col, ha="center")

    # Total FL
    fl_sum = sum(fl for _, fl, _ in matches["fsm"]) + sum(fl for _, fl, _ in matches["ism"])
    fl_col = MAG_FSM if fl_sum > 0 else "#ccc"
    ax2.text(header_x[6], y, str(fl_sum) if fl_sum > 0 else "0",
             fontsize=FS, va="center", fontweight="bold", color=fl_col, ha="center")

    # Key PB isoform
    if matches["fsm"]:
        top_pb = max(matches["fsm"], key=lambda x: x[1])
        ax2.text(header_x[7], y,
                 f"{top_pb[0]} (FL={top_pb[1]})",
                 fontsize=FS_XS, va="center", color=MAG_FSM)
    elif matches["filtered"]:
        ax2.text(header_x[7], y,
                 f"filtered ({len(matches['filtered'])} IntraPriming)",
                 fontsize=FS_XS, va="center", color="#999", style="italic")
    else:
        ax2.text(header_x[7], y, "no match", fontsize=FS_XS,
                 va="center", color="#999", style="italic")

# Summary box
ax2.text(0, -1.2,
         "Isoform 2 (ENST00000863709.1, APPRIS P1): strongest IsoSeq support — "
         "7 FSM matches including dominant PB.4795.4 (FL=418)\n"
         "Isoform 1 (MANE Select): no IsoSeq reads matched — may not be expressed "
         "in this tissue/context\n"
         "Isoforms 3 & 4: only IntraPriming artifacts — no genuine IsoSeq support",
         fontsize=FS_SM, va="top", color="#555",
         bbox=dict(boxstyle="round,pad=0.4", facecolor="#f0f0f8",
                   edgecolor=MAG_FSM, alpha=0.8))

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    plt.tight_layout()

for ext in ("png", "pdf", "svg", "tiff"):
    out = os.path.join(RESULTS_DIR, f"lrrfip1_6_reference_isoforms_mapping.{ext}")
    if ext == "tiff":
        fig.savefig(out, dpi=600, bbox_inches="tight",
                    pil_kwargs={"compression": "tiff_lzw"})
    else:
        fig.savefig(out, dpi=600, bbox_inches="tight")
plt.close()

# ============================================================
# Presentation-style summary
# ============================================================
print(f"\n  {'─'*60}")
print(f"  │  LRRFIP1 — 6 Reference Isoforms IsoSeq Mapping (Magma)")
print(f"  {'─'*60}\n")

print(f"  ┌{'─'*72}┐")
print(f"  │  {'Isoform':<11} {'ENST ID':<23} {'bp':>5} {'FSM':>4} {'ISM':>4} {'FL':>5}  {'Status':<18} │")
print(f"  ├{'─'*72}┤")
for ref in ref_isoforms:
    enst = ref["enst"]
    m = pb_matches[enst]
    n_fsm = len(m["fsm"])
    n_ism = len(m["ism"])
    fl = sum(fl for _, fl, _ in m["fsm"]) + sum(fl for _, fl, _ in m["ism"])
    n_filt = len(m["filtered"])
    if n_fsm > 0 or n_ism > 0:
        status = f"{n_fsm} FSM + {n_ism} ISM"
    elif n_filt > 0:
        status = f"{n_filt} filt (IntraPrm)"
    else:
        status = "no match"
    print(f"  │  {ref['name']:<11} {enst:<23} {ref['bp']:>5} {n_fsm:>4} {n_ism:>4} "
          f"{fl:>5}  {status:<18} │")
print(f"  └{'─'*72}┘")

print(f"\n  ▸ Figures saved to: {RESULTS_DIR}/")
print(f"    lrrfip1_6_reference_isoforms_mapping.{{png,pdf,svg,tiff}}")
print(f"\n  Done.\n")
