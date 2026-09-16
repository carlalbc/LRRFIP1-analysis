"""
LRRFIP1 burden signal — variant composition by region mask.
Shows case-only, shared, and control-only variants per mask
with consequence breakdown and burden test OR/p.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

sns.set_theme(style="white", font_scale=1.1)

RESULTS_DIR = os.path.dirname(os.path.abspath(__file__)) if "__file__" in dir() else "."

# ── Data ──────────────────────────────────────────────────────────────────────
cons_data = {
    'LoF':           {'case': 1, 'shared': 0, 'ctrl': 47},
    'Missense':      {'case': 1, 'shared': 1, 'ctrl': 230},
    'Synonymous':    {'case': 1, 'shared': 3, 'ctrl': 118},
    'Splice region': {'case': 0, 'shared': 1, 'ctrl': 106},
    'In-frame indel':{'case': 0, 'shared': 0, 'ctrl': 5},
    'UTR':           {'case': 5, 'shared': 2, 'ctrl': 458},
    'Intronic':      {'case': 118, 'shared': 253, 'ctrl': 29636},
}

CDS_CONS = ['LoF', 'Missense', 'Synonymous', 'Splice region', 'In-frame indel']
UTR_CONS = ['UTR']
INTRONIC_CONS = ['Intronic']
CONS_ORDER = ['LoF', 'Missense', 'Synonymous', 'Splice region', 'UTR', 'Intronic']

CONS_COLORS = {
    'LoF':            '#011f4b',
    'Missense':       '#03396c',
    'Synonymous':     '#6497b1',
    'Splice region':  '#343d46',
    'In-frame indel': '#b3cde0',
    'UTR':            '#4c516d',
    'Intronic':       '#a8b2bd',
}

masks = [
    {'label': 'CDS', 'cons': CDS_CONS,
     'OR': 3.43, 'CI': '0.6–18.7', 'p': 1.5e-1, 'sig': False},
    {'label': 'CDS + UTR', 'cons': CDS_CONS + UTR_CONS,
     'OR': 2.66, 'CI': '1.3–5.3', 'p': 1.5e-2, 'sig': True},
    {'label': 'UTR + promoter', 'cons': UTR_CONS,
     'OR': 2.45, 'CI': '1.1–5.4', 'p': 5.0e-2, 'sig': True},
    {'label': 'Intronic', 'cons': INTRONIC_CONS,
     'OR': 0.99, 'CI': '0.8–1.2', 'p': 9.2e-1, 'sig': False},
    {'label': 'CDS + UTR + intronic', 'cons': CDS_CONS + UTR_CONS + INTRONIC_CONS,
     'OR': 1.02, 'CI': '0.8–1.2', 'p': 8.5e-1, 'sig': False},
    {'label': 'All regions', 'cons': CDS_CONS + UTR_CONS + INTRONIC_CONS,
     'OR': 1.02, 'CI': '0.9–1.2', 'p': 8.3e-1, 'sig': False},
]

for m in masks:
    for status in ['case', 'shared', 'ctrl']:
        counts = {c: cons_data[c][status] if c in m['cons'] else 0 for c in CONS_ORDER}
        m[f'{status}_counts'] = counts
        m[f'{status}_total'] = sum(counts.values())

# ── Figure ────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(16, 10.5))

BAR_WIDTH = 100
row_h = 0.38
sub_gap = 0.06
group_gap = 0.55
statuses = [('case', 'case-only'), ('shared', 'shared'), ('ctrl', 'control-only')]
STATUS_COLORS = {'case': '#011f4b', 'shared': '#03396c', 'ctrl': '#6497b1'}

y_positions = {}
for i in range(len(masks)):
    y_base = -(i * (3 * row_h + 2 * sub_gap + group_gap))
    y_positions[i] = {
        'case': y_base,
        'shared': y_base - row_h - sub_gap,
        'ctrl': y_base - 2 * (row_h + sub_gap),
        'center': y_base - row_h - sub_gap,
    }

for i, m in enumerate(masks):
    for status, status_label in statuses:
        y = y_positions[i][status]
        total = m[f'{status}_total']
        if total == 0:
            ax.barh(y, BAR_WIDTH, height=row_h, color='none',
                    edgecolor='#ddd', linewidth=0.5, zorder=2)
            ax.text(BAR_WIDTH + 2, y, f"n = 0  {status_label}",
                    ha='left', va='center', fontsize=12, color='#bbb')
            continue

        left = 0
        for c in CONS_ORDER:
            val = m[f'{status}_counts'][c]
            if val == 0:
                continue
            pct = val / total * BAR_WIDTH
            ax.barh(y, pct, height=row_h, left=left,
                    color=CONS_COLORS[c], alpha=1.0,
                    edgecolor='white', linewidth=0.8, zorder=3)
            if pct >= 8:
                ax.text(left + pct / 2, y, f"{pct:.0f}%",
                        ha='center', va='center', fontsize=12,
                        color='white', fontweight='bold', zorder=4)
            left += pct

        ax.text(BAR_WIDTH + 2, y, f"n = {total:,}  {status_label}",
                ha='left', va='center', fontsize=12,
                color=STATUS_COLORS[status],
                fontweight='bold' if status == 'case' else 'normal')

    yc = y_positions[i]['center']
    ax.text(-2, yc, m['label'], ha='right', va='center',
            fontsize=13, fontweight='bold', color='#222')

    is_sig = m['sig']
    or_color = '#c0392b' if is_sig else '#888'
    p_str = f"p = {m['p']:.1e}" if m['p'] < 0.01 else f"p = {m['p']:.2f}"
    sig_mark = "  *" if is_sig else ""
    arrow = "  ▲" if m['OR'] > 1 and is_sig else ""

    or_x = BAR_WIDTH + 38
    ax.text(or_x, yc + 0.12,
            f"OR = {m['OR']:.2f}  ({m['CI']})",
            ha='left', va='center', fontsize=12,
            fontweight='bold' if is_sig else 'normal', color=or_color)
    ax.text(or_x, yc - 0.28,
            f"{p_str}{sig_mark}{arrow}",
            ha='left', va='center', fontsize=12, color=or_color,
            fontweight='bold' if is_sig else 'normal')

# Title
ax.set_title("LRRFIP1 — Burden signal variants by region mask\n"
             "Consequence composition by carrier status  ·  MAC ≤ 5 singleton  ·  ADD burden",
             fontsize=14, fontweight='bold', pad=12, linespacing=1.6)

# Bottom annotation
bottom_y = y_positions[5]['ctrl'] - row_h - 0.7
ax.text(BAR_WIDTH / 2, bottom_y,
        "Signal concentrated in coding + regulatory regions — diluted by ~30,000 intronic variants (OR ≈ 1.0)",
        ha='center', va='top', fontsize=11, color='#1B5E20',
        fontweight='bold', style='italic')

# Legend
legend_elements = [
    mpatches.Patch(color='#011f4b', label='LoF'),
    mpatches.Patch(color='#03396c', label='Missense'),
    mpatches.Patch(color='#6497b1', label='Synonymous'),
    mpatches.Patch(color='#343d46', label='Splice region'),
    mpatches.Patch(color='#4c516d', label='UTR'),
    mpatches.Patch(color='#a8b2bd', label='Intronic'),
]
leg = ax.legend(handles=legend_elements, loc='upper center', fontsize=12,
          frameon=True, framealpha=0.95, edgecolor='#ccc', ncol=6,
          bbox_to_anchor=(0.42, 1.0))
leg.get_frame().set_linewidth(0.8)

# Cleanup
ax.set_xlim(-3, BAR_WIDTH + 85)
ax.set_ylim(bottom_y - 0.4, y_positions[0]['case'] + row_h + 0.8)
ax.set_yticks([])
ax.set_xticks([])
for spine in ax.spines.values():
    spine.set_visible(False)

plt.tight_layout()

# Uncomment to save:
# for ext in ('png', 'pdf', 'svg'):
#                 dpi=200, bbox_inches='tight', facecolor='white')

plt.show()
