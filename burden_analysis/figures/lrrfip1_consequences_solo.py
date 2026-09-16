"""
LRRFIP1 solo VEP consequences figure.
Panel A: Log-scale grouped bars (case-only / shared / control-only)
Panel B: Proportional stacked bars (coding+regulatory only, excluding intronic)
         Stacked bottom→top: LoF, Splice, Missense, In-frame indel, UTR, Synonymous
         Cobalt blue palette (darkest = most severe)
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_theme(style="white", font_scale=1.1)

RESULTS_DIR = os.path.dirname(os.path.abspath(__file__)) if "__file__" in dir() else "."

# ── DATA ──────────────────────────────────────────────────────────────────────
categories = ['LoF', 'Missense', 'Synonymous', 'Splice\nregion', 'UTR', 'In-frame\nindel', 'Intronic']
case_only =  [1,   1,   1,   0,   5,   0,   118]
shared =     [0,   1,   3,   1,   2,   0,   253]
ctrl_only =  [47,  230, 118, 106, 458, 5,   29636]

fig, axes = plt.subplots(1, 2, figsize=(16, 8))

# ── Panel A: Consequences by carrier status (log scale) ───────────────────────
ax = axes[0]
x = np.arange(len(categories))
width = 0.30

bars_cu = ax.bar(x - width, case_only, width, color='#003fa8', edgecolor='white', linewidth=0.5, label='Case-only')
bars_sh = ax.bar(x, shared, width, color='#668ccb', edgecolor='white', linewidth=0.5, label='Shared')
bars_ko = ax.bar(x + width, ctrl_only, width, color='#ccd9ee', edgecolor='white', linewidth=0.5, label='Control-only')

ax.set_yscale('log')
ax.set_ylim(0.5, 50000)
ax.set_xticks(x)
ax.set_xticklabels(categories, fontsize=12, rotation=45, ha='right')
ax.set_ylabel('Number of unique variants (log scale)', fontsize=16)
ax.set_title('LRRFIP1 - VEP consequences\nby carrier status', fontsize=18, fontweight='bold')
ax.legend(fontsize=12, loc='upper left')
ax.tick_params(axis='x', pad=1)
ax.tick_params(axis='y', pad=2)

for bars in [bars_cu, bars_sh, bars_ko]:
    for bar in bars:
        h = bar.get_height()
        if h > 0:
            ax.text(bar.get_x() + bar.get_width()/2, h*1.1, str(int(h)),
                   ha='center', va='bottom', fontsize=10, color='#444', fontweight='bold')

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# ── Panel B: Proportional view (coding+regulatory only) ───────────────────────
ax2 = axes[1]

# Original data: [LoF, Missense, Synonymous, Splice, UTR, In-frame indel]
cu_focus_orig = [1, 1, 1, 0, 5, 0]
sh_focus_orig = [0, 1, 3, 1, 2, 0]
ko_focus_orig = [47, 230, 118, 106, 458, 5]

# Reorder for stacking (bottom→top): LoF, Splice, Missense, In-frame indel, UTR, Synonymous
stack_order = [0, 3, 1, 5, 4, 2]  # indices into original data

cu_focus = [cu_focus_orig[i] for i in stack_order]
sh_focus = [sh_focus_orig[i] for i in stack_order]
ko_focus = [ko_focus_orig[i] for i in stack_order]

cu_total = sum(cu_focus)
sh_total = sum(sh_focus)
ko_total = sum(ko_focus)

cu_pct = [v/cu_total*100 if cu_total > 0 else 0 for v in cu_focus]
sh_pct = [v/sh_total*100 if sh_total > 0 else 0 for v in sh_focus]
ko_pct = [v/ko_total*100 if ko_total > 0 else 0 for v in ko_focus]

x2 = np.arange(3)
width2 = 0.55

# Cobalt blue palette: darkest (LoF) → lightest (Synonymous)
colors_stack = [
    '#003fa8',   # LoF — darkest cobalt
    '#1a5bbf',   # Splice region
    '#3365b9',   # Missense
    '#668ccb',   # In-frame indel
    '#5580cc',   # UTR
    '#99b2dc',   # Synonymous — lightest
]

labels_stack = ['LoF', 'Splice region', 'Missense', 'In-frame indel', 'UTR', 'Synonymous']

groups = [cu_pct, sh_pct, ko_pct]
group_labels = [f'Case-only\n(n={cu_total})', f'Shared\n(n={sh_total})', f'Control-only\n(n={ko_total})']

for i, (label, color) in enumerate(zip(labels_stack, colors_stack)):
    vals = [g[i] for g in groups]
    ax2.bar(x2, vals, width2, bottom=[sum(groups[j][:i]) for j in range(3)],
            color=color, edgecolor='white', linewidth=0.5, label=label)

ax2.set_xticks(x2)
ax2.set_xticklabels(group_labels, fontsize=11)
ax2.set_ylabel('Proportion of coding+regulatory variants (%)', fontsize=16)
ax2.set_ylim(0, 105)
ax2.set_title('LRRFIP1 - Coding+regulatory composition\n(excluding intronic)', fontsize=18, fontweight='bold')

# Legend: inverted so it reads top-to-bottom matching the visual stack
handles, lbls = ax2.get_legend_handles_labels()
ax2.legend(handles[::-1], lbls[::-1], fontsize=12, loc='upper right', bbox_to_anchor=(1.35, 1))

ax2.tick_params(axis='x', pad=2)
ax2.tick_params(axis='y', pad=2)

# Percentage labels — show if >= 3%
for j in range(3):
    bottom = 0
    for i in range(len(labels_stack)):
        val = groups[j][i]
        if val >= 3:
            txt_color = 'white' if colors_stack[i] in ['#003fa8', '#1a5bbf', '#3365b9'] else 'black'
            ax2.text(x2[j], bottom + val/2, f'{val:.0f}%', ha='center', va='center',
                    fontsize=12, fontweight='bold', color=txt_color)
        bottom += val

ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)

# Key finding
ax2.text(0, -18, 'Case-only: 63% UTR + 13% LoF + 13% Missense - signal concentrated in regulatory + damaging variants',
         fontsize=12, color='#003fa8', style='italic')

plt.tight_layout()

# Uncomment to save:
# for ext in ('png', 'pdf', 'svg'):
#                 dpi=200, bbox_inches='tight', facecolor='white')

plt.show()
