import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import os, glob, warnings
warnings.filterwarnings("ignore")

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
plt.rcParams["font.size"] = 9

OUTDIR = "./data/salmon_quant_iPSC"
IPSC_ROOT = OUTDIR
FIBRO_ROOT = "./data/salmon_quant_fibroblast"
TX2GENE = os.path.join(OUTDIR, "tx2gene_gencode49.tsv")

tx2g = pd.read_csv(TX2GENE, sep="\t", header=None, names=["transcript_id", "gene_name"])
tx_map = dict(zip(tx2g["transcript_id"], tx2g["gene_name"]))

LINKED_GENES = ["COL6A3", "LRRFIP1", "MLPH", "RAB17", "RBM44"]
# PRLH excluded (0 in both tissues)

def load_gene_tpm(quant_files, tx_map, label_func):
    records = []
    for qf in sorted(quant_files):
        sample = label_func(qf)
        df = pd.read_csv(qf, sep="\t")
        df["gene"] = df["Name"].map(tx_map)
        df = df.dropna(subset=["gene"])
        gene_tpm = df.groupby("gene")["TPM"].sum().reset_index()
        gene_tpm["sample"] = sample
        records.append(gene_tpm)
    return pd.concat(records, ignore_index=True)

# ── Load iPSC ──
ipsc_files = glob.glob(os.path.join(IPSC_ROOT, "*/quant.sf"))
ipsc_data = load_gene_tpm(ipsc_files, tx_map, lambda f: os.path.basename(os.path.dirname(f)))
ipsc_data["genotype"] = ipsc_data["sample"].apply(lambda s: "Hetero" if s.startswith("hetero") else "WT")

# Exclude 2 outlier samples
lrrfip1_tpm = ipsc_data[ipsc_data["gene"] == "LRRFIP1"].set_index("sample")["TPM"]
bad_samples = lrrfip1_tpm[lrrfip1_tpm < 10].index.tolist()
ipsc_data = ipsc_data[~ipsc_data["sample"].isin(bad_samples)]

# ── Load fibroblast ──
fibro_files = glob.glob(os.path.join(FIBRO_ROOT, "*/quant.sf"))
fibro_data = load_gene_tpm(fibro_files, tx_map, lambda f: os.path.basename(os.path.dirname(f)))

def get_individual(s):
    if "33016" in s: return "III-13"
    elif "33017" in s: return "II-8.5"
    elif "33018" in s: return "III-12"
    elif "33019" in s: return "II-8"
    return s

fibro_data["individual"] = fibro_data["sample"].apply(get_individual)
fibro_data["genotype"] = fibro_data["individual"].map({
    "III-13": "WT", "II-8.5": "WT", "III-12": "Affected", "II-8": "Affected"
})

# ───────────────────────────────────────────────────────
# FIGURE 1: Linked genes — iPSC vs Fibroblast grouped bars
# ───────────────────────────────────────────────────────
fig1, ax1 = plt.subplots(figsize=(18/2.54, 12/2.54))

x = np.arange(len(LINKED_GENES))
bar_w = 0.35

# Compute means and SDs per tissue
ipsc_means, ipsc_sds = [], []
fibro_means, fibro_sds = [], []
for g in LINKED_GENES:
    ig = ipsc_data[ipsc_data["gene"] == g]["TPM"]
    fg = fibro_data[fibro_data["gene"] == g]["TPM"]
    ipsc_means.append(ig.mean() if len(ig) > 0 else 0)
    ipsc_sds.append(ig.std() if len(ig) > 0 else 0)
    fibro_means.append(fg.mean() if len(fg) > 0 else 0)
    fibro_sds.append(fg.std() if len(fg) > 0 else 0)

ipsc_means = np.array(ipsc_means)
ipsc_sds = np.array(ipsc_sds)
fibro_means = np.array(fibro_means)
fibro_sds = np.array(fibro_sds)

bars_ipsc = ax1.bar(x - bar_w/2, ipsc_means, bar_w, yerr=ipsc_sds, 
                     color="#5e4fa2", alpha=0.85, capsize=3, error_kw={"linewidth": 0.8},
                     edgecolor="white", linewidth=0.5, label="iPSC")
bars_fibro = ax1.bar(x + bar_w/2, fibro_means, bar_w, yerr=fibro_sds,
                      color="#f47d4d", alpha=0.85, capsize=3, error_kw={"linewidth": 0.8},
                      edgecolor="white", linewidth=0.5, label="Fibroblast")

# Annotate values on top
for i, (iv, fv) in enumerate(zip(ipsc_means, fibro_means)):
    if iv > 0.5:
        ax1.text(i - bar_w/2, iv + ipsc_sds[i] + max(fibro_means)*0.02, f"{iv:.1f}",
                 ha="center", va="bottom", fontsize=7, color="#5e4fa2")
    else:
        ax1.text(i - bar_w/2, iv + ipsc_sds[i] + max(fibro_means)*0.02, f"{iv:.1f}",
                 ha="center", va="bottom", fontsize=6.5, color="#5e4fa2")
    if fv > 0.5:
        ax1.text(i + bar_w/2, fv + fibro_sds[i] + max(fibro_means)*0.02, f"{fv:.1f}",
                 ha="center", va="bottom", fontsize=7, color="#c0392b")
    else:
        ax1.text(i + bar_w/2, fv + fibro_sds[i] + max(fibro_means)*0.02, f"{fv:.1f}",
                 ha="center", va="bottom", fontsize=6.5, color="#c0392b")

ax1.set_xticks(x)
ax1.set_xticklabels([f"$\it{{{g}}}$" for g in LINKED_GENES], fontsize=9)
ax1.set_ylabel("Mean TPM (± SD)", fontsize=9)
ax1.spines["top"].set_visible(False)
ax1.spines["right"].set_visible(False)
ax1.legend(fontsize=8, frameon=True, framealpha=0.9, loc="upper right")
ax1.set_title("Linked gene expression — iPSC vs Fibroblast\n(Salmon, GENCODE v49 iso6 index)",
              fontsize=10, fontweight="bold")

plt.tight_layout()
for ext in ["png", "pdf", "svg"]:
    fig1.savefig(os.path.join(OUTDIR, f"linked_genes_barplot.{ext}"),
                 dpi=300, bbox_inches="tight", facecolor="white")
    print(f"Saved: linked_genes_barplot.{ext}")
plt.close()

# ───────────────────────────────────────────────────────
# FIGURE 2: NR3C1 expression — per sample, iPSC vs Fibroblast
# ───────────────────────────────────────────────────────
fig2, (ax2a, ax2b) = plt.subplots(1, 2, figsize=(22/2.54, 13/2.54),
                                    gridspec_kw={"width_ratios": [1.3, 1], "wspace": 0.45})

# --- Panel A: NR3C1 per clone/individual ---
# iPSC: per clone means
ipsc_nr3c1 = ipsc_data[ipsc_data["gene"] == "NR3C1"].copy()
ipsc_nr3c1["clone"] = ipsc_nr3c1["sample"].apply(lambda s: s.split("_Repl")[0])
clone_means = ipsc_nr3c1.groupby(["clone", "genotype"])["TPM"].agg(["mean", "std"]).reset_index()
clone_means["std"] = clone_means["std"].fillna(0)

wt_clones = sorted([c for c in clone_means["clone"] if c.startswith("wt")])
het_clones = sorted([c for c in clone_means["clone"] if c.startswith("hetero")])
ipsc_clone_order = wt_clones + het_clones

clone_labels_ipsc = {
    "wt_1": "WT-1", "wt_2": "WT-2", "wt_3": "WT-3",
    "hetero_1": "Het-1", "hetero_2": "Het-2", "hetero_3": "Het-3",
}

# Fibroblast: per individual means
fibro_nr3c1 = fibro_data[fibro_data["gene"] == "NR3C1"].copy()
indiv_means = fibro_nr3c1.groupby(["individual", "genotype"])["TPM"].agg(["mean", "std"]).reset_index()
indiv_means["std"] = indiv_means["std"].fillna(0)
fibro_indiv_order = ["III-13", "II-8.5", "III-12", "II-8"]

# Combined x-axis labels — pedigree IDs only (genotype shown in sub-headers above)
all_labels = ([clone_labels_ipsc.get(c, c) for c in ipsc_clone_order] +
              fibro_indiv_order)
n_ipsc = len(ipsc_clone_order)
n_fibro = len(fibro_indiv_order)
gap = 1.0
x_ipsc = np.arange(n_ipsc)
x_fibro = np.arange(n_fibro) + n_ipsc + gap
x_all = np.concatenate([x_ipsc, x_fibro])

# Plot iPSC bars
ipsc_vals = []
ipsc_errs = []
ipsc_colors = []
for c in ipsc_clone_order:
    row = clone_means[clone_means["clone"] == c]
    ipsc_vals.append(row["mean"].values[0] if len(row) > 0 else 0)
    ipsc_errs.append(row["std"].values[0] if len(row) > 0 else 0)
    ipsc_colors.append("#888888" if c.startswith("wt") else "#2874A6")

ax2a.bar(x_ipsc, ipsc_vals, 0.65, yerr=ipsc_errs, color=ipsc_colors, alpha=0.8,
         capsize=3, error_kw={"linewidth": 0.8}, edgecolor="white", linewidth=0.5)

# Plot fibroblast bars
fibro_vals = []
fibro_errs = []
fibro_colors = []
for ind in fibro_indiv_order:
    row = indiv_means[indiv_means["individual"] == ind]
    fibro_vals.append(row["mean"].values[0] if len(row) > 0 else 0)
    fibro_errs.append(row["std"].values[0] if len(row) > 0 else 0)
    geno = "WT" if ind in ["III-13", "II-8.5"] else "Affected"
    fibro_colors.append("#888888" if geno == "WT" else "#2874A6")

ax2a.bar(x_fibro, fibro_vals, 0.65, yerr=fibro_errs, color=fibro_colors, alpha=0.8,
         capsize=3, error_kw={"linewidth": 0.8}, edgecolor="white", linewidth=0.5)

# Separator
sep_x = n_ipsc + gap/2 - 0.5
ax2a.axvline(sep_x, color="black", linewidth=1, linestyle="--", alpha=0.5)

# Also within iPSC: WT vs Het separator
ax2a.axvline(len(wt_clones) - 0.5, color="black", linewidth=0.5, linestyle=":", alpha=0.4)
# Within fibroblast
ax2a.axvline(x_fibro[0] + 1.5, color="black", linewidth=0.5, linestyle=":", alpha=0.4)

ax2a.set_xticks(x_all)
ax2a.set_xticklabels(all_labels, fontsize=6.5, rotation=0, ha="center")
ax2a.set_ylabel("$\it{NR3C1}$ TPM", fontsize=9)
ax2a.spines["top"].set_visible(False)
ax2a.spines["right"].set_visible(False)

# Tissue + genotype headers above bars (use axes fraction for y, data coords for x)
ymax = max(max(fibro_vals) + max(fibro_errs), max(ipsc_vals) + max(ipsc_errs)) * 1.25
ax2a.set_ylim(0, ymax)

ax2a.text(np.mean(x_ipsc), ymax * 0.95, "iPSC", ha="center", fontsize=9, fontweight="bold")
ax2a.text(np.mean(x_fibro), ymax * 0.95, "Fibroblast", ha="center", fontsize=9, fontweight="bold")

# WT / variant sub-labels
ax2a.text(np.mean(x_ipsc[:len(wt_clones)]), ymax * 0.88, "WT", ha="center", fontsize=7, color="#888888")
ax2a.text(np.mean(x_ipsc[len(wt_clones):]), ymax * 0.88, "Het", ha="center", fontsize=7, color="#2874A6")
ax2a.text(np.mean(x_fibro[:2]), ymax * 0.88, "healthy", ha="center", fontsize=7, color="#888888")
ax2a.text(np.mean(x_fibro[2:]), ymax * 0.88, "affected", ha="center", fontsize=7, color="#2874A6")

ax2a.text(-0.12, 1.05, "A", transform=ax2a.transAxes, fontsize=14, fontweight="bold", va="top")

# --- Panel B: Tissue comparison (mean ± SD, simple) ---
tissues = ["iPSC", "Fibroblast"]
tissue_means = [np.mean(ipsc_vals), np.mean(fibro_vals)]
tissue_sds = [ipsc_nr3c1["TPM"].std(), fibro_nr3c1["TPM"].std()]

bars = ax2b.bar([0, 1], tissue_means, 0.6, yerr=tissue_sds,
                color=["#5e4fa2", "#f47d4d"], alpha=0.85, capsize=4,
                error_kw={"linewidth": 1}, edgecolor="white", linewidth=0.5)

# Individual data points
np.random.seed(42)
for val in ipsc_nr3c1["TPM"].values:
    ax2b.scatter(0 + np.random.uniform(-0.15, 0.15), val, color="#5e4fa2", 
                 s=15, alpha=0.5, zorder=5, edgecolor="white", linewidth=0.3)
for val in fibro_nr3c1["TPM"].values:
    ax2b.scatter(1 + np.random.uniform(-0.15, 0.15), val, color="#f47d4d",
                 s=15, alpha=0.5, zorder=5, edgecolor="white", linewidth=0.3)

# Annotate means
for i, (m, s) in enumerate(zip(tissue_means, tissue_sds)):
    ax2b.text(i, m + s + 3, f"{m:.1f}", ha="center", va="bottom", fontsize=9, fontweight="bold")

# Fold change annotation
fc = tissue_means[1] / tissue_means[0] if tissue_means[0] > 0 else 0
ax2b.annotate(f"{fc:.0f}×", xy=(0.5, max(tissue_means)*0.7), fontsize=14, fontweight="bold",
              ha="center", color="#c0392b")

ax2b.set_xticks([0, 1])
ax2b.set_xticklabels(tissues, fontsize=9)
ax2b.set_ylabel("$\it{NR3C1}$ TPM", fontsize=9)
ax2b.spines["top"].set_visible(False)
ax2b.spines["right"].set_visible(False)

ax2b.text(-0.15, 1.05, "B", transform=ax2b.transAxes, fontsize=14, fontweight="bold", va="top")

fig2.suptitle("$\it{NR3C1}$ expression — iPSC vs Fibroblast\n(Salmon, GENCODE v49 iso6 index)",
              fontsize=10, fontweight="bold", y=1.06)

for ext in ["png", "pdf", "svg"]:
    fig2.savefig(os.path.join(OUTDIR, f"NR3C1_expression.{ext}"),
                 dpi=300, bbox_inches="tight", facecolor="white")
    print(f"Saved: NR3C1_expression.{ext}")
plt.close()

# Print summary
print(f"\nNR3C1 summary:")
print(f"  iPSC mean: {np.mean(ipsc_vals):.1f} TPM")
print(f"  Fibroblast mean: {np.mean(fibro_vals):.1f} TPM")
print(f"  Fold change: {fc:.1f}×")
