"""
NR3C1 expression figure — iPSC vs Fibroblast
Panel A: per clone/individual bars
Panel B: tissue summary with fold change
Journal specs: 600 dpi, 8pt font, max 18 cm width, TIFF
TPM reported as mean-of-individual/clone-means
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os, glob, warnings
warnings.filterwarnings("ignore")

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Helvetica", "Arial", "DejaVu Sans"]
plt.rcParams["font.size"] = 8

OUTDIR = "./data/salmon_quant_iPSC"
IPSC_ROOT = OUTDIR
FIBRO_ROOT = "./data/salmon_quant_fibroblast"
TX2GENE = os.path.join(OUTDIR, "tx2gene_gencode49.tsv")

tx2g = pd.read_csv(TX2GENE, sep="\t", header=None, names=["transcript_id", "gene_name"])
tx_map = dict(zip(tx2g["transcript_id"], tx2g["gene_name"]))

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

def get_individual(s):
    if "33016" in s: return "III-13"
    elif "33017" in s: return "II-8.5"
    elif "33018" in s: return "III-12"
    elif "33019" in s: return "II-8"
    return s

# Load iPSC
ipsc_files = glob.glob(os.path.join(IPSC_ROOT, "*/quant.sf"))
ipsc_data = load_gene_tpm(ipsc_files, tx_map, lambda f: os.path.basename(os.path.dirname(f)))
ipsc_data["genotype"] = ipsc_data["sample"].apply(lambda s: "Hetero" if s.startswith("hetero") else "WT")

lrrfip1_tpm = ipsc_data[ipsc_data["gene"] == "LRRFIP1"].set_index("sample")["TPM"]
bad_samples = lrrfip1_tpm[lrrfip1_tpm < 10].index.tolist()
ipsc_data = ipsc_data[~ipsc_data["sample"].isin(bad_samples)]

# Load fibroblast
fibro_files = glob.glob(os.path.join(FIBRO_ROOT, "*/quant.sf"))
fibro_data = load_gene_tpm(fibro_files, tx_map, lambda f: os.path.basename(os.path.dirname(f)))
fibro_data["individual"] = fibro_data["sample"].apply(get_individual)
fibro_data["genotype"] = fibro_data["individual"].map({
    "III-13": "WT", "II-8.5": "WT", "III-12": "Affected", "II-8": "Affected"
})

# NR3C1 data
ipsc_nr3c1 = ipsc_data[ipsc_data["gene"] == "NR3C1"].copy()
ipsc_nr3c1["clone"] = ipsc_nr3c1["sample"].apply(lambda s: s.split("_Repl")[0])

fibro_nr3c1 = fibro_data[fibro_data["gene"] == "NR3C1"].copy()

# Clone-level means (iPSC)
clone_means = ipsc_nr3c1.groupby(["clone", "genotype"])["TPM"].agg(["mean", "std"]).reset_index()
clone_means["std"] = clone_means["std"].fillna(0)
wt_clones = sorted([c for c in clone_means["clone"] if c.startswith("wt")])
het_clones = sorted([c for c in clone_means["clone"] if c.startswith("hetero")])
ipsc_clone_order = wt_clones + het_clones
clone_labels_ipsc = {"wt_1": "WT-1", "wt_2": "WT-2", "wt_3": "WT-3",
                     "hetero_1": "Het-1", "hetero_2": "Het-2", "hetero_3": "Het-3"}

# Individual-level means (fibroblast) — Gen order
fibro_indiv_order = ["II-8.5", "II-8", "III-13", "III-12"]
indiv_means = fibro_nr3c1.groupby(["individual", "genotype"])["TPM"].agg(["mean", "std"]).reset_index()
indiv_means["std"] = indiv_means["std"].fillna(0)

# ── Figure ──
fig_w = 18 / 2.54
fig_h = 10 / 2.54

fig, (ax2a, ax2b) = plt.subplots(1, 2, figsize=(fig_w, fig_h),
                                   gridspec_kw={"width_ratios": [1.3, 1], "wspace": 0.45})

# --- Panel A: Per clone/individual ---
all_labels = ([clone_labels_ipsc.get(c, c) for c in ipsc_clone_order] + fibro_indiv_order)
n_ipsc = len(ipsc_clone_order)
n_fibro = len(fibro_indiv_order)
gap = 1.0
x_ipsc = np.arange(n_ipsc)
x_fibro = np.arange(n_fibro) + n_ipsc + gap
x_all = np.concatenate([x_ipsc, x_fibro])

ipsc_vals, ipsc_errs, ipsc_colors = [], [], []
for c in ipsc_clone_order:
    row = clone_means[clone_means["clone"] == c]
    ipsc_vals.append(row["mean"].values[0] if len(row) > 0 else 0)
    ipsc_errs.append(row["std"].values[0] if len(row) > 0 else 0)
    ipsc_colors.append("#888888" if c.startswith("wt") else "#2874A6")

ax2a.bar(x_ipsc, ipsc_vals, 0.65, yerr=ipsc_errs, color=ipsc_colors, alpha=0.8,
         capsize=3, error_kw={"linewidth": 0.8}, edgecolor="white", linewidth=0.5)

fibro_vals, fibro_errs, fibro_colors = [], [], []
for ind in fibro_indiv_order:
    row = indiv_means[indiv_means["individual"] == ind]
    fibro_vals.append(row["mean"].values[0] if len(row) > 0 else 0)
    fibro_errs.append(row["std"].values[0] if len(row) > 0 else 0)
    geno = "WT" if ind in ["III-13", "II-8.5"] else "Affected"
    fibro_colors.append("#888888" if geno == "WT" else "#2874A6")

ax2a.bar(x_fibro, fibro_vals, 0.65, yerr=fibro_errs, color=fibro_colors, alpha=0.8,
         capsize=3, error_kw={"linewidth": 0.8}, edgecolor="white", linewidth=0.5)

# Separators
sep_x = n_ipsc + gap/2 - 0.5
ax2a.axvline(sep_x, color="black", linewidth=1, linestyle="--", alpha=0.5)
ax2a.axvline(len(wt_clones) - 0.5, color="black", linewidth=0.5, linestyle=":", alpha=0.4)
ax2a.axvline(x_fibro[0] + 1.5, color="black", linewidth=0.5, linestyle=":", alpha=0.4)

ax2a.set_xticks(x_all)
ax2a.set_xticklabels(all_labels, fontsize=6, rotation=0, ha="center")
ax2a.set_ylabel("$\\it{NR3C1}$ TPM", fontsize=8)
ax2a.spines["top"].set_visible(False)
ax2a.spines["right"].set_visible(False)

ymax = max(max(fibro_vals) + max(fibro_errs), max(ipsc_vals) + max(ipsc_errs)) * 1.25
ax2a.set_ylim(0, ymax)

ax2a.text(np.mean(x_ipsc), ymax * 0.95, "iPSC", ha="center", fontsize=8, fontweight="bold")
ax2a.text(np.mean(x_fibro), ymax * 0.95, "Fibroblast", ha="center", fontsize=8, fontweight="bold")

ax2a.text(np.mean(x_ipsc[:len(wt_clones)]), ymax * 0.88, "WT", ha="center", fontsize=6.5, color="#888888")
ax2a.text(np.mean(x_ipsc[len(wt_clones):]), ymax * 0.88, "Het", ha="center", fontsize=6.5, color="#2874A6")
ax2a.text(np.mean(x_fibro[:2]), ymax * 0.88, "healthy", ha="center", fontsize=6.5, color="#888888")
ax2a.text(np.mean(x_fibro[2:]), ymax * 0.88, "affected", ha="center", fontsize=6.5, color="#2874A6")

ax2a.text(-0.12, 1.05, "A", transform=ax2a.transAxes, fontsize=12, fontweight="bold", va="top")

# --- Panel B: Tissue comparison ---
# Mean-of-individual/clone-means
ipsc_tissue_mean = np.mean(ipsc_vals)
fibro_tissue_mean = np.mean(fibro_vals)
tissues = ["iPSC", "Fibroblast"]
tissue_means = [ipsc_tissue_mean, fibro_tissue_mean]
tissue_sds = [ipsc_nr3c1["TPM"].std(), fibro_nr3c1["TPM"].std()]

bars = ax2b.bar([0, 1], tissue_means, 0.6, yerr=tissue_sds,
                color=["#5e4fa2", "#f47d4d"], alpha=0.85, capsize=4,
                error_kw={"linewidth": 1}, edgecolor="white", linewidth=0.5)

np.random.seed(42)
for val in ipsc_nr3c1["TPM"].values:
    ax2b.scatter(0 + np.random.uniform(-0.15, 0.15), val, color="#5e4fa2",
                 s=12, alpha=0.5, zorder=5, edgecolor="white", linewidth=0.3)
for val in fibro_nr3c1["TPM"].values:
    ax2b.scatter(1 + np.random.uniform(-0.15, 0.15), val, color="#f47d4d",
                 s=12, alpha=0.5, zorder=5, edgecolor="white", linewidth=0.3)

for i, (m, s) in enumerate(zip(tissue_means, tissue_sds)):
    ax2b.text(i, m + s + 3, f"{m:.1f}", ha="center", va="bottom", fontsize=8, fontweight="bold")

fc = tissue_means[1] / tissue_means[0] if tissue_means[0] > 0 else 0
ax2b.annotate(f"{fc:.0f}×", xy=(0.5, max(tissue_means)*0.7), fontsize=12, fontweight="bold",
              ha="center", color="#c0392b")

ax2b.set_xticks([0, 1])
ax2b.set_xticklabels(tissues, fontsize=8)
ax2b.set_ylabel("$\\it{NR3C1}$ TPM", fontsize=8)
ax2b.spines["top"].set_visible(False)
ax2b.spines["right"].set_visible(False)

ax2b.text(-0.15, 1.05, "B", transform=ax2b.transAxes, fontsize=12, fontweight="bold", va="top")

fig.suptitle("$\\it{NR3C1}$ expression — iPSC vs Fibroblast\n(Salmon, GENCODE v49 iso6 index)",
             fontsize=9, fontweight="bold", y=1.04)

for ext in ["png", "pdf", "svg"]:
    fig.savefig(os.path.join(OUTDIR, f"NR3C1_expression_journal.{ext}"),
                dpi=600, bbox_inches="tight", facecolor="white")
    print(f"Saved: NR3C1_expression_journal.{ext}")
fig.savefig(os.path.join(OUTDIR, "NR3C1_expression_journal.tiff"),
            dpi=600, bbox_inches="tight", facecolor="white",
            pil_kwargs={"compression": "tiff_lzw"})
print("Saved: NR3C1_expression_journal.tiff")

plt.close()

print(f"\nNR3C1 summary:")
print(f"  iPSC mean (of clone means): {ipsc_tissue_mean:.1f} TPM")
print(f"  Fibroblast mean (of individual means): {fibro_tissue_mean:.1f} TPM")
print(f"  Fold change: {fc:.1f}×")
