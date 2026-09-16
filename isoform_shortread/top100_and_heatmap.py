import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import os, glob, warnings
warnings.filterwarnings("ignore")

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans"]
plt.rcParams["font.size"] = 9

OUTDIR = "./data/salmon_quant_iPSC"
IPSC_ROOT = OUTDIR
FIBRO_ROOT = "./data/salmon_quant_fibroblast"
TX2GENE = os.path.join(OUTDIR, "tx2gene_gencode49.tsv")

# Genes of interest for heatmap
HEATMAP_GENES = ["COL6A3", "LRRFIP1", "NR3C1", "MLPH", "PRLH", "RAB17", "RBM44"]

# ── tx2gene ──
tx2g = pd.read_csv(TX2GENE, sep="\t", header=None, names=["transcript_id", "gene_name"])
tx_map = dict(zip(tx2g["transcript_id"], tx2g["gene_name"]))

def load_gene_tpm(quant_files, tx_map, label_func):
    """Load quant.sf files, sum TPM per gene, return DataFrame."""
    records = []
    for qf in sorted(quant_files):
        sample = label_func(qf)
        df = pd.read_csv(qf, sep="\t")
        df["gene"] = df["Name"].map(tx_map)
        # Drop unmapped transcripts (shouldn't be many)
        df = df.dropna(subset=["gene"])
        gene_tpm = df.groupby("gene")["TPM"].sum().reset_index()
        gene_tpm["sample"] = sample
        records.append(gene_tpm)
    return pd.concat(records, ignore_index=True)

# ── iPSC ──
ipsc_files = glob.glob(os.path.join(IPSC_ROOT, "*/quant.sf"))
ipsc_data = load_gene_tpm(ipsc_files, tx_map, lambda f: os.path.basename(os.path.dirname(f)))

# Add genotype
ipsc_data["genotype"] = ipsc_data["sample"].apply(lambda s: "Hetero" if s.startswith("hetero") else "WT")

# Exclude the 2 outlier samples (low LRRFIP1)
lrrfip1_tpm = ipsc_data[ipsc_data["gene"] == "LRRFIP1"].set_index("sample")["TPM"]
bad_samples = lrrfip1_tpm[lrrfip1_tpm < 10].index.tolist()
print(f"iPSC: excluding {len(bad_samples)} samples: {bad_samples}")
ipsc_data = ipsc_data[~ipsc_data["sample"].isin(bad_samples)]

# Gene-level mean TPM across all good samples
ipsc_mean = ipsc_data.groupby("gene")["TPM"].mean().reset_index()
ipsc_mean.columns = ["gene", "mean_TPM"]
ipsc_mean = ipsc_mean.sort_values("mean_TPM", ascending=False).reset_index(drop=True)
ipsc_mean["rank"] = range(1, len(ipsc_mean) + 1)

# Also compute per-genotype means
ipsc_wt = ipsc_data[ipsc_data["genotype"] == "WT"].groupby("gene")["TPM"].mean().reset_index()
ipsc_wt.columns = ["gene", "WT_mean_TPM"]
ipsc_het = ipsc_data[ipsc_data["genotype"] == "Hetero"].groupby("gene")["TPM"].mean().reset_index()
ipsc_het.columns = ["gene", "Het_mean_TPM"]
ipsc_sd = ipsc_data.groupby("gene")["TPM"].std().reset_index()
ipsc_sd.columns = ["gene", "sd_TPM"]

ipsc_table = ipsc_mean.merge(ipsc_wt, on="gene").merge(ipsc_het, on="gene").merge(ipsc_sd, on="gene")
ipsc_top100 = ipsc_table.head(100)

print(f"\niPSC top 5 genes:")
print(ipsc_top100[["rank", "gene", "mean_TPM", "WT_mean_TPM", "Het_mean_TPM"]].head())

# ── Fibroblasts ──
fibro_files = glob.glob(os.path.join(FIBRO_ROOT, "*/quant.sf"))

def fibro_label(f):
    s = os.path.basename(os.path.dirname(f))
    return s

fibro_data = load_gene_tpm(fibro_files, tx_map, fibro_label)

# Add individual/genotype
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

fibro_mean = fibro_data.groupby("gene")["TPM"].mean().reset_index()
fibro_mean.columns = ["gene", "mean_TPM"]
fibro_mean = fibro_mean.sort_values("mean_TPM", ascending=False).reset_index(drop=True)
fibro_mean["rank"] = range(1, len(fibro_mean) + 1)

fibro_wt = fibro_data[fibro_data["genotype"] == "WT"].groupby("gene")["TPM"].mean().reset_index()
fibro_wt.columns = ["gene", "WT_mean_TPM"]
fibro_aff = fibro_data[fibro_data["genotype"] == "Affected"].groupby("gene")["TPM"].mean().reset_index()
fibro_aff.columns = ["gene", "Aff_mean_TPM"]
fibro_sd = fibro_data.groupby("gene")["TPM"].std().reset_index()
fibro_sd.columns = ["gene", "sd_TPM"]

fibro_table = fibro_mean.merge(fibro_wt, on="gene").merge(fibro_aff, on="gene").merge(fibro_sd, on="gene")
fibro_top100 = fibro_table.head(100)

print(f"\nFibroblast top 5 genes:")
print(fibro_top100[["rank", "gene", "mean_TPM", "WT_mean_TPM", "Aff_mean_TPM"]].head())

# ── Check heatmap genes ──
print("\n=== Heatmap genes ===")
for g in HEATMAP_GENES:
    ipsc_row = ipsc_table[ipsc_table["gene"] == g]
    fibro_row = fibro_table[fibro_table["gene"] == g]
    ipsc_val = ipsc_row["mean_TPM"].values[0] if len(ipsc_row) > 0 else "NOT FOUND"
    fibro_val = fibro_row["mean_TPM"].values[0] if len(fibro_row) > 0 else "NOT FOUND"
    ipsc_rank = ipsc_row["rank"].values[0] if len(ipsc_row) > 0 else "N/A"
    fibro_rank = fibro_row["rank"].values[0] if len(fibro_row) > 0 else "N/A"
    print(f"  {g}: iPSC={ipsc_val:.1f} TPM (rank {ipsc_rank}), Fibro={fibro_val:.1f} TPM (rank {fibro_rank})" if isinstance(ipsc_val, float) and isinstance(fibro_val, float) else f"  {g}: iPSC={ipsc_val}, Fibro={fibro_val}")

# Also check RAB17 vs RABB17
for alt in ["RABB17", "RAB7A", "RAB7B"]:
    found = ipsc_table[ipsc_table["gene"] == alt]
    if len(found) > 0:
        print(f"  NOTE: {alt} also found: {found['mean_TPM'].values[0]:.1f} TPM")

# ── Save Excel ──
xlsx_path = os.path.join(OUTDIR, "top100_expressed_genes.xlsx")
with pd.ExcelWriter(xlsx_path, engine="openpyxl") as writer:
    # iPSC sheet
    ipsc_out = ipsc_top100[["rank", "gene", "mean_TPM", "WT_mean_TPM", "Het_mean_TPM", "sd_TPM"]].copy()
    ipsc_out.columns = ["Rank", "Gene", "Mean TPM", "WT Mean TPM", "Hetero Mean TPM", "SD"]
    ipsc_out.to_excel(writer, sheet_name="iPSC Top 100", index=False)
    
    # Fibroblast sheet
    fibro_out = fibro_top100[["rank", "gene", "mean_TPM", "WT_mean_TPM", "Aff_mean_TPM", "sd_TPM"]].copy()
    fibro_out.columns = ["Rank", "Gene", "Mean TPM", "WT Mean TPM", "Affected Mean TPM", "SD"]
    fibro_out.to_excel(writer, sheet_name="Fibroblast Top 100", index=False)
    
    # Heatmap genes sheet (both tissues)
    heatmap_rows = []
    for g in HEATMAP_GENES:
        ipsc_row = ipsc_table[ipsc_table["gene"] == g]
        fibro_row = fibro_table[fibro_table["gene"] == g]
        heatmap_rows.append({
            "Gene": g,
            "iPSC Mean TPM": ipsc_row["mean_TPM"].values[0] if len(ipsc_row) > 0 else 0,
            "iPSC WT": ipsc_row["WT_mean_TPM"].values[0] if len(ipsc_row) > 0 else 0,
            "iPSC Het": ipsc_row["Het_mean_TPM"].values[0] if len(ipsc_row) > 0 else 0,
            "iPSC Rank": int(ipsc_row["rank"].values[0]) if len(ipsc_row) > 0 else "N/A",
            "Fibro Mean TPM": fibro_row["mean_TPM"].values[0] if len(fibro_row) > 0 else 0,
            "Fibro WT": fibro_row["WT_mean_TPM"].values[0] if len(fibro_row) > 0 else 0,
            "Fibro Affected": fibro_row["Aff_mean_TPM"].values[0] if len(fibro_row) > 0 else 0,
            "Fibro Rank": int(fibro_row["rank"].values[0]) if len(fibro_row) > 0 else "N/A",
        })
    hm_df = pd.DataFrame(heatmap_rows)
    hm_df.to_excel(writer, sheet_name="Linked Genes", index=False)

print(f"\nSaved: {xlsx_path}")

# ── Heatmap figure ──
# Build matrix: genes × (iPSC_WT, iPSC_Het, Fibro_WT, Fibro_Aff)
hm_matrix = []
for g in HEATMAP_GENES:
    ipsc_row = ipsc_table[ipsc_table["gene"] == g]
    fibro_row = fibro_table[fibro_table["gene"] == g]
    hm_matrix.append([
        ipsc_row["WT_mean_TPM"].values[0] if len(ipsc_row) > 0 else 0,
        ipsc_row["Het_mean_TPM"].values[0] if len(ipsc_row) > 0 else 0,
        fibro_row["WT_mean_TPM"].values[0] if len(fibro_row) > 0 else 0,
        fibro_row["Aff_mean_TPM"].values[0] if len(fibro_row) > 0 else 0,
    ])

hm_arr = np.array(hm_matrix)
col_labels = ["iPSC\nWT", "iPSC\nHetero", "Fibroblast\nWT", "Fibroblast\nAffected"]

# Log2(TPM+1) for better visualization
hm_log = np.log2(hm_arr + 1)

fig, ax = plt.subplots(figsize=(14/2.54, 12/2.54))

# Custom colormap
cmap = plt.cm.YlOrRd

im = ax.imshow(hm_log, cmap=cmap, aspect="auto")

# Annotate with actual TPM values
for i in range(len(HEATMAP_GENES)):
    for j in range(4):
        val = hm_arr[i, j]
        text_color = "white" if hm_log[i, j] > (hm_log.max() * 0.65) else "black"
        ax.text(j, i, f"{val:.1f}", ha="center", va="center", fontsize=7.5,
                color=text_color, fontweight="bold" if val > 50 else "normal")

ax.set_xticks(range(4))
ax.set_xticklabels(col_labels, fontsize=8)
ax.set_yticks(range(len(HEATMAP_GENES)))
ax.set_yticklabels([f"$\it{{{g}}}$" for g in HEATMAP_GENES], fontsize=9)

# Separator between iPSC and fibroblast columns
ax.axvline(1.5, color="white", linewidth=2.5)

# Colorbar
cbar = fig.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
cbar.set_label("log₂(TPM + 1)", fontsize=8)
cbar.ax.tick_params(labelsize=7)

ax.set_title("Expression of linked genes\niPSC vs Fibroblast (Salmon, GENCODE v49 iso6)",
             fontsize=10, fontweight="bold", pad=10)

plt.tight_layout()

for ext in ["png", "pdf", "svg"]:
    fig.savefig(os.path.join(OUTDIR, f"linked_genes_heatmap.{ext}"),
                dpi=300, bbox_inches="tight", facecolor="white")
    print(f"Saved: linked_genes_heatmap.{ext}")

plt.close()
