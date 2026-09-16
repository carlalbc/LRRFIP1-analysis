#!/usr/bin/env bash
#===============================================================================
# Regenie SKAT-O Pipeline — Scoliosis Rare Variant Analysis
# 6 candidate genes on chr2, UK Biobank WGS (490,549 samples)
#
# Uses GENCODE GTF functional sub-regions for mask definitions.
# AAF bins: 0.001, 0.01. Tests: skato, skat, burden, acatv.
# Separate step 2 runs at different AAF thresholds for comparison.
#
# Instance: mem3_ssd1_v2_x16 (16 cores, 128GB RAM, SSD)
# Disk:     Request ≥1TB when launching
#
# Run each phase sequentially. The script is written in blocks so we can
# restart from any phase if something fails.
#===============================================================================
set -euo pipefail

#-------------------------------------------------------------------------------
# CONFIGURATION — edit these to match our DNAnexus environment
#-------------------------------------------------------------------------------

THREADS=16
PROJECT="scoliosis - Dec 17, 2024"

# Array genotype files (field 22418)
ARRAY_PREFIX="ukb22418"

# Phenotype and covariate files
PHENO_FILE="pheno_500k.regenie.txt"
PHENO_COL="scoliosis"
COVAR_FILE="covariates_p22009_500k_imputed.regenie.txt"
COVAR_COLS="p21022,p22009_a1,p22009_a2,p22009_a3,p22009_a4,p22009_a5"
CAT_COVAR="p31"

# pVCF blocks for the 6-gene region (chr2:237,324,003–237,842,808)
BLOCK_START=11865
BLOCK_END=11892

# Output prefixes
ARRAY_MERGED="ukb22418_allChrs"
WGS_RAW="wgs_6genes_raw"
WGS_CLEAN="wgs_6genes_clean"
STEP1_OUT="step1_scoliosis"
STEP2_PREFIX="scoliosis_skato"

# GENCODE GTF for functional annotation
GTF_FILE="gencode.v49.annotation.gtf"

# Annotation prefix
ANNOT_PREFIX="regenie_annotations"

#===============================================================================
# PHASE 1a: Download & merge array genotypes (all 22 autosomes)
#===============================================================================
echo "===== PHASE 1a: Array genotype merge ====="

# Download all per-chromosome array files
for chr in $(seq 1 22); do
    dx download "project-GxQkyKQJJPBxY3K442z9Qp4V:/Bulk/Genotype Results/Genotype calls/${ARRAY_PREFIX}_c${chr}_b0_v2.bed" &
    dx download "project-GxQkyKQJJPBxY3K442z9Qp4V:/Bulk/Genotype Results/Genotype calls/${ARRAY_PREFIX}_c${chr}_b0_v2.bim" &
    dx download "project-GxQkyKQJJPBxY3K442z9Qp4V:/Bulk/Genotype Results/Genotype calls/${ARRAY_PREFIX}_c${chr}_b0_v2.fam" &
done
wait
echo "Array downloads complete."

# Build merge list (chromosomes 2-22)
rm -f merge_list.txt
for chr in $(seq 2 22); do
    echo "${ARRAY_PREFIX}_c${chr}_b0_v2.bed ${ARRAY_PREFIX}_c${chr}_b0_v2.bim ${ARRAY_PREFIX}_c${chr}_b0_v2.fam" >> merge_list.txt
done

# Merge all chromosomes into one file
/opt/conda/bin/plink \
    --bed ${ARRAY_PREFIX}_c1_b0_v2.bed \
    --bim ${ARRAY_PREFIX}_c1_b0_v2.bim \
    --fam ${ARRAY_PREFIX}_c1_b0_v2.fam \
    --merge-list merge_list.txt \
    --make-bed \
    --out ${ARRAY_MERGED}

echo "Array merge complete: ${ARRAY_MERGED}.bed/bim/fam"

#===============================================================================
# PHASE 1b: Download pVCF blocks for 6-gene region (28 blocks)
# Sequential downloads to avoid overwhelming the instance.
# Run with nohup from a terminal:
#   nohup bash regenie_pipeline.sh > pipeline.log 2>&1 &
#===============================================================================
echo "===== PHASE 1b: pVCF download ====="

for b in $(seq ${BLOCK_START} ${BLOCK_END}); do
    echo "$(date) — Downloading block ${b} (vcf.gz + tbi)..."
    dx download "project-GxQkyKQJJPBxY3K442z9Qp4V:/Bulk/GATK and GraphTyper WGS/GraphTyper population level WGS variants, pVCF format - 500k release/chr2/ukb23374_c2_b${b}_v1.vcf.gz"
    dx download "project-GxQkyKQJJPBxY3K442z9Qp4V:/Bulk/GATK and GraphTyper WGS/GraphTyper population level WGS variants, pVCF format - 500k release/chr2/ukb23374_c2_b${b}_v1.vcf.gz.tbi"
    echo "$(date) — Block ${b} done."
done
echo "pVCF downloads complete (blocks ${BLOCK_START}-${BLOCK_END})."

#===============================================================================
# PHASE 1c: Merge pVCFs → clean BED (allele-specific IDs, no * alleles)
#===============================================================================
echo "===== PHASE 1c: pVCF merge & BED creation ====="

# Build file list for bcftools (brace expansion doesn't work with variables)
for b in $(seq ${BLOCK_START} ${BLOCK_END}); do
    echo "ukb23374_c2_b${b}_v1.vcf.gz"
done > vcf_list.txt

# Concatenate all blocks (--naive: fast, no decompress/recompress)
bcftools concat --naive --file-list vcf_list.txt -Oz -o 6genes_region.vcf.gz
tabix -p vcf 6genes_region.vcf.gz

# Free disk: remove individual blocks (~308GB)
for b in $(seq ${BLOCK_START} ${BLOCK_END}); do
    rm -f "ukb23374_c2_b${b}_v1.vcf.gz" "ukb23374_c2_b${b}_v1.vcf.gz.tbi"
done
echo "Individual blocks removed to free disk."

# Convert to BED with allele-specific IDs (chr:pos:ref:alt)
# NOTE: plink2 puts ALT in A1 (BIM col 5), REF in A2 (BIM col 6)
plink2 --vcf 6genes_region.vcf.gz \
    --max-alleles 2 \
    --set-all-var-ids @:#:\$r:\$a \
    --new-id-max-allele-len 500 \
    --make-bed \
    --threads ${THREADS} \
    --out ${WGS_RAW}

# Remove * alleles (deletion placeholders from multiallelic split)
awk '$5 != "*" && $6 != "*" {print $2}' ${WGS_RAW}.bim > keep_vars.txt
plink2 --bfile ${WGS_RAW} \
    --extract keep_vars.txt \
    --make-bed \
    --threads ${THREADS} \
    --out ${WGS_CLEAN}

# Free disk: remove raw BED (keeping only clean)
rm -f ${WGS_RAW}.bed ${WGS_RAW}.bim ${WGS_RAW}.fam

# Verify no duplicate variant IDs
n_dups=$(awk '{print $2}' ${WGS_CLEAN}.bim | sort | uniq -d | wc -l)
if [ "$n_dups" -gt 0 ]; then
    echo "WARNING: ${n_dups} duplicate variant IDs found in ${WGS_CLEAN}.bim"
    awk '{print $2}' ${WGS_CLEAN}.bim | sort | uniq -d > duplicate_ids.txt
    echo "Duplicate IDs written to duplicate_ids.txt — inspect before continuing."
    exit 1
fi
echo "Clean BED created: ${WGS_CLEAN}.bed/bim/fam (no * alleles, no duplicates)"

#===============================================================================
# PHASE 1d: Compute allele frequencies
#===============================================================================
echo "===== PHASE 1d: Allele frequencies ====="

plink2 --bfile ${WGS_CLEAN} \
    --freq counts \
    --threads ${THREADS} \
    --out wgs_6genes_freq

echo "Frequency file: wgs_6genes_freq.acount"

#===============================================================================
# SAFETY CHECK: Verify BIM A1/A2 orientation
# plink2 VCF→BED: A1 = ALT, A2 = REF. Do NOT use --ref-first in step 2.
#===============================================================================
echo "===== BIM ORIENTATION CHECK ====="
echo "First 5 variants (expect A1=ALT, A2=REF):"
head -5 ${WGS_CLEAN}.bim | awk '{print "  ID="$2, "A1="$5, "A2="$6}'
echo ""
echo "If A1 looks like REF alleles, STOP and investigate before continuing."
echo "For plink2 VCF conversions, A1=ALT is standard. Do NOT use --ref-first."

#===============================================================================
# PHASE 2: Regenie Step 1 — LOCO predictions from common array variants
#===============================================================================
echo "===== PHASE 2: Regenie Step 1 ====="

# QC filter array data for step 1
plink2 --bfile ${ARRAY_MERGED} \
    --maf 0.01 \
    --mac 100 \
    --geno 0.1 \
    --hwe 1e-15 \
    --mind 0.1 \
    --write-snplist \
    --write-samples --no-id-header \
    --threads ${THREADS} \
    --out qc_pass

echo "QC-passed variants: $(wc -l < qc_pass.snplist)"
echo "QC-passed samples:  $(wc -l < qc_pass.id)"

# Run step 1 — full in-memory mode (128GB available, no --lowmem needed)
regenie --step 1 \
    --bed ${ARRAY_MERGED} \
    --extract qc_pass.snplist \
    --keep qc_pass.id \
    --phenoFile ${PHENO_FILE} \
    --phenoCol ${PHENO_COL} \
    --covarFile ${COVAR_FILE} \
    --covarColList ${COVAR_COLS} \
    --catCovarList ${CAT_COVAR} \
    --bt \
    --bsize 1000 \
    --threads ${THREADS} \
    --out ${STEP1_OUT}

echo "Step 1 complete: ${STEP1_OUT}_pred.list"

#===============================================================================
# PHASE 3: Build annotation files (GENCODE GTF functional sub-regions)
#===============================================================================
echo "===== PHASE 3: Build annotation files ====="

python3 build_annotations.py \
    --bim ${WGS_CLEAN}.bim \
    --gtf ${GTF_FILE} \
    --out-prefix ${ANNOT_PREFIX}

echo "Annotation files created:"
echo "  ${ANNOT_PREFIX}.annot"
echo "  ${ANNOT_PREFIX}.setlist"
echo "  ${ANNOT_PREFIX}.masks"
echo ""
echo "Annotation summary:"
wc -l ${ANNOT_PREFIX}.annot ${ANNOT_PREFIX}.setlist ${ANNOT_PREFIX}.masks

#===============================================================================
# PHASE 4a: Regenie Step 2 — Primary run (AAF bins 0.001, 0.01)
# Broadest useful thresholds: captures rare + low-frequency variants.
# NO --ref-first (plink2 VCF→BED: A1=ALT, not REF)
#===============================================================================
echo "===== PHASE 4a: Regenie Step 2 — Primary (maxAAF 0.01) ====="

/opt/conda/envs/regenie_env/bin/regenie --step 2 \
    --bed ${WGS_CLEAN} \
    --phenoFile ${PHENO_FILE} \
    --phenoCol ${PHENO_COL} \
    --covarFile ${COVAR_FILE} \
    --covarColList ${COVAR_COLS} \
    --catCovarList ${CAT_COVAR} \
    --bt --firth --approx \
    --pred ${STEP1_OUT}_pred.list \
    --anno-file ${ANNOT_PREFIX}.annot \
    --set-list ${ANNOT_PREFIX}.setlist \
    --mask-def ${ANNOT_PREFIX}.masks \
    --aaf-bins 0.001,0.01 \
    --vc-tests skato,skat,acato,acatv,skato-acat \
    --vc-maxAAF 0.01 \
    --vc-MACthr 5 \
    --minMAC 1 \
    --bsize 400 \
    --check-burden-files --strict-check-burden \
    --write-mask-snplist \
    --threads ${THREADS} \
    --out ${STEP2_PREFIX}_primary

echo "Primary run complete: ${STEP2_PREFIX}_primary_scoliosis.regenie"

#===============================================================================
# PHASE 4b: Regenie Step 2 — Ultra-rare only (AAF bins 0.0001, 0.001)
# Focused on very rare variants only.
#===============================================================================
echo "===== PHASE 4b: Regenie Step 2 — Ultra-rare (maxAAF 0.001) ====="

/opt/conda/envs/regenie_env/bin/regenie --step 2 \
    --bed ${WGS_CLEAN} \
    --phenoFile ${PHENO_FILE} \
    --phenoCol ${PHENO_COL} \
    --covarFile ${COVAR_FILE} \
    --covarColList ${COVAR_COLS} \
    --catCovarList ${CAT_COVAR} \
    --bt --firth --approx \
    --pred ${STEP1_OUT}_pred.list \
    --anno-file ${ANNOT_PREFIX}.annot \
    --set-list ${ANNOT_PREFIX}.setlist \
    --mask-def ${ANNOT_PREFIX}.masks \
    --aaf-bins 0.0001,0.001 \
    --vc-tests skato,skat,acato,acatv,skato-acat \
    --vc-maxAAF 0.001 \
    --vc-MACthr 5 \
    --minMAC 1 \
    --bsize 400 \
    --check-burden-files --strict-check-burden \
    --threads ${THREADS} \
    --out ${STEP2_PREFIX}_ultrarare

echo "Ultra-rare run complete: ${STEP2_PREFIX}_ultrarare_scoliosis.regenie"

#===============================================================================
# PHASE 4c: Regenie Step 2 — Singleton-focused (AAF bins 0.00001, 0.0001)
# Most conservative: singletons and near-singletons only.
#===============================================================================
echo "===== PHASE 4c: Regenie Step 2 — Singletons (maxAAF 0.0001) ====="

/opt/conda/envs/regenie_env/bin/regenie --step 2 \
    --bed ${WGS_CLEAN} \
    --phenoFile ${PHENO_FILE} \
    --phenoCol ${PHENO_COL} \
    --covarFile ${COVAR_FILE} \
    --covarColList ${COVAR_COLS} \
    --catCovarList ${CAT_COVAR} \
    --bt --firth --approx \
    --pred ${STEP1_OUT}_pred.list \
    --anno-file ${ANNOT_PREFIX}.annot \
    --set-list ${ANNOT_PREFIX}.setlist \
    --mask-def ${ANNOT_PREFIX}.masks \
    --aaf-bins 0.00001,0.0001 \
    --vc-tests skato,skat,acato,acatv,skato-acat \
    --vc-maxAAF 0.0001 \
    --vc-MACthr 5 \
    --minMAC 1 \
    --bsize 400 \
    --check-burden-files --strict-check-burden \
    --threads ${THREADS} \
    --out ${STEP2_PREFIX}_singletons

echo "Singleton run complete: ${STEP2_PREFIX}_singletons_scoliosis.regenie"

#===============================================================================
# PHASE 5: Quick results summary
#===============================================================================
echo ""
echo "==============================================================================="
echo "===== PIPELINE COMPLETE ====="
echo "==============================================================================="
echo ""
echo "Results files:"
echo "  Primary (AAF≤0.01):     ${STEP2_PREFIX}_primary_scoliosis.regenie"
echo "  Ultra-rare (AAF≤0.001): ${STEP2_PREFIX}_ultrarare_scoliosis.regenie"
echo "  Singletons (AAF≤0.0001):${STEP2_PREFIX}_singletons_scoliosis.regenie"
echo ""
echo "Gene sets tested (functional sub-regions per gene):"
echo "  COL6A3, MLPH, PRLH, RAB17, LRRFIP1_unique, RBM44_unique"
echo "  LRRFIP1_RBM44_overlap"
echo "  LRRFIP1_whole, RBM44_whole (non-independent, report separately)"
echo ""
echo "Masks tested:"
echo "  M1_cds          = CDS only"
echo "  M2_coding       = CDS + UTR5 + UTR3"
echo "  M3_exonic       = CDS + UTR5 + UTR3 (alias)"
echo "  M4_regulatory   = UTR5 + UTR3 + promoter"
echo "  M5_intronic     = intronic only"
echo "  M6_genic        = gene body (no promoter)"
echo "  M7_all          = everything"
echo ""
echo "Tests: SKAT-O, SKAT, Burden, ACAT-V"
echo ""
echo "Quick peek at top hits (primary run, sorted by LOG10P):"
awk 'NR>1 && $12 != "NA" {print $3, $8, $12}' \
    ${STEP2_PREFIX}_primary_scoliosis.regenie \
    | sort -k3 -rn | head -20
