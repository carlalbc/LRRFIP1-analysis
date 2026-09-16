#!/bin/zsh
set -euo pipefail

# =============================================================
# IsoQuant LRRFIP1: FLNC reads (PacBio recommended workflow)
#
# This follows the PacBio-recommended pipeline:
#   FLNC reads → pbmm2 align → IsoQuant (skip isoseq cluster)
#
# BEFORE RUNNING:
#   1. Find and copy per-sample FLNC BAMs into $FLNC_DIR below
#      Look in SMRT Link output for files like:
#        - refine/BioSample_1.flnc.bam
#        - or refine.flnc.bam (if single sample)
#      They should be the output of `isoseq refine`
#
#   2. Update the FLNC_FILES array below with actual filenames
#
#   3. Install pbmm2 if not already available:
#        conda install -c bioconda pbmm2
#      (or use minimap2 — set USE_MINIMAP2=true below)
#
# Strategy: align FLNC to chr2-only reference (saves RAM),
#   extract LRRFIP1 region, run IsoQuant
# =============================================================

# ── Configuration ──
BASE=./data/pacbio
WORKDIR=$BASE/isoquant_lrrfip1
REF=$BASE/pigeon_v49_rerun/inputs/human_GRCh38_no_alt_analysis_set.fasta
GTF=$BASE/pigeon_v49_rerun/inputs/gencode.v49.annotation.gtf
CHR2_REF=$WORKDIR/chr2.fasta
CHR2_GTF=$WORKDIR/chr2_lrrfip1.gtf
REGION="chr2:237600000-237850000"
THREADS=4

# Set to true if pbmm2 is not available
USE_MINIMAP2=false

# ── FLNC input files ──
# UPDATE THESE with your actual FLNC file paths!
FLNC_DIR=$BASE/isoquant_lrrfip1/flnc_inputs
LABELS=(Son_WT Mother_WT Son_aff Father_aff)

# Option A: Per-sample FLNC BAMs (from isoseq refine, one per BioSample)
# Uncomment and update these paths:
# FLNC_FILES=(
#     "$FLNC_DIR/BioSample_1.flnc.bam"
#     "$FLNC_DIR/BioSample_2.flnc.bam"
#     "$FLNC_DIR/BioSample_3.flnc.bam"
#     "$FLNC_DIR/BioSample_4.flnc.bam"
# )
# INPUT_TYPE="bam"

# Option B: Per-sample FLNC FASTAs (converted from BAMs)
# FLNC_FILES=(
#     "$FLNC_DIR/BioSample_1.flnc.fasta"
#     "$FLNC_DIR/BioSample_2.flnc.fasta"
#     "$FLNC_DIR/BioSample_3.flnc.fasta"
#     "$FLNC_DIR/BioSample_4.flnc.fasta"
# )
# INPUT_TYPE="fasta"

# Option C: Single pooled FLNC BAM (will need to be split by read group / barcode)
# FLNC_POOLED="$FLNC_DIR/pooled.flnc.bam"
# INPUT_TYPE="pooled_bam"

# ── For now, exit if not configured ──
if [ -z "${FLNC_FILES+x}" ] && [ -z "${FLNC_POOLED+x}" ]; then
    echo "ERROR: No FLNC input files configured!"
    echo ""
    echo "Edit this script and uncomment one of the FLNC_FILES options above."
    echo ""
    echo "To find your FLNC files, look in SMRT Link for 'isoseq refine' output:"
    echo "  - Per-sample: BioSample_N.flnc.bam"
    echo "  - Pooled: something.flnc.bam"
    echo ""
    echo "If you only have BAMs and need FASTAs:"
    echo "  samtools fasta input.flnc.bam > output.flnc.fasta"
    echo ""
    echo "If you have a pooled BAM that needs splitting by sample:"
    echo "  samtools split -f '%!.flnc.bam' pooled.flnc.bam"
    exit 1
fi

# ── Prep chr2 reference (reuse if exists) ──
if [ ! -f "$CHR2_REF" ]; then
    echo "Extracting chr2 from reference..."
    samtools faidx "$REF" chr2 > "$CHR2_REF"
    samtools faidx "$CHR2_REF"
fi

if [ ! -f "$CHR2_GTF" ]; then
    echo "Extracting LRRFIP1 region GTF..."
    grep "^chr2" "$GTF" | awk -F'\t' '$4 >= 237500000 && $5 <= 237900000' > "$CHR2_GTF"
fi

mkdir -p $WORKDIR/flnc_aligned $WORKDIR/flnc_lrrfip1_reads

# Keep computer awake
caffeinate -s &
CAFF_PID=$!
trap "kill $CAFF_PID 2>/dev/null" EXIT

echo ""
echo "============================================"
echo " Step 1: Align FLNC reads to chr2"
echo "============================================"

for i in 1 2 3 4; do
    LABEL=${LABELS[$i]}
    FLNC=${FLNC_FILES[$i]}
    BAM="$WORKDIR/flnc_aligned/${LABEL}.bam"
    LRRFIP1_FASTA="$WORKDIR/flnc_lrrfip1_reads/${LABEL}_lrrfip1.fasta"

    if [ -f "$LRRFIP1_FASTA" ] && [ -s "$LRRFIP1_FASTA" ]; then
        echo "  $LABEL: LRRFIP1 FASTA already exists, skipping."
        continue
    fi

    echo ""
    echo "  ── $LABEL ──"
    echo "  Input: $(basename "$FLNC")"

    if [ "$USE_MINIMAP2" = true ]; then
        echo "  Aligning with minimap2..."
        minimap2 -ax splice:hq -t $THREADS --secondary=no "$CHR2_REF" "$FLNC" 2>/dev/null \
            | samtools sort -@ 2 -o "$BAM"
    else
        echo "  Aligning with pbmm2..."
        pbmm2 align "$CHR2_REF" "$FLNC" "$BAM" \
            --preset ISOSEQ --sort -j $THREADS --log-level WARN
    fi
    samtools index "$BAM"

    # Extract LRRFIP1 region reads as FASTA with UNIQUE IDs
    N_REGION=$(samtools view -c "$BAM" $REGION)
    echo "  Found $N_REGION reads in LRRFIP1 region"

    # Use proper read names (FLNC IDs are already globally unique — movie/zmw/ccs format)
    # But prefix with label just to be safe
    samtools view "$BAM" $REGION \
        | awk -v label="$LABEL" '{print ">"label"__"$1"\n"$10}' \
        | awk '/^>/ { if (!(seen[$0]++)) { id=$0; getline seq; print id"\n"seq } }' \
        > "$LRRFIP1_FASTA"

    N_FASTA=$(grep -c "^>" "$LRRFIP1_FASTA" || echo 0)
    echo "  Wrote $N_FASTA unique reads to $(basename $LRRFIP1_FASTA)"

    # Clean up large BAM
    rm -f "$BAM" "${BAM}.bai"
    echo "  Cleaned up aligned BAM"
done

echo ""
echo "Read counts per sample:"
TOTAL=0
for f in $WORKDIR/flnc_lrrfip1_reads/*_lrrfip1.fasta; do
    n=$(grep -c "^>" "$f" || echo 0)
    TOTAL=$((TOTAL + n))
    echo "  $(basename $f _lrrfip1.fasta): $n"
done
echo "  Total: $TOTAL"

# Verify no ID collisions
N_COLLISIONS=$(cat $WORKDIR/flnc_lrrfip1_reads/*_lrrfip1.fasta | grep "^>" | sort | uniq -d | wc -l)
if [ "$N_COLLISIONS" -gt 0 ]; then
    echo "WARNING: $N_COLLISIONS ID collisions detected!"
else
    echo "  No ID collisions — good."
fi

echo ""
echo "============================================"
echo " Step 2: Run IsoQuant"
echo "============================================"

isoquant \
    --reference "$CHR2_REF" \
    --genedb "$CHR2_GTF" \
    --fastq \
        "$WORKDIR/flnc_lrrfip1_reads/Son_WT_lrrfip1.fasta" \
        "$WORKDIR/flnc_lrrfip1_reads/Mother_WT_lrrfip1.fasta" \
        "$WORKDIR/flnc_lrrfip1_reads/Son_aff_lrrfip1.fasta" \
        "$WORKDIR/flnc_lrrfip1_reads/Father_aff_lrrfip1.fasta" \
    --labels "Son_WT" "Mother_WT" "Son_aff" "Father_aff" \
    --data_type pacbio_ccs \
    --fl_data \
    --polya_trimmed all \
    --sqanti_output \
    --check_canonical \
    --threads $THREADS \
    --force \
    --output "$WORKDIR/isoquant_output_flnc"

echo ""
echo "============================================"
echo " DONE!"
echo "============================================"
echo ""
echo "Gene counts (FLNC-based):"
cat "$WORKDIR/isoquant_output_flnc/OUT/OUT.gene_grouped_file_name_counts.tsv"
echo ""
echo "Transcript counts:"
cat "$WORKDIR/isoquant_output_flnc/OUT/OUT.transcript_grouped_file_name_counts.tsv"
echo ""
echo "For reference — HQ transcript run (v3):"
echo "  Father_aff=102, Mother_WT=117, Son_WT=108, Son_aff=104"
echo ""
echo "Pigeon FSM totals:"
echo "  Son_WT=166, Mother_WT=165, Son_aff=143, Father_aff=202"
