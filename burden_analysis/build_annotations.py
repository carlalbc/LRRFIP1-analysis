#!/usr/bin/env python3
"""
Build regenie annotation, set-list, and mask files for scoliosis SKAT-O analysis.
Uses GENCODE GTF-based functional sub-region annotations.

Strategy:
  - Uses MANE Select / Ensembl canonical transcript per gene from GENCODE v49
  - All 6 genes get functional annotation (CDS, UTR5, UTR3, intronic, promoter)
  - LRRFIP1 and RBM44 base sets EXCLUDE the overlap zone (237,798,389-237,813,682)
    so each base set contains only variants unique to that gene
  - Overlap zone tested as a separate set (LRRFIP1_RBM44_overlap)
  - Whole-gene aggregates (LRRFIP1_whole, RBM44_whole) include overlap for sensitivity
  - Multiple masks test different functional groupings

Gene sets tested:
  - COL6A3, MLPH, PRLH, RAB17: full gene + promoter
  - LRRFIP1: gene body EXCLUDING overlap zone (unique variants only)
  - RBM44: gene body EXCLUDING overlap zone (unique variants only)
  - LRRFIP1_RBM44_overlap: overlap zone only (15.3 kb)
  - LRRFIP1_whole: full LRRFIP1 gene INCLUDING overlap (sensitivity analysis)
  - RBM44_whole: full RBM44 gene INCLUDING overlap (sensitivity analysis)

Annotation categories:
  CDS       = coding sequence (from canonical transcript CDS features)
  UTR5      = 5' UTR
  UTR3      = 3' UTR
  intronic  = within gene body but not exonic (using canonical transcript exons)
  promoter  = 1kb upstream of gene start (strand-aware)

Masks (applied to every gene; empty combos auto-skipped by regenie):
  M1_cds          = CDS only
  M2_coding       = CDS + UTR5 + UTR3
  M3_exonic       = CDS + UTR5 + UTR3  (same as M2, alias for clarity)
  M4_regulatory   = UTR5 + UTR3 + promoter
  M5_intronic     = intronic only
  M6_genic        = CDS + UTR5 + UTR3 + intronic (gene body, no promoter)
  M7_all          = CDS + UTR5 + UTR3 + intronic + promoter

Usage:
    python3 build_annotations.py \\
        --bim wgs_6genes_clean.bim \\
        --gtf gencode.v49.annotation.gtf \\
        --out-prefix regenie_annotations
"""

import argparse
import re
import sys
from collections import OrderedDict


# =============================================================================
# Gene definitions — Ensembl IDs and genomic boundaries (GRCh38, GENCODE v49)
# =============================================================================

GENE_DEFS = OrderedDict([
    ("COL6A3",  {"ensembl": "ENSG00000163359", "gene_start": 237_324_003, "gene_end": 237_414_405, "strand": "-"}),
    ("MLPH",    {"ensembl": "ENSG00000115648", "gene_start": 237_485_428, "gene_end": 237_555_322, "strand": "+"}),
    ("PRLH",    {"ensembl": "ENSG00000071677", "gene_start": 237_559_555, "gene_end": 237_567_288, "strand": "+"}),
    ("RAB17",   {"ensembl": "ENSG00000124839", "gene_start": 237_574_322, "gene_end": 237_601_623, "strand": "-"}),
    ("LRRFIP1", {"ensembl": "ENSG00000124831", "gene_start": 237_627_554, "gene_end": 237_813_682, "strand": "+"}),
    ("RBM44",   {"ensembl": "ENSG00000177483", "gene_start": 237_798_389, "gene_end": 237_842_808, "strand": "+"}),
])

# Overlap zone between LRRFIP1 and RBM44
OVERLAP_START = 237_798_389
OVERLAP_END   = 237_813_682

# Promoter size (bp upstream of gene start, strand-aware)
PROMOTER_SIZE = 1000

# Annotation categories
CATEGORIES = ["CDS", "UTR5", "UTR3", "intronic", "promoter"]


# =============================================================================
# GTF parsing
# =============================================================================

def parse_gtf_for_genes(gtf_path, gene_defs):
    """
    Parse GENCODE GTF and extract MANE Select / canonical transcript features
    for each gene.

    Returns dict: gene_name -> {
        "cds": [(start, end), ...],
        "utr5": [(start, end), ...],
        "utr3": [(start, end), ...],
        "exons": [(start, end), ...],
        "transcript_id": str,
    }
    """
    # Build lookup: ensembl_id (without version) -> gene_name
    eid_to_gene = {}
    for gname, gdef in gene_defs.items():
        eid_to_gene[gdef["ensembl"]] = gname

    # First pass: collect all features for our genes
    raw_features = {gname: [] for gname in gene_defs}

    print("Parsing GTF file...")
    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#") or not line.startswith("chr2\t"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue

            attrs = fields[8]
            attr_dict = {}
            for m in re.finditer(r'(\w+) "([^"]+)"', attrs):
                attr_dict[m.group(1)] = m.group(2)

            gid = attr_dict.get("gene_id", "").split(".")[0]
            if gid not in eid_to_gene:
                continue

            gname = eid_to_gene[gid]
            raw_features[gname].append({
                "feature": fields[2],
                "start": int(fields[3]),
                "end": int(fields[4]),
                "strand": fields[6],
                "transcript_id": attr_dict.get("transcript_id", ""),
                "tag": attrs,
            })

    # Second pass: find canonical transcript and extract features per gene
    result = {}
    for gname, gdef in gene_defs.items():
        feats = raw_features[gname]

        # Find MANE Select, fall back to Ensembl_canonical
        canonical_tid = None
        mane_tid = None
        for f in feats:
            if f["feature"] == "transcript":
                if "MANE_Select" in f["tag"]:
                    mane_tid = f["transcript_id"]
                if "Ensembl_canonical" in f["tag"]:
                    canonical_tid = f["transcript_id"]

        use_tid = mane_tid or canonical_tid
        if not use_tid:
            print(f"  WARNING: No canonical transcript found for {gname}")
            result[gname] = {"cds": [], "utr5": [], "utr3": [], "exons": [], "transcript_id": None}
            continue

        print(f"  {gname}: using {use_tid}" + (" (MANE Select)" if use_tid == mane_tid else " (Ensembl canonical)"))

        t_feats = [f for f in feats if f["transcript_id"] == use_tid]
        cds = sorted([(f["start"], f["end"]) for f in t_feats if f["feature"] == "CDS"])
        exons = sorted([(f["start"], f["end"]) for f in t_feats if f["feature"] == "exon"])
        utrs = sorted([(f["start"], f["end"]) for f in t_feats if f["feature"] == "UTR"])

        # Classify UTRs as 5' or 3' based on strand and position relative to CDS
        strand = gdef["strand"]
        utr5, utr3 = [], []
        if cds:
            cds_min = min(s for s, e in cds)
            cds_max = max(e for s, e in cds)
            for s, e in utrs:
                if strand == "+":
                    if e < cds_min:
                        utr5.append((s, e))
                    else:
                        utr3.append((s, e))
                else:  # minus strand
                    if s > cds_max:
                        utr5.append((s, e))
                    else:
                        utr3.append((s, e))

        result[gname] = {
            "cds": cds,
            "utr5": sorted(utr5),
            "utr3": sorted(utr3),
            "exons": exons,
            "transcript_id": use_tid,
        }

        # Print summary
        print(f"           CDS: {len(cds)} intervals ({sum(e-s+1 for s,e in cds)} bp)")
        print(f"           5'UTR: {len(utr5)} intervals ({sum(e-s+1 for s,e in utr5)} bp)")
        print(f"           3'UTR: {len(utr3)} intervals ({sum(e-s+1 for s,e in utr3)} bp)")
        print(f"           Exons: {len(exons)} intervals ({sum(e-s+1 for s,e in exons)} bp)")

    return result


# =============================================================================
# Interval helpers
# =============================================================================

def merge_intervals(intervals):
    """Merge overlapping/adjacent intervals."""
    if not intervals:
        return []
    merged = []
    for s, e in sorted(intervals):
        if merged and s <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
        else:
            merged.append((s, e))
    return merged


def subtract_intervals(base, remove):
    """Subtract remove intervals from base intervals."""
    result = []
    for bs, be in base:
        pieces = [(bs, be)]
        for rs, re_ in remove:
            new_pieces = []
            for ps, pe in pieces:
                if rs > pe or re_ < ps:
                    new_pieces.append((ps, pe))
                else:
                    if ps < rs:
                        new_pieces.append((ps, rs - 1))
                    if pe > re_:
                        new_pieces.append((re_ + 1, pe))
            pieces = new_pieces
        result.extend(pieces)
    return sorted(result)


def point_in_intervals(pos, intervals):
    """Check if a position falls in any interval."""
    for s, e in intervals:
        if s <= pos <= e:
            return True
        if pos < s:
            break  # intervals are sorted
    return False


# =============================================================================
# Build functional regions per gene
# =============================================================================

def build_gene_regions(gene_defs, gtf_data):
    """
    Build functional sub-region interval lists for each gene.

    Returns dict: gene_name -> {
        "CDS": [(start, end), ...],
        "UTR5": [...],
        "UTR3": [...],
        "intronic": [...],
        "promoter": [...],
        "full_region": (start, end),  # gene body + promoter
    }
    """
    regions = {}

    for gname, gdef in gene_defs.items():
        gstart = gdef["gene_start"]
        gend = gdef["gene_end"]
        strand = gdef["strand"]
        gdata = gtf_data[gname]

        # Promoter: 1kb upstream (strand-aware)
        if strand == "+":
            prom_start = gstart - PROMOTER_SIZE
            prom_end = gstart - 1
        else:
            prom_start = gend + 1
            prom_end = gend + PROMOTER_SIZE

        promoter = [(prom_start, prom_end)]

        # CDS, UTR from canonical transcript
        cds = merge_intervals(gdata["cds"])
        utr5 = merge_intervals(gdata["utr5"])
        utr3 = merge_intervals(gdata["utr3"])

        # All exonic = union of CDS + UTR (they come from exon features)
        all_exonic = merge_intervals(cds + utr5 + utr3)

        # Intronic = gene body minus all exonic
        gene_body = [(gstart, gend)]
        intronic = subtract_intervals(gene_body, all_exonic)

        # Full analysis region = promoter + gene body
        full_start = min(prom_start, gstart)
        full_end = max(prom_end, gend)

        regions[gname] = {
            "CDS": cds,
            "UTR5": utr5,
            "UTR3": utr3,
            "intronic": intronic,
            "promoter": promoter,
            "full_region": (full_start, full_end),
        }

    return regions


# =============================================================================
# Variant assignment
# =============================================================================

def parse_bim(bim_path):
    """Read plink .bim file, return list of (chr, var_id, pos) tuples."""
    variants = []
    with open(bim_path) as f:
        for line in f:
            fields = line.strip().split("\t")
            chrom = fields[0]
            var_id = fields[1]
            pos = int(fields[3])
            variants.append((chrom, var_id, pos))
    return variants


def _get_category(pos, region):
    """Determine functional category for a position in a gene's regions."""
    if point_in_intervals(pos, region["CDS"]):
        return "CDS"
    elif point_in_intervals(pos, region["UTR5"]):
        return "UTR5"
    elif point_in_intervals(pos, region["UTR3"]):
        return "UTR3"
    elif point_in_intervals(pos, region["intronic"]):
        return "intronic"
    elif point_in_intervals(pos, region["promoter"]):
        return "promoter"
    return "intronic"  # fallback


def assign_variants(variants, gene_defs, regions):
    """
    Assign each chr2 variant to gene sets with functional annotations.

    For each gene, a variant within the gene's full region (gene body + promoter)
    gets assigned the most specific functional category:
        CDS > UTR5 > UTR3 > intronic > promoter

    A variant can appear in multiple genes (e.g., LRRFIP1 overlap zone variants
    also appear in the overlap set).

    Returns:
        annotations: list of (var_id, gene_name, category) for annotation file
        set_variants: dict of gene_name -> list of var_ids for set-list
    """
    # Build ordered gene list: main genes + overlap + whole-gene sets
    gene_names = list(gene_defs.keys())

    annotations = []
    set_variants = OrderedDict((gn, []) for gn in gene_names)

    # Also track overlap and whole-gene sets
    set_variants["LRRFIP1_RBM44_overlap"] = []
    set_variants["LRRFIP1_whole"] = []
    set_variants["RBM44_whole"] = []

    for chrom, var_id, pos in variants:
        if chrom not in ("2", "chr2"):
            continue

        # Check each gene
        in_overlap = OVERLAP_START <= pos <= OVERLAP_END

        for gname in gene_names:
            r = regions[gname]
            full_start, full_end = r["full_region"]

            if pos < full_start or pos > full_end:
                continue

            # Determine functional category (priority order)
            cat = _get_category(pos, r)

            # For LRRFIP1 and RBM44 base sets: EXCLUDE overlap-zone variants
            # so that base sets are truly "unique" to each gene.
            # Overlap variants go into the dedicated overlap set instead.
            if in_overlap and gname in ("LRRFIP1", "RBM44"):
                continue

            annotations.append((var_id, gname, cat))
            set_variants[gname].append(var_id)

        # Overlap set: variants in the LRRFIP1/RBM44 overlap zone
        if OVERLAP_START <= pos <= OVERLAP_END:
            # These get a generic annotation — masks will handle them
            annotations.append((var_id, "LRRFIP1_RBM44_overlap", "CDS"
                                if point_in_intervals(pos, regions["LRRFIP1"]["CDS"]) or
                                   point_in_intervals(pos, regions["RBM44"]["CDS"])
                                else "intronic"))
            set_variants["LRRFIP1_RBM44_overlap"].append(var_id)

        # Whole-gene aggregates — need their own annotation entries
        # (regenie matches gene_name in annotation to gene_name in set-list)
        lrrfip1_full = regions["LRRFIP1"]["full_region"]
        if lrrfip1_full[0] <= pos <= lrrfip1_full[1]:
            set_variants["LRRFIP1_whole"].append(var_id)
            # Reuse the category already determined for LRRFIP1
            lcat = _get_category(pos, regions["LRRFIP1"])
            annotations.append((var_id, "LRRFIP1_whole", lcat))
        rbm44_full = regions["RBM44"]["full_region"]
        if rbm44_full[0] <= pos <= rbm44_full[1]:
            set_variants["RBM44_whole"].append(var_id)
            rcat = _get_category(pos, regions["RBM44"])
            annotations.append((var_id, "RBM44_whole", rcat))

    return annotations, set_variants


# =============================================================================
# Output files
# =============================================================================

def write_annotation_file(annotations, out_path):
    """
    Write regenie annotation file.
    Format: variant_id <tab> gene_name <tab> annotation_category
    """
    with open(out_path, "w") as f:
        for var_id, gene_name, cat in annotations:
            f.write(f"{var_id}\t{gene_name}\t{cat}\n")


def write_set_list(set_variants, regions, out_path):
    """
    Write regenie set-list file.
    Format: gene_name <tab> chrom <tab> pos <tab> comma-separated variant_ids
    """
    with open(out_path, "w") as f:
        for set_name, var_ids in set_variants.items():
            if not var_ids:
                print(f"  WARNING: Set '{set_name}' has 0 variants — skipping.")
                continue

            # Get start position for set
            if set_name in regions:
                pos = regions[set_name]["full_region"][0]
            elif set_name == "LRRFIP1_RBM44_overlap":
                pos = OVERLAP_START
            elif set_name == "LRRFIP1_whole":
                pos = regions["LRRFIP1"]["full_region"][0]
            elif set_name == "RBM44_whole":
                pos = regions["RBM44"]["full_region"][0]
            else:
                pos = 0

            f.write(f"{set_name}\t2\t{pos}\t{','.join(var_ids)}\n")


def write_mask_file(out_path):
    """
    Write regenie mask definition file.
    Multiple masks for different functional groupings.

    Regenie applies every mask to every gene. If a mask's categories
    don't match any variants in a gene, regenie silently skips it.
    """
    masks = [
        ("M1_cds",        "CDS"),
        ("M2_coding",     "CDS,UTR5,UTR3"),
        ("M3_regulatory", "UTR5,UTR3,promoter"),
        ("M4_intronic",   "intronic"),
        ("M5_genic",      "CDS,UTR5,UTR3,intronic"),
        ("M6_all",        "CDS,UTR5,UTR3,intronic,promoter"),
    ]
    with open(out_path, "w") as f:
        for mask_name, categories in masks:
            f.write(f"{mask_name}\t{categories}\n")
    return masks


# =============================================================================
# Summary
# =============================================================================

def print_summary(annotations, set_variants, gene_defs, regions):
    """Print summary of variant assignments."""
    from collections import Counter

    print("\n" + "=" * 70)
    print("VARIANT ASSIGNMENT SUMMARY")
    print("=" * 70)

    for gname in gene_defs:
        n = len(set_variants[gname])
        r = regions[gname]
        full_s, full_e = r["full_region"]

        # Count by category
        cats = Counter(cat for vid, gn, cat in annotations if gn == gname)
        print(f"\n  {gname} ({full_s:,d}-{full_e:,d}): {n:,d} variants")
        for cat in CATEGORIES:
            if cats[cat] > 0:
                print(f"    {cat:12s}: {cats[cat]:>8,d}")

    # Overlap and whole-gene sets
    print(f"\n  --- Additional sets ---")
    for sn in ["LRRFIP1_RBM44_overlap", "LRRFIP1_whole", "RBM44_whole"]:
        n = len(set_variants[sn])
        print(f"  {sn:30s}: {n:>8,d} variants")

    # Validation for LRRFIP1
    lrrfip1_n = len(set_variants["LRRFIP1"])
    rbm44_n = len(set_variants["RBM44"])
    overlap_n = len(set_variants["LRRFIP1_RBM44_overlap"])
    whole_l = len(set_variants["LRRFIP1_whole"])
    whole_r = len(set_variants["RBM44_whole"])

    print(f"\n  --- Validation ---")
    print(f"  LRRFIP1 (unique, excl overlap) = {lrrfip1_n:,d}")
    print(f"  RBM44 (unique, excl overlap)   = {rbm44_n:,d}")
    print(f"  Overlap set                    = {overlap_n:,d}")
    print(f"  LRRFIP1_whole                  = {whole_l:,d} (should = LRRFIP1 + overlap = {lrrfip1_n + overlap_n:,d}: {'OK' if whole_l == lrrfip1_n + overlap_n else 'MISMATCH'})")
    print(f"  RBM44_whole                    = {whole_r:,d} (should = RBM44 + overlap = {rbm44_n + overlap_n:,d}: {'OK' if whole_r == rbm44_n + overlap_n else 'MISMATCH'})")
    print("=" * 70)


# =============================================================================
# Main
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Build regenie annotation files with GENCODE functional sub-regions"
    )
    parser.add_argument("--bim", required=True,
                        help="Path to plink .bim file (wgs_6genes_clean.bim)")
    parser.add_argument("--gtf", required=True,
                        help="Path to GENCODE GTF file (gencode.v49.annotation.gtf)")
    parser.add_argument("--out-prefix", required=True,
                        help="Output prefix (creates .annot, .setlist, .masks)")
    args = parser.parse_args()

    # Step 1: Parse GTF for canonical transcripts
    gtf_data = parse_gtf_for_genes(args.gtf, GENE_DEFS)

    # Step 2: Build functional regions per gene
    print("\nBuilding functional sub-regions...")
    regions = build_gene_regions(GENE_DEFS, gtf_data)

    # Print region summary
    for gname in GENE_DEFS:
        r = regions[gname]
        total_bp = sum(
            sum(e - s + 1 for s, e in r[cat])
            for cat in CATEGORIES
        )
        gene_bp = GENE_DEFS[gname]["gene_end"] - GENE_DEFS[gname]["gene_start"] + 1
        print(f"  {gname}: {total_bp:,d} bp covered (gene body {gene_bp:,d} + promoter {PROMOTER_SIZE})")

    # Step 3: Parse BIM
    print(f"\nReading BIM file: {args.bim}")
    variants = parse_bim(args.bim)
    chr2_count = sum(1 for c, _, _ in variants if c in ("2", "chr2"))
    print(f"  Total variants: {len(variants):,d}  (chr2: {chr2_count:,d})")

    # Step 4: Assign variants
    print("\nAssigning variants to gene sets with functional annotations...")
    annotations, set_variants = assign_variants(variants, GENE_DEFS, regions)

    print_summary(annotations, set_variants, GENE_DEFS, regions)

    # Step 5: Write output files
    annot_path = f"{args.out_prefix}.annot"
    setlist_path = f"{args.out_prefix}.setlist"
    mask_path = f"{args.out_prefix}.masks"

    print(f"\nWriting annotation file: {annot_path}")
    write_annotation_file(annotations, annot_path)

    print(f"Writing set-list file:   {setlist_path}")
    write_set_list(set_variants, regions, setlist_path)

    print(f"Writing mask file:       {mask_path}")
    masks = write_mask_file(mask_path)
    for mname, cats in masks:
        print(f"  {mname:20s} = {cats}")

    print("\nDone. Files ready for regenie --step 2.")
    print("\nRecommended regenie step 2 command:")
    print(f"""
regenie --step 2 \\
    --bed wgs_6genes_clean \\
    --ref-first \\
    --phenoFile pheno_500k.regenie.txt \\
    --phenoCol scoliosis \\
    --covarFile covariates_p22009_500k_imputed.regenie.txt \\
    --covarColList p21022,p22009_a1,p22009_a2,p22009_a3,p22009_a4,p22009_a5 \\
    --catCovarList p31 \\
    --bt --firth --approx \\
    --pred step1_scoliosis_pred.list \\
    --anno-file {annot_path} \\
    --set-list {setlist_path} \\
    --mask-def {mask_path} \\
    --aaf-bins 0.001,0.01 \\
    --vc-tests skato,skat,acatv,acato \\
    --vc-maxAAF 0.01 \\
    --vc-MACthr 5 \\
    --minMAC 1 \\
    --check-burden-files --strict-check-burden \\
    --write-mask-snplist \\
    --bsize 400 \\
    --threads 16 \\
    --out scoliosis_skato
""")


if __name__ == "__main__":
    main()
