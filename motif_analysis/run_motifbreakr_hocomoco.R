## --- motifbreakR analysis (hg38) — per-SNP CSVs + robust PDFs ---
## Transcription factor binding site disruption analysis using HOCOMOCO motifs.
## Input TSV must have columns: chr  pos  ref  alt  (1-based hg38)

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(Biostrings)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(motifbreakR)
  library(BiocParallel)
  library(grid)        # fallback plotting
  library(gridExtra)   # fallback table
  library(seqLogo)     # logos
})

## ---- EDIT THESE ----
snps_file   <- "./data/motifbreakr/snps_hg38.tsv"
out_dir     <- "./data/motifbreakr/output"
thr         <- 0.85          # score threshold (use when use_filterp = FALSE)
use_filterp <- TRUE          # TRUE to use p-value filtering instead of score
p_thresh    <- 1e-4          # only used if use_filterp = TRUE
top_k       <- 8             # max motifs per PDF page
## --------------------

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## 1) Read TSV and normalize
dt <- fread(snps_file)
setnames(dt, tolower(names(dt)))
stopifnot(all(c("chr","pos","ref","alt") %in% names(dt)))

dt[, `:=`(
  chr = ifelse(grepl("^chr", chr), chr, paste0("chr", chr)),
  pos = as.integer(pos),
  ref = toupper(ref),
  alt = toupper(alt)
)]
dt <- dt[, .(alt = unlist(strsplit(alt, ","))), by = .(chr, pos, ref)]

## 2) Quick REF vs hg38 sanity check (vectorized)
gr_chk  <- GRanges(seqnames = dt$chr, ranges = IRanges(dt$pos, width = 1))
ref_gen <- as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38, gr_chk))
ref_ok  <- data.table(ok = toupper(ref_gen) == dt$ref, ref_gen = toupper(ref_gen))

if (any(!ref_ok$ok)) {
  bad <- cbind(dt[!ref_ok$ok], ref_genome = ref_ok$ref_gen[!ref_ok$ok])
  print(bad)
  stop("Some REF alleles don't match hg38. Fix these first (off-by-one or wrong build).")
}

## 3) Write proper BED (0-based start, half-open end)
bed <- dt[, .(
  chrom  = chr,
  start  = pos - 1L,
  end    = pos,
  name   = sprintf("%s:%d:%s:%s", chr, pos, ref, alt),
  score  = 0L,
  strand = "+"
)]
bed_path <- file.path(out_dir, "tmp_snps.bed")
fwrite(bed, bed_path, sep = "\t", col.names = FALSE)

## 4) Build variants from BED (motifbreakR helper)
variants <- motifbreakR::snps.from.file(
  file = bed_path,
  search.genome = BSgenome.Hsapiens.UCSC.hg38,
  format = "bed"
)
if (length(variants) == 0L) stop("No variants parsed from BED.")

## 5) Motifs (HOCOMOCO bundle preferred)
motifs <- tryCatch({
  data(hocomoco, package = "motifbreakR")
  hocomoco
}, error = function(e) {
  data(motifbreakR_motif, package = "motifbreakR")
  motifbreakR_motif
})

## Fast lookup for fallback (metadata about motifs)
mot_meta <- as.data.frame(mcols(motifs))
mot_meta$idx <- seq_along(motifs)

## 6) Run motifbreakR
param <- BiocParallel::SerialParam()  # safest inside RStudio
message(sprintf("Running motifbreakR on %d SNP(s) with %d motif(s)...",
                length(variants), length(motifs)))

res <- motifbreakR::motifbreakR(
  snpList   = variants,
  pwmList   = motifs,
  threshold = if (use_filterp) p_thresh else thr,
  method    = "ic",
  filterp   = use_filterp,
  BPPARAM   = param
)

## 7) Helper to make a clean data.table (fix list-cols)
to_dt <- function(grx) {
  if (length(grx) == 0L) return(data.table())
  snpvec <- names(grx); names(grx) <- NULL
  df <- as.data.frame(grx); data.table::setDT(df)
  data.table::setnames(df, tolower(names(df)))
  df[, snp := snpvec]
  if (!"pctref"    %in% names(df) && "scoreref" %in% names(df)) df[, pctref := scoreref]
  if (!"pctalt"    %in% names(df) && "scorealt" %in% names(df)) df[, pctalt := scorealt]
  if (!"delta_pct" %in% names(df)) df[, delta_pct := pctalt - pctref]

  # unwrap common list-cols into scalars
  unlist1 <- function(x) if (is.list(x)) vapply(x, function(y) if (length(y)) y[[1]] else NA, FUN.VALUE=character(1)) else x
  for (nm in intersect(c("genesymbol","providername","providerid","seqmatch"), names(df))) {
    df[[nm]] <- unlist1(df[[nm]])
  }
  for (nm in intersect(c("motifpos","pctref","pctalt","delta_pct"), names(df))) {
    if (is.list(df[[nm]])) df[[nm]] <- vapply(df[[nm]], function(y) if (length(y)) as.numeric(y[[1]]) else NA_real_, FUN.VALUE = 0.0)
  }
  df[]
}

## 8) Fallback logo drawer (robust)
draw_fallback_logos <- function(tab_top, sid, page_txt) {
  n <- nrow(tab_top)
  if (n == 0L) {
    grid.newpage()
    grid.text(sprintf("No motif hits for %s", sid),
              gp = gpar(fontsize = 12, fontface = "bold"))
    return(invisible())
  }

  # tiny helpers
  unlist1c <- function(x) if (is.list(x)) if (length(x[[1]])) as.character(x[[1]]) else NA_character_ else as.character(x)
  unlist1n <- function(x) if (is.list(x)) if (length(x[[1]])) as.numeric(x[[1]]) else NA_real_ else as.numeric(x)

  grid.newpage()
  grid.text(page_txt, y = unit(0.98, "npc"),
            gp = gpar(fontsize = 12, fontface = "bold"))

  for (i in seq_len(n)) {
    pushViewport(viewport(y = unit(1 - (i - 0.5)/max(1, n), "npc"),
                          height = unit(1/max(1, n), "npc"),
                          x = unit(0.02, "npc"), width = unit(0.96, "npc"),
                          just = c("left", "center")))

    rowi <- tab_top[i]

    gene <- unlist1c(rowi$genesymbol)
    pname <- unlist1c(rowi$providername)
    pid   <- unlist1c(rowi$providerid)
    eff   <- unlist1c(rowi$effect)
    pref  <- unlist1n(rowi$pctref)
    palt  <- unlist1n(rowi$pctalt)
    dlt   <- unlist1n(rowi$delta_pct)
    mpos  <- unlist1n(rowi$motifpos)

    # Find motif matrix by providerName/id
    hit <- integer(0)
    if (!is.na(pname) && !is.na(pid)) {
      hit <- mot_meta$idx[
        tolower(mot_meta$providerName) == tolower(pname) &
          tolower(mot_meta$providerId) == tolower(pid)
      ]
    }
    if (length(hit) == 0L && !is.na(pname)) {
      hit <- mot_meta$idx[tolower(mot_meta$providerName) == tolower(pname)]
    }

    # ASCII title line to avoid device encoding issues
    grid.text(sprintf("%s  |  %s  |  effect=%s  |  Ref=%.3f  Alt=%.3f  d=%.3f",
                      ifelse(is.na(gene)  || gene  == "",  "?", gene),
                      ifelse(is.na(pname) || pname == "", "?", pname),
                      ifelse(is.na(eff)   || eff   == "", "?", eff),
                      pref, palt, dlt),
              x = unit(0, "npc"), y = unit(0.95, "npc"), just = "left",
              gp = gpar(cex = 0.8))

    if (length(hit) > 0L) {
      mat <- as.matrix(motifs[[hit[1]]])  # PPM
      pushViewport(viewport(x = unit(0.35, "npc"), width = unit(0.6, "npc"),
                            y = unit(0.05, "npc"), height = unit(0.75, "npc"),
                            just = c("left","bottom")))
      ok <- TRUE
      tryCatch({
        seqLogo::seqLogo(seqLogo::makePWM(mat))
      }, error = function(e) {
        ok <<- FALSE
      })
      if (ok && !is.na(mpos) && ncol(mat) > 0) {
        xline <- (mpos - 0.5) / ncol(mat)
        if (is.finite(xline) && xline >= 0 && xline <= 1)
          grid.lines(x = unit(c(xline, xline), "npc"),
                     y = unit(c(0, 1), "npc"),
                     gp = gpar(lwd = 1))
      }
      popViewport()
    } else {
      grid.text("(matrix not found for logo)", x = unit(0.35, "npc"),
                y = unit(0.5, "npc"), just = "left", gp = gpar(cex = 0.8))
    }
    popViewport()
  }
}

## 9) Per-SNP CSVs + PDFs (chunked; always at least one page)
all_snps <- unique(names(res))

if (length(all_snps) == 0L) {
  message("No motif hits. ",
          if (use_filterp) "Try a larger p-value (e.g., 1e-3)." else "Try a lower score (e.g., thr=0.80).")
} else {
  for (sid in all_snps) {
    sub <- res[names(res) == sid]
    tab <- to_dt(sub)

    ## CSV sorted by |delta_pct|
    if (nrow(tab)) {
      tab[, abs_delta := abs(delta_pct)]
      data.table::setorder(tab, -abs_delta)
      tab[, abs_delta := NULL]
      csv_path <- file.path(out_dir, paste0(gsub("[:>_]", "-", sid), "_hits.csv"))
      fwrite(tab, csv_path)
    }

    ## PDF (chunking avoids "Too many stacks")
    pdf_path <- file.path(out_dir, paste0(gsub("[:>_]", "-", sid), ".pdf"))
    pdf(pdf_path, width = 8.5, height = 7)

    if (length(sub) > 0L) {
      ord <- if (nrow(tab)) order(abs(tab$delta_pct), decreasing = TRUE) else seq_along(sub)
      sub_sorted <- sub[ord]

      n <- length(sub_sorted)
      pages <- max(1L, ceiling(n / max(1L, as.integer(top_k))))

      plotted_any <- FALSE
      for (p in seq_len(pages)) {
        idx <- ((p - 1L) * top_k + 1L) : min(p * top_k, n)
        sub_chunk <- sub_sorted[idx]
        page_txt <- sprintf(
          "%s  |  %s  |  HOCOMOCO  |  page %d/%d",
          sid, if (use_filterp) paste0("p<= ", p_thresh) else paste0("thr=", thr), p, pages
        )

        ok <- TRUE
        tryCatch({
          motifbreakR::plotMB(
            results = sub_chunk,
            rsid    = sid,
            effect  = c("strong","weak","neutral"),
            main    = page_txt
          )
          plotted_any <- TRUE
        }, error = function(e) {
          ok <<- FALSE
        })

        if (!ok) {
          # Build a chunk table and use robust fallback logos
          tab_chunk <- if (nrow(tab)) {
            tab[idx, .(genesymbol, providername, providerid, effect,
                       motifpos, pctref, pctalt, delta_pct, seqmatch)]
          } else data.table()
          err <- try(draw_fallback_logos(tab_chunk, sid, page_txt), silent = TRUE)
          if (inherits(err, "try-error")) {
            # Guaranteed readable table page
            grid.newpage()
            grid.text(page_txt, y = unit(0.98, "npc"),
                      gp = gpar(fontsize = 12, fontface = "bold"))
            if (nrow(tab_chunk)) {
              tg <- gridExtra::tableGrob(
                tab_chunk[, .(genesymbol, providername, effect, pctref, pctalt, delta_pct)],
                rows = NULL
              )
              pushViewport(viewport(y = 0.45, height = 0.9))
              grid.draw(tg)
              popViewport()
            } else {
              grid.text(sprintf("No plottable panels for %s", sid))
            }
          }
          plotted_any <- TRUE
        }
      }

      if (!plotted_any) {
        grid.newpage()
        grid.text(sprintf("No plottable panels for %s", sid),
                  gp = gpar(fontsize = 12, fontface = "bold"))
      }

    } else {
      grid.newpage()
      grid.text(sprintf("No motif hits for %s", sid),
                gp = gpar(fontsize = 12, fontface = "bold"))
    }

    dev.off()
    message("Wrote: ", basename(pdf_path),
            if (exists("csv_path")) paste0(" & ", basename(csv_path)) else "")
  }
}

message("Done. Outputs in: ", out_dir)
