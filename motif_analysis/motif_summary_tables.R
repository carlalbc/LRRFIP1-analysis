## --- Summarize motifbreakR hit CSVs into boolean matrix + counts tables ---
## Reads per-SNP _hits.csv files produced by run_motifbreakr_hocomoco.R
## Outputs: TSV summary tables and per-SNP Excel workbook

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})

## ------------------ Settings ------------------
in_dir  <- "./data/motifbreakr/output"
thr     <- 0.85                 # threshold to call a hit
keep_all_transitions <- FALSE   # TRUE = show all 12 REF>ALT columns
## ----------------------------------------------

message("Reading motifbreakR hit tables from: ", in_dir)

csvs <- list.files(in_dir, pattern = "_hits\\.csv$", full.names = TRUE)
if (!length(csvs)) stop("No *_hits.csv files found in: ", in_dir)

allele_from_filename <- function(fp) {
  b <- basename(fp)
  m <- str_match(b, "-([ACGTN])-([ACGTN])_hits\\.csv$")
  if (is.na(m[1,2]) || is.na(m[1,3])) return(NA_character_)
  paste0(m[1,2], ">", m[1,3])
}

long <- rbindlist(lapply(csvs, function(fp) {
  dt <- fread(fp)
  if (!nrow(dt)) return(NULL)
  setnames(dt, tolower(names(dt)))
  dt[, allele := allele_from_filename(fp)]
  dt[, .(
    genesymbol = as.character(genesymbol),
    pctref     = suppressWarnings(as.numeric(pctref)),
    pctalt     = suppressWarnings(as.numeric(pctalt)),
    delta_pct  = suppressWarnings(as.numeric(delta_pct)),
    effect     = as.character(effect),
    allele     = as.character(allele)
  )]
}), use.names = TRUE, fill = TRUE)

long <- long[!is.na(genesymbol) & nzchar(genesymbol)]
if (!nrow(long)) stop("All hit tables are empty after filtering.")

present_alleles <- sort(unique(long$allele))
present_alleles <- present_alleles[!is.na(present_alleles)]

all_alleles_order <- c(
  "A>C","A>G","A>T",
  "C>A","C>G","C>T",
  "G>A","G>C","G>T",
  "T>A","T>C","T>G"
)
allele_cols <- if (keep_all_transitions) all_alleles_order else intersect(all_alleles_order, present_alleles)

## mark WT/ALT at the chosen threshold
long[, mark_wt  := pctref >= thr]
long[, mark_alt := pctalt >= thr]

wt_marks  <- long[mark_wt  == TRUE, .(Wildtype = TRUE), by = .(genesymbol)]
alt_marks <- long[mark_alt == TRUE, .(mark     = TRUE), by = .(genesymbol, allele)]

## boolean matrix
all_genes <- sort(unique(long$genesymbol))
mat_bool  <- data.table(`MOTIF gene symbol` = all_genes)

mat_bool[, Wildtype := FALSE]
if (nrow(wt_marks)) {
  mat_bool[`MOTIF gene symbol` %in% wt_marks$genesymbol, Wildtype := TRUE]
}

for (ac in allele_cols) mat_bool[, (ac) := FALSE]
if (nrow(alt_marks)) {
  for (ac in allele_cols) {
    genes_hit <- alt_marks[allele == ac, unique(genesymbol)]
    if (length(genes_hit)) {
      mat_bool[`MOTIF gene symbol` %in% genes_hit, (ac) := TRUE]
    }
  }
}

## counts table
counts <- long[, .(Wildtype_count = sum(pctref >= thr, na.rm = TRUE)), by = .(genesymbol)]
setnames(counts, "genesymbol", "MOTIF gene symbol")

for (ac in allele_cols) {
  n_ac <- long[allele == ac, .(cnt = sum(pctalt >= thr, na.rm = TRUE)), by = .(genesymbol)]
  if (nrow(n_ac)) {
    setnames(n_ac, c("genesymbol","cnt"), c("MOTIF gene symbol", paste0(ac, "_count")))
    counts <- merge(counts, n_ac, by = "MOTIF gene symbol", all.x = TRUE)
  } else {
    counts[, (paste0(ac, "_count")) := 0L]
  }
}
for (j in names(counts)) if (grepl("_count$", j)) counts[is.na(get(j)), (j) := 0L]

## convert booleans -> "x"/""
mat_x <- copy(mat_bool)
for (j in setdiff(names(mat_x), "MOTIF gene symbol")) {
  mat_x[is.na(get(j)), (j) := FALSE]
  mat_x[, (j) := ifelse(get(j), "x", "")]
}

setcolorder(mat_x, c("MOTIF gene symbol", "Wildtype", allele_cols))

## write TSVs
out_bool   <- file.path(in_dir, "motif_summary_matrix.tsv")
out_counts <- file.path(in_dir, "motif_summary_counts.tsv")

fwrite(mat_x,   out_bool,   sep = "\t", quote = FALSE, na = "")
fwrite(counts,  out_counts, sep = "\t", quote = FALSE, na = "")

message("\nPreview (matrix):")
print(mat_x)

message("\nPreview (counts):")
print(counts)

cat("\nWrote:\n -", out_bool, "\n -", out_counts, "\n")


## =====================================================================
## Per-SNP Excel workbook with counts + boolean matrix per variant
## =====================================================================

suppressPackageStartupMessages({
  library(openxlsx)
})

out_xlsx <- file.path(in_dir, "motif_summaries_by_snp.xlsx")

## Re-read and stack all hit CSVs
parse_from_filename <- function(fp) {
  b <- basename(fp)
  m <- str_match(b, "^(chr[0-9XYM]+)-(\\d+)-([ACGTN])-([ACGTN])_hits\\.csv$")
  if (any(is.na(m))) return(NULL)
  list(
    chr   = m[2],
    pos   = as.integer(m[3]),
    ref   = m[4],
    alt   = m[5],
    snp_id_pretty = sprintf("%s:%s %s>%s", m[2], m[3], m[4], m[5]),
    sheet_name    = sprintf("%s_%s_%s_%s", m[2], m[3], m[4], m[5]) |> gsub("[:>-]", "_", x = _)
  )
}

long_list <- lapply(csvs, function(fp) {
  p <- parse_from_filename(fp); if (is.null(p)) return(NULL)
  dt <- fread(fp)
  if (!nrow(dt)) return(NULL)
  setnames(dt, tolower(names(dt)))
  keep <- intersect(
    c("genesymbol","pctref","pctalt","delta_pct","effect","snp","seqmatch"),
    names(dt)
  )
  dt <- dt[, ..keep]
  dt[, `:=`(
    snp_file   = fp,
    snp_pretty = p$snp_id_pretty,
    allele     = paste0(substr(p$ref,1,1), ">", substr(p$alt,1,1))
  )]
  for (nm in c("pctref","pctalt","delta_pct")) if (nm %in% names(dt)) dt[, (nm) := as.numeric(get(nm))]
  dt[]
})

long2 <- rbindlist(long_list, use.names = TRUE, fill = TRUE)
long2 <- long2[!is.na(genesymbol) & nzchar(genesymbol)]
long2[, genesymbol := as.character(genesymbol)]

## Per-SNP summary builder
build_per_snp <- function(df, thr, snp_pretty) {
  df[, wt_hit  := !is.na(pctref) & pctref >= thr]
  df[, alt_hit := !is.na(pctalt) & pctalt >= thr]
  allele_str <- unique(df$allele)
  if (length(allele_str) != 1L) allele_str <- sort(allele_str)[1]

  counts <- df[, .(
    Wildtype_count = sum(wt_hit,  na.rm = TRUE),
    Alt_count      = sum(alt_hit, na.rm = TRUE)
  ), by = .(genesymbol)]
  setorder(counts, -Alt_count, -Wildtype_count, genesymbol)
  alt_col <- paste0(allele_str, "_count")
  setnames(counts, "Alt_count", alt_col)

  bool <- df[, .(
    Wildtype = any(wt_hit,  na.rm = TRUE),
    Alt      = any(alt_hit, na.rm = TRUE)
  ), by = .(genesymbol)]
  setorder(bool, -Alt, -Wildtype, genesymbol)
  setnames(bool, "Alt", allele_str)

  for (nm in c("Wildtype", allele_str)) {
    bool[, (nm) := ifelse(isTRUE(get(nm)), "x", "")]
  }

  counts[, SNP := snp_pretty]; setcolorder(counts, c("SNP", "genesymbol", "Wildtype_count", alt_col))
  setnames(counts, c("genesymbol","Wildtype_count"), c("MOTIF gene symbol","Wildtype_count"))
  bool[,   SNP := snp_pretty]; setcolorder(bool,   c("SNP", "genesymbol", "Wildtype", allele_str))
  setnames(bool,   c("genesymbol","Wildtype"), c("MOTIF gene symbol","Wildtype"))

  list(counts = counts, bool = bool)
}

## Build workbook
wb <- createWorkbook()
all_counts <- list()
all_bool   <- list()

by_snp <- split(long2, long2$snp_pretty)
for (snp_pretty in names(by_snp)) {
  df <- by_snp[[snp_pretty]]
  first_fp <- df$snp_file[1]
  p <- parse_from_filename(first_fp)
  sheet <- substr(p$sheet_name, 1, 31)

  tabs <- build_per_snp(df, thr = thr, snp_pretty = snp_pretty)
  cnt  <- tabs$counts
  boo  <- tabs$bool

  addWorksheet(wb, sheet)
  writeData(wb, sheet, paste0("SNP: ", snp_pretty), startRow = 1, startCol = 1)
  writeDataTable(wb, sheet, cnt, startRow = 3, startCol = 1, tableStyle = "TableStyleMedium9", withFilter = TRUE)
  start_row_2 <- 5 + nrow(cnt)
  writeDataTable(wb, sheet, boo, startRow = start_row_2, startCol = 1, tableStyle = "TableStyleMedium2", withFilter = TRUE)

  all_counts[[snp_pretty]] <- cnt
  all_bool[[snp_pretty]]   <- boo
}

## Overview sheets
overview_counts <- rbindlist(all_counts, use.names = TRUE, fill = TRUE)
overview_bool   <- rbindlist(all_bool,   use.names = TRUE, fill = TRUE)

addWorksheet(wb, "OVERVIEW_counts")
writeDataTable(wb, "OVERVIEW_counts", overview_counts, tableStyle = "TableStyleMedium9", withFilter = TRUE)

addWorksheet(wb, "OVERVIEW_bool")
writeDataTable(wb, "OVERVIEW_bool", overview_bool, tableStyle = "TableStyleMedium2", withFilter = TRUE)

saveWorkbook(wb, out_xlsx, overwrite = TRUE)

cat("Wrote Excel:\n -", out_xlsx, "\n")
