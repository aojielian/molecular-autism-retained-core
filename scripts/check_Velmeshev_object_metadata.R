#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(data.table)
})

option_list <- list(
  make_option("--rds", type = "character",
              default = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/Velmeshev_scRNA/Velmeshev_Object.rds"),
  make_option("--outdir", type = "character",
              default = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Velmeshev_metadata_check")
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(opt$outdir, "00_check_log.txt")
log_msg <- function(...) {
  msg <- sprintf(...)
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

safe_fwrite <- function(x, file) {
  if (grepl("\\.gz$", file)) {
    fwrite(x, file = file, sep = "\t", quote = FALSE, na = "NA", compress = "gzip")
  } else {
    fwrite(x, file = file, sep = "\t", quote = FALSE, na = "NA")
  }
}

pick_first_existing <- function(nms, candidates) {
  hit <- candidates[candidates %in% nms]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

is_dx_like <- function(x) {
  x <- unique(tolower(trimws(as.character(x))))
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(FALSE)
  patterns <- c("asd", "autism", "control", "ctrl", "ctl", "case", "normal", "neurotyp", "hc")
  any(vapply(patterns, function(p) any(grepl(p, x)), logical(1)))
}

guess_columns <- function(meta) {
  nms <- colnames(meta)

  cluster_candidates <- c(
    "cluster_use", "seurat_clusters", "cluster", "Cluster", "subcluster", "ident",
    "celltype", "CellType", "cell_type", "annotation", "Annotation",
    "subclass", "class", "broad_cell_type"
  )
  donor_candidates <- c(
    "donor_id", "donor", "individualID", "individual_id", "subject_id", "subject",
    "projid", "sample_id", "orig.ident", "sample", "patient", "Patient"
  )
  dx_candidates <- c(
    "Diagnosis", "diagnosis", "Dx", "dx", "Group", "group", "condition", "Condition",
    "phenotype", "Phenotype", "disease", "Disease"
  )

  cluster_col <- pick_first_existing(nms, cluster_candidates)
  donor_col   <- pick_first_existing(nms, donor_candidates)
  dx_col      <- pick_first_existing(nms, dx_candidates)

  if (is.na(dx_col)) {
    dx_hits <- nms[vapply(nms, function(nm) is_dx_like(meta[[nm]]), logical(1))]
    if (length(dx_hits) > 0) dx_col <- dx_hits[1]
  }

  data.table(
    guessed_cluster_col = cluster_col,
    guessed_donor_col = donor_col,
    guessed_dx_col = dx_col
  )
}

summarize_column <- function(vec) {
  x <- as.character(vec)
  x_non_na <- x[!is.na(x) & x != ""]
  uniq <- unique(x_non_na)
  data.table(
    class = class(vec)[1],
    n_unique = length(uniq),
    n_na = sum(is.na(vec)),
    example_values = paste(head(uniq, 10), collapse = " | ")
  )
}

obj <- readRDS(opt$rds)
log_msg("Loaded: %s", opt$rds)
log_msg("Object class: %s", paste(class(obj), collapse = ", "))

assays <- tryCatch(Assays(obj), error = function(e) character())
reductions <- tryCatch(Reductions(obj), error = function(e) character())
default_assay <- tryCatch(DefaultAssay(obj), error = function(e) NA_character_)

meta <- as.data.table(obj@meta.data, keep.rownames = "cell")

# ------------------------------------------------------------
# Basic object summary
# ------------------------------------------------------------
basic_summary <- data.table(
  object_class = paste(class(obj), collapse = ", "),
  n_cells = nrow(meta),
  n_metadata_columns = ncol(meta) - 1,
  default_assay = default_assay,
  assays = paste(assays, collapse = ", "),
  reductions = paste(reductions, collapse = ", ")
)
safe_fwrite(basic_summary, file.path(opt$outdir, "01_basic_summary.tsv"))

# Feature dimensions per assay
assay_summary <- rbindlist(lapply(assays, function(a) {
  nc <- tryCatch(ncol(obj[[a]]), error = function(e) NA_integer_)
  nr <- tryCatch(nrow(obj[[a]]), error = function(e) NA_integer_)
  data.table(assay = a, n_features = nr, n_cells = nc)
}), fill = TRUE)
safe_fwrite(assay_summary, file.path(opt$outdir, "02_assay_summary.tsv"))

# ------------------------------------------------------------
# Metadata column summary
# ------------------------------------------------------------
meta_cols <- setdiff(colnames(meta), "cell")
meta_summary <- rbindlist(lapply(meta_cols, function(nm) {
  sm <- summarize_column(meta[[nm]])
  data.table(column = nm, sm)
}), fill = TRUE)
safe_fwrite(meta_summary, file.path(opt$outdir, "03_metadata_column_summary.tsv"))

# ------------------------------------------------------------
# Suggested columns
# ------------------------------------------------------------
guesses <- guess_columns(meta)
safe_fwrite(guesses, file.path(opt$outdir, "04_suggested_columns.tsv"))

# ------------------------------------------------------------
# Candidate column counts
# ------------------------------------------------------------
candidate_cols <- unique(na.omit(c(
  guesses$guessed_cluster_col,
  guesses$guessed_donor_col,
  guesses$guessed_dx_col
)))

if (length(candidate_cols) > 0) {
  for (nm in candidate_cols) {
    dt <- meta[, .N, by = .(value = as.character(get(nm)))][order(-N)]
    safe_fwrite(dt, file.path(opt$outdir, sprintf("05_counts__%s.tsv", nm)))
  }
}

# ------------------------------------------------------------
# Heuristic extra inspection
# ------------------------------------------------------------
# likely cluster-like columns = low/moderate cardinality and repeated values
cluster_like <- meta_summary[
  n_unique >= 2 & n_unique <= 100 &
    !grepl("^nCount_|^nFeature_|percent|scrublet|mt|mito|RNA|ATAC|SCT", column, ignore.case = TRUE)
][order(n_unique, column)]

safe_fwrite(cluster_like, file.path(opt$outdir, "06_cluster_like_columns.tsv"))

# likely dx-like columns by content
dx_like_cols <- meta_cols[vapply(meta_cols, function(nm) is_dx_like(meta[[nm]]), logical(1))]
dx_like_dt <- data.table(dx_like_columns = dx_like_cols)
safe_fwrite(dx_like_dt, file.path(opt$outdir, "07_dx_like_columns.tsv"))

# likely donor-like columns = high cardinality but fewer than cells
donor_like <- meta_summary[
  n_unique >= 5 & n_unique < nrow(meta) &
    !grepl("^cluster|seurat_clusters|celltype|annotation|class|subclass", column, ignore.case = TRUE)
][order(-n_unique, column)]
safe_fwrite(donor_like, file.path(opt$outdir, "08_donor_like_columns.tsv"))

# ------------------------------------------------------------
# Print key messages
# ------------------------------------------------------------
log_msg("Default assay: %s", default_assay)
log_msg("Assays: %s", paste(assays, collapse = ", "))
log_msg("Reductions: %s", paste(reductions, collapse = ", "))
log_msg("n_cells: %d", nrow(meta))
log_msg("n_metadata_columns: %d", ncol(meta) - 1)
log_msg("Suggested cluster col: %s", guesses$guessed_cluster_col)
log_msg("Suggested donor col: %s", guesses$guessed_donor_col)
log_msg("Suggested dx col: %s", guesses$guessed_dx_col)

cat("\n===== KEY SUMMARY =====\n")
print(basic_summary)
cat("\n===== SUGGESTED COLUMNS =====\n")
print(guesses)

if (length(candidate_cols) > 0) {
  cat("\n===== TOP COUNTS OF SUGGESTED COLUMNS =====\n")
  for (nm in candidate_cols) {
    cat("\n---", nm, "---\n")
    print(head(meta[, .N, by = .(value = as.character(get(nm)))][order(-N)], 20))
  }
}
