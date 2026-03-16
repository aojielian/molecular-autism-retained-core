#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(optparse)
})


project_dir_default <- Sys.getenv("MA_PROJECT_DIR", ".")
option_list <- list(
  make_option("--root", type = "character", default = project_dir_default, help = "Project directory [default: %default]"),
  make_option("--pkg4", type = "character", default = file.path(project_dir_default, "results", "Package4_CrossCohortValidation_Velmeshev"), help = "Package4 directory [default: %default]"),
  make_option("--pkg4v2", type = "character", default = file.path(project_dir_default, "results", "Package4_CrossCohortValidation_Velmeshev_v2"), help = "Package4 v2 directory [default: %default]"),
  make_option("--pkg4b", type = "character", default = file.path(project_dir_default, "results", "Package4b_DirectionConsistencySummary_v2_manual"), help = "Package4b directory [default: %default]"),
  make_option("--pkg5", type = "character", default = file.path(project_dir_default, "results", "Package5_DirectionalConcordance_GlobalSummary"), help = "Package5 directory [default: %default]"),
  make_option("--outdir", type = "character", default = file.path(project_dir_default, "results", "Supplementary_Figures", "S3_Validation_QC"), help = "Output directory [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
root <- opt$root
pkg4 <- opt$pkg4
pkg4v2 <- opt$pkg4v2
pkg4b <- opt$pkg4b
pkg5 <- opt$pkg5
outdir <- opt$outdir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

read_tsv_auto <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
  if (grepl("\\.gz$", path)) {
    fread(cmd = paste("gzip -dc", shQuote(path)), sep = "\t", header = TRUE, data.table = TRUE)
  } else {
    fread(path, sep = "\t", header = TRUE, data.table = TRUE)
  }
}

find_col <- function(dt, candidates, required = TRUE) {
  hit <- candidates[candidates %in% names(dt)]
  if (length(hit) > 0) return(hit[1])
  if (required) stop("Cannot find column. Candidates tried: ", paste(candidates, collapse = ", "),
                     "\nAvailable: ", paste(names(dt), collapse = ", "))
  return(NULL)
}

find_metric_col <- function(dt, preferred, exclude = c("sample","Subject","subject","dx","dx_std","Diagnosis","diagnosis","group","Group",
                                                       "n_total","n_candidate","n_reference","prop_reference",
                                                       "SeqBatch","Sex","Ancestry","Age","Age_sqd","PMI","RIN")) {
  hit <- preferred[preferred %in% names(dt)]
  if (length(hit) > 0) return(hit[1])
  num_cols <- names(dt)[vapply(dt, is.numeric, logical(1))]
  num_cols <- setdiff(num_cols, exclude)
  if (length(num_cols) == 1) return(num_cols[1])
  if (length(num_cols) > 1) {
    # Prefer columns with informative names
    ord <- order(!grepl("candidate|reference|shift|prop|ratio|score", num_cols, ignore.case = TRUE), num_cols)
    return(num_cols[ord][1])
  }
  stop("Cannot infer metric column. Available numeric columns: ", paste(names(dt)[vapply(dt, is.numeric, logical(1))], collapse = ", "))
}

std_dx <- function(x) {
  x <- as.character(x)
  ifelse(x %in% c("ASD","Case","case","Autism"), "ASD",
         ifelse(x %in% c("Control","CTL","control"), "Control", x))
}

fmt_num <- function(x, digits = 3) {
  if (length(x) == 0 || all(is.na(x))) return("NA")
  sprintf(paste0("%.", digits, "f"), as.numeric(x[1]))
}

fmt_p <- function(x) {
  if (length(x) == 0 || all(is.na(x))) return("NA")
  x <- as.numeric(x[1])
  if (!is.finite(x)) return("NA")
  if (x < 1e-4) return(formatC(x, format = "e", digits = 2))
  sprintf("%.3f", x)
}

theme_clean <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    )
}

col_control <- "#C9C9C9"
col_asd <- "#D7B0B7"

# ---------- Read core inputs ----------
umap <- read_tsv_auto(file.path(pkg4, "06_validation_microglia_umap.tsv.gz"))
dx_counts <- read_tsv_auto(file.path(pkg4v2, "03_validation_microglia_dx_counts.tsv"))
all_micro <- read_tsv_auto(file.path(pkg4v2, "14_validation_sample_scores_all_microglia.tsv"))
score_prop <- read_tsv_auto(file.path(pkg4v2, "15_validation_sample_score_defined_proportions.tsv"))
score_shift <- read_tsv_auto(file.path(pkg4v2, "16_validation_sample_score_defined_shift.tsv"))

# ---------- Standardize / detect columns ----------
umap_x <- find_col(umap, c("UMAP_1","UMAP1","umap_1","x"))
umap_y <- find_col(umap, c("UMAP_2","UMAP2","umap_2","y"))
subcluster_col <- find_col(umap, c("val_subcluster","subcluster","Subcluster","cluster","seurat_clusters","validation_subcluster"))
umap[[subcluster_col]] <- as.character(umap[[subcluster_col]])

dx_col_counts <- find_col(dx_counts, c("dx","Diagnosis","diagnosis","group","Group","dx_std"))
n_col_counts <- find_col(dx_counts, c("N","n","count","cells"))
dx_counts[, dx_std := std_dx(get(dx_col_counts))]

dx_col_all <- find_col(all_micro, c("dx","Diagnosis","diagnosis","group","Group","dx_std"))
all_micro[, dx_std := std_dx(get(dx_col_all))]
all_val_col <- find_col(all_micro, c("candidate_minus_reference","candidate_minus_reference_score","score","candidate_reference","candidate_minus_reference_all_microglia"))

dx_col_prop <- find_col(score_prop, c("dx","Diagnosis","diagnosis","group","Group","dx_std"))
score_prop[, dx_std := std_dx(get(dx_col_prop))]
prop_val_col <- find_metric_col(score_prop, c("prop_candidate","candidate_proportion","candidate_prop","proportion","score_defined_candidate_proportion"))

dx_col_shift <- find_col(score_shift, c("dx","Diagnosis","diagnosis","group","Group","dx_std"))
score_shift[, dx_std := std_dx(get(dx_col_shift))]
shift_val_col <- find_metric_col(score_shift, c("candidate_minus_reference_shift","candidate_vs_reference_shift","candidate_reference_shift","shift","score_defined_shift"))

# ---------- Stats lookup ----------
summary_sources <- c(
  file.path(pkg4b, "03_direction_consistency_summary_pretty.tsv"),
  file.path(pkg4b, "04_direction_consistency_summary_manuscript.tsv")
)
summary_sources <- summary_sources[file.exists(summary_sources)]
summary_dt <- rbindlist(lapply(summary_sources, function(f) {
  dt <- tryCatch(read_tsv_auto(f), error = function(e) NULL)
  if (is.null(dt)) return(NULL)
  dt[, source_file__ := basename(f)]
  dt
}), fill = TRUE)

extract_stat <- function(row, candidates) {
  if (is.null(row) || nrow(row) == 0) return(NA_real_)
  hit <- intersect(candidates, names(row))
  if (length(hit) == 0) return(NA_real_)
  as.numeric(row[[hit[1]]][1])
}

find_metric_row <- function(metric_values) {
  if (nrow(summary_dt) == 0) return(NULL)
  # prefer exact match in metric_label column if available
  if ("metric_label" %in% names(summary_dt)) {
    for (mv in metric_values) {
      hit <- summary_dt[metric_label == mv]
      if (nrow(hit) > 0) return(hit[1])
    }
  }
  # fallback: scan rows as text
  rowtxt <- apply(summary_dt, 1, function(z) paste(z, collapse = " | "))
  for (mv in metric_values) {
    idx <- which(grepl(tolower(mv), tolower(rowtxt), fixed = TRUE))
    if (length(idx) > 0) return(summary_dt[idx[1]])
  }
  return(NULL)
}

make_label <- function(row) {
  beta <- extract_stat(row, c("beta_ASD_vs_Control","beta","observed_effect","effect_case_minus_control","effect","aligned_effect"))
  pval <- extract_stat(row, c("lm_p","p_value","nominal_p","wilcox_p","p","pval"))
  fdr  <- extract_stat(row, c("BH_FDR","FDR","fdr","padj","adj_p"))
  if (is.finite(fdr)) {
    paste0("beta = ", fmt_num(beta), "\nP = ", fmt_p(pval), "\nFDR = ", fmt_p(fdr))
  } else {
    paste0("beta = ", fmt_num(beta), "\nP = ", fmt_p(pval))
  }
}

row_C <- find_metric_row(c("All microglia: candidate - reference", "candidate - reference"))
row_D <- find_metric_row(c("Score-defined cells: candidate proportion", "candidate proportion"))
row_E <- find_metric_row(c("Score-defined cells: candidate vs reference shift", "candidate vs reference shift"))

lab_C <- make_label(row_C)
lab_D <- make_label(row_D)
lab_E <- make_label(row_E)

# ---------- Panels ----------
pA <- ggplot(umap, aes(x = .data[[umap_x]], y = .data[[umap_y]], color = .data[[subcluster_col]])) +
  geom_point(size = 0.6, alpha = 0.85, show.legend = FALSE) +
  theme_clean(12) +
  labs(title = "A. Validation microglia UMAP by subcluster", x = "UMAP1", y = "UMAP2")

pB <- ggplot(dx_counts, aes(x = dx_std, y = .data[[n_col_counts]], fill = dx_std)) +
  geom_col(width = 0.62, color = NA) +
  scale_fill_manual(values = c("Control" = col_control, "ASD" = "#C77E86")) +
  theme_clean(12) +
  theme(legend.position = "none") +
  labs(title = "B. Validation microglial cell counts by diagnosis",
       x = NULL, y = "Number of cells")

make_violin_panel <- function(dt, ycol, ylab, title, ann_label) {
  y_max <- max(dt[[ycol]], na.rm = TRUE)
  y_min <- min(dt[[ycol]], na.rm = TRUE)
  y_pos <- y_max - 0.12 * (y_max - y_min)
  ggplot(dt, aes(x = dx_std, y = .data[[ycol]], fill = dx_std)) +
    geom_violin(trim = FALSE, width = 0.95, color = NA) +
    geom_boxplot(width = 0.18, outlier.shape = NA, fill = "white") +
    geom_jitter(width = 0.08, size = 1.8, alpha = 0.9) +
    scale_fill_manual(values = c("Control" = col_control, "ASD" = col_asd)) +
    theme_clean(12) +
    theme(legend.position = "none") +
    annotate("label", x = 1.5, y = y_pos, label = ann_label, size = 3.2) +
    labs(title = title, x = NULL, y = ylab)
}

pC <- make_violin_panel(all_micro, all_val_col, "Score",
                        "C. All-microglia candidate-minus-reference", lab_C)
pD <- make_violin_panel(score_prop, prop_val_col, "Proportion",
                        "D. Score-defined candidate proportion", lab_D)
pE <- make_violin_panel(score_shift, shift_val_col, "Shift",
                        "E. Score-defined candidate-versus-reference shift", lab_E)

fig <- (pA | pB) / (pC | pD | pE) +
  plot_layout(widths = c(1.15, 0.95), heights = c(1, 1))

ggsave(file.path(outdir, "Supplementary_Figure_S2_validation_QC.png"), fig, width = 16, height = 11.5, dpi = 300)
ggsave(file.path(outdir, "Supplementary_Figure_S2_validation_QC.pdf"), fig, width = 16, height = 11.5)

message("Written to: ", outdir)
