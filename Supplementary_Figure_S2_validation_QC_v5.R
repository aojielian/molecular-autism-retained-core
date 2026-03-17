#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

root <- "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision"
pkg4  <- file.path(root, "Package4_CrossCohortValidation_Velmeshev")
pkg4v2 <- file.path(root, "Package4_CrossCohortValidation_Velmeshev_v2")
pkg4b <- file.path(root, "Package4b_DirectionConsistencySummary_v2_manual")
pkg5 <- file.path(root, "Package5_DirectionalConcordance_GlobalSummary")

outdir <- file.path(root, "Supplementary_Figures_v2", "S2_Validation_QC_v2")
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
  file.path(pkg4b, "04_direction_consistency_summary_manuscript.tsv"),
  file.path(pkg4b, "03_direction_consistency_summary_pretty.tsv"),
  file.path(pkg5, "01_direction_concordance_master.tsv"),
  file.path(pkg5, "03_metrics_ranked_by_nominal_p.tsv"),
  file.path(pkg5, "04_metrics_ranked_by_aligned_effect.tsv")
)
summary_sources <- summary_sources[file.exists(summary_sources)]
summary_dt <- rbindlist(lapply(summary_sources, function(f) {
  dt <- tryCatch(read_tsv_auto(f), error = function(e) NULL)
  if (is.null(dt)) return(NULL)
  dt[, source_file__ := basename(f)]
  dt
}), fill = TRUE)

find_stats_row <- function(keywords) {
  if (nrow(summary_dt) == 0) return(NULL)
  txt <- apply(summary_dt, 1, function(z) paste(z, collapse = " | "))
  score <- sapply(txt, function(s) {
    s2 <- tolower(s)
    sum(sapply(keywords, function(k) grepl(tolower(k), s2, fixed = TRUE)))
  })
  if (max(score) <= 0) return(NULL)
  summary_dt[which.max(score)]
}

extract_stat <- function(row, candidates) {
  if (is.null(row) || nrow(row) == 0) return(NA_real_)
  hit <- intersect(candidates, names(row))
  if (length(hit) == 0) return(NA_real_)
  as.numeric(row[[hit[1]]][1])
}

make_label <- function(row, beta_override = NA_real_, p_override = NA_real_, fdr_override = NA_real_) {
  beta <- if (!is.na(beta_override)) beta_override else extract_stat(row, c("beta_ASD_vs_Control","beta","observed_effect","effect","aligned_effect"))
  pval <- if (!is.na(p_override)) p_override else extract_stat(row, c("p_value","nominal_p","p","pval"))
  fdr  <- if (!is.na(fdr_override)) fdr_override else extract_stat(row, c("BH_FDR","FDR","fdr","padj","adj_p"))
  paste0("beta = ", fmt_num(beta), "\nP = ", fmt_p(pval), "\nFDR = ", fmt_p(fdr))
}

row_C <- find_stats_row(c("all microglia", "candidate-minus-reference"))
row_D <- find_stats_row(c("candidate proportion", "score-defined"))
row_E <- find_stats_row(c("candidate-reference shift", "score-defined"))

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
