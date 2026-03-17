#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(optparse)
})

# =========================================================
# Configuration
# =========================================================

project_dir_default <- Sys.getenv("MA_PROJECT_DIR", ".")
option_list <- list(
  make_option("--indir", type = "character", default = file.path(project_dir_default, "results", "Package3_Discovery_DonorAwareDiseaseAssociation"), help = "Package3 output directory [default: %default]"),
  make_option("--outdir", type = "character", default = file.path(project_dir_default, "results", "Figure4_Pseudobulk_Remodeling"), help = "Output directory [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
cfg <- list(
  indir = opt$indir,
  outdir = opt$outdir,
  pooled_file = "13_pseudobulk_all_retained_ASD_vs_Control.tsv.gz",
  cluster0_file = "15_pseudobulk_cluster0_ASD_vs_Control.tsv.gz",
  cluster2_file = "17_pseudobulk_cluster2_ASD_vs_Control.tsv.gz",
  fdr_cut = 0.10,
  p_cut = 0.05,
  lfc_label_min = 0.35,
  case_color = "#C97B84",
  ctrl_color = "#7A8796",
  neutral_color = "#CFCFCF"
)

dir.create(cfg$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cfg$outdir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cfg$outdir, "tables"), recursive = TRUE, showWarnings = FALSE)

# =========================================================
# Helpers
# =========================================================
read_tsv_auto <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  if (grepl("\\.gz$", path)) {
    fread(cmd = paste("gzip -dc", shQuote(path)), sep = "\t", header = TRUE, data.table = TRUE)
  } else {
    fread(path, sep = "\t", header = TRUE, data.table = TRUE)
  }
}

pick_col <- function(dt, candidates, required = TRUE) {
  hit <- candidates[candidates %in% colnames(dt)]
  if (length(hit) > 0) return(hit[1])
  if (required) stop("Cannot find required column. Available columns: ", paste(colnames(dt), collapse = ", "))
  return(NULL)
}

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, digits = digits, format = "f"))
}

prep_de <- function(dt, analysis_name) {
  dt <- as.data.table(copy(dt))
  gene_col <- pick_col(dt, c("gene", "Gene", "symbol"))
  fc_col   <- pick_col(dt, c("logFC", "log2FC", "avg_log2FC", "avg_logFC"))
  p_col    <- pick_col(dt, c("PValue", "p_val", "pvalue", "P.Value"))
  fdr_col  <- pick_col(dt, c("FDR", "padj", "adj.P.Val", "qvalue"))

  setnames(dt, c(gene_col, fc_col, p_col, fdr_col), c("gene", "logFC", "PValue", "FDR"))
  dt[, gene := as.character(gene)]
  dt <- dt[!is.na(gene) & !is.na(logFC) & !is.na(PValue) & !is.na(FDR)]
  dt[, neglog10p := -log10(pmax(PValue, 1e-300))]
  dt[, direction := fifelse(FDR < cfg$fdr_cut & logFC > 0, "ASD_up",
                     fifelse(FDR < cfg$fdr_cut & logFC < 0, "Control_up", "NS"))]
  dt[, analysis := analysis_name]
  dt
}

make_stats <- function(dt) {
  data.table(
    analysis = unique(dt$analysis)[1],
    n_genes = nrow(dt),
    n_fdr10_up = sum(dt$FDR < cfg$fdr_cut & dt$logFC > 0, na.rm = TRUE),
    n_fdr10_down = sum(dt$FDR < cfg$fdr_cut & dt$logFC < 0, na.rm = TRUE),
    min_fdr = min(dt$FDR, na.rm = TRUE),
    top_up_gene = if (sum(dt$logFC > 0, na.rm = TRUE) > 0) dt[logFC > 0][order(FDR, -abs(logFC))]$gene[1] else NA_character_,
    top_down_gene = if (sum(dt$logFC < 0, na.rm = TRUE) > 0) dt[logFC < 0][order(FDR, -abs(logFC))]$gene[1] else NA_character_
  )
}

choose_labels <- function(dt, must_label = character(), max_labels = NULL) {
  dt <- copy(dt)
  must_label <- unique(must_label[must_label %in% dt$gene])
  dt[, label_flag := FALSE]
  if (length(must_label) > 0) dt[gene %in% must_label, label_flag := TRUE]

  if (is.null(max_labels)) max_labels <- length(must_label)
  remain_slots <- max(0, max_labels - sum(dt$label_flag))
  extra <- character(0)
  if (remain_slots > 0) {
    extra <- dt[!label_flag & abs(logFC) >= cfg$lfc_label_min][order(FDR, -abs(logFC))][seq_len(min(.N, remain_slots)), gene]
  }
  if (length(extra) > 0) dt[gene %in% extra, label_flag := TRUE]
  dt[label_flag == TRUE]
}

add_labels <- function(p, lab_dt, size = 3.1) {
  if (nrow(lab_dt) == 0) return(p)
  if (requireNamespace("ggrepel", quietly = TRUE)) {
    p + ggrepel::geom_text_repel(
      data = lab_dt,
      aes(x = logFC, y = neglog10p, label = gene),
      size = size,
      box.padding = 0.25,
      point.padding = 0.15,
      min.segment.length = 0,
      seed = 1,
      max.overlaps = Inf
    )
  } else {
    p + geom_text(
      data = lab_dt,
      aes(x = logFC, y = neglog10p, label = gene),
      size = size,
      check_overlap = TRUE,
      vjust = -0.4
    )
  }
}

plot_volcano <- function(dt, title, must_label, max_labels = length(must_label)) {
  lab_dt <- choose_labels(dt, must_label = must_label, max_labels = max_labels)
  x_lim <- max(abs(dt$logFC), na.rm = TRUE)
  x_lim <- max(1.2, x_lim * 1.12)
  y_lim <- max(dt$neglog10p, na.rm = TRUE)
  y_lim <- max(4.0, y_lim * 1.08)

  stats <- make_stats(dt)
  ann <- paste0(
    "ASD-up genes at FDR < 0.10: ", stats$n_fdr10_up,
    "\nControl-up genes at FDR < 0.10: ", stats$n_fdr10_down,
    "\nMinimum FDR: ", fmt_num(stats$min_fdr, 3)
  )

  p <- ggplot(dt, aes(x = logFC, y = neglog10p)) +
    geom_point(aes(color = direction), size = 1.2, alpha = 0.85) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey45") +
    geom_hline(yintercept = -log10(cfg$p_cut), linetype = 3, color = "grey60") +
    scale_color_manual(values = c(ASD_up = cfg$case_color, Control_up = cfg$ctrl_color, NS = cfg$neutral_color)) +
    coord_cartesian(xlim = c(-x_lim, x_lim), ylim = c(0, y_lim), clip = "off") +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    ) +
    labs(
      title = title,
      x = expression(log[2]~fold-change~"[ASD vs Control]"),
      y = expression(-log[10](P))
    ) +
    annotate(
      "label",
      x = -0.98 * x_lim,
      y = 0.98 * y_lim,
      hjust = 0,
      vjust = 1,
      size = 3.0,
      label = ann
    )

  add_labels(p, lab_dt, size = 3.1)
}

# =========================================================
# Load data
# =========================================================
pooled <- prep_de(
  read_tsv_auto(file.path(cfg$indir, cfg$pooled_file)),
  analysis_name = "pooled_retained"
)
cluster0 <- prep_de(
  read_tsv_auto(file.path(cfg$indir, cfg$cluster0_file)),
  analysis_name = "cluster0"
)
cluster2 <- prep_de(
  read_tsv_auto(file.path(cfg$indir, cfg$cluster2_file)),
  analysis_name = "cluster2"
)

# Save stats
stats_dt <- rbindlist(list(make_stats(pooled), make_stats(cluster0), make_stats(cluster2)), fill = TRUE)
fwrite(stats_dt, file.path(cfg$outdir, "tables", "Figure4_DE_summary_stats.tsv"), sep = "\t")

# Save top tables for convenience
fwrite(pooled[order(FDR, -abs(logFC))][1:min(.N, 200)], file.path(cfg$outdir, "tables", "Figure4_pooled_top200.tsv"), sep = "\t")
fwrite(cluster0[order(FDR, -abs(logFC))][1:min(.N, 200)], file.path(cfg$outdir, "tables", "Figure4_cluster0_top200.tsv"), sep = "\t")
fwrite(cluster2[order(FDR, -abs(logFC))][1:min(.N, 200)], file.path(cfg$outdir, "tables", "Figure4_cluster2_top200.tsv"), sep = "\t")

# =========================================================
# Label sets aligned with manuscript narrative
# =========================================================
labels_pooled <- c("SPP1", "FCGR2A", "ACSL1", "HIF1A", "P4HA1", "SOAT1", "HSPA1A", "MYO1E")
labels_cluster0 <- c("SPP1", "STAT3", "CHST15", "ZC3HAV1", "P4HA1", "HIF1A", "LY96")
labels_cluster2 <- c("SPP1", "ZC3HAV1", "B2M", "MYO1E", "ASAH1", "P4HA1", "SKAP2")

# =========================================================
# Plot
# =========================================================
pA <- plot_volcano(
  pooled,
  title = "A. Pooled retained core microglia show ASD-associated remodeling",
  must_label = labels_pooled,
  max_labels = length(labels_pooled)
)

pB <- plot_volcano(
  cluster0,
  title = "B. Cluster 0 also carries ASD-associated remodeling",
  must_label = labels_cluster0,
  max_labels = length(labels_cluster0)
)

pC <- plot_volcano(
  cluster2,
  title = "C. Cluster 2 also carries ASD-associated remodeling",
  must_label = labels_cluster2,
  max_labels = length(labels_cluster2)
)

fig4 <- pA / pB / pC + plot_layout(heights = c(1, 1, 1))

# =========================================================
# Save
# =========================================================
ggsave(file.path(cfg$outdir, "plots", "Figure4_combined.pdf"), fig4, width = 11.5, height = 16.0)
ggsave(file.path(cfg$outdir, "plots", "Figure4_combined.png"), fig4, width = 11.5, height = 16.0, dpi = 300)

ggsave(file.path(cfg$outdir, "plots", "Figure4A_pooled_volcano.pdf"), pA, width = 11.0, height = 5.2)
ggsave(file.path(cfg$outdir, "plots", "Figure4B_cluster0_volcano.pdf"), pB, width = 11.0, height = 5.2)
ggsave(file.path(cfg$outdir, "plots", "Figure4C_cluster2_volcano.pdf"), pC, width = 11.0, height = 5.2)

ggsave(file.path(cfg$outdir, "plots", "Figure4A_pooled_volcano.png"), pA, width = 11.0, height = 5.2, dpi = 300)
ggsave(file.path(cfg$outdir, "plots", "Figure4B_cluster0_volcano.png"), pB, width = 11.0, height = 5.2, dpi = 300)
ggsave(file.path(cfg$outdir, "plots", "Figure4C_cluster2_volcano.png"), pC, width = 11.0, height = 5.2, dpi = 300)

cat("Figure 4 done. Output directory:\n", cfg$outdir, "\n", sep = "")
