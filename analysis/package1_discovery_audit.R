#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(Matrix)
  library(patchwork)
  library(optparse)
})

options(stringsAsFactors = FALSE)

DEFAULT_PROJECT_DIR <- normalizePath(Sys.getenv("MA_PROJECT_DIR", "."), mustWork = FALSE)

option_list <- list(
  make_option("--microglia_rds", type = "character", default = Sys.getenv("MA_MICROGLIA_RDS", NA_character_),
              help = "Path to discovery microglia Seurat object (.rds) [required or set MA_MICROGLIA_RDS]"),
  make_option("--outdir", type = "character",
              default = file.path(DEFAULT_PROJECT_DIR, "results", "Package1_Figure1_DiscoveryAudit"),
              help = "Output directory [default: %default]"),
  make_option("--assay", type = "character", default = "RNA"),
  make_option("--cluster_col", type = "character", default = NA_character_),
  make_option("--donor_col", type = "character", default = NA_character_),
  make_option("--dx_col", type = "character", default = NA_character_),
  make_option("--top_n", type = "integer", default = 50),
  make_option("--min_pct", type = "double", default = 0.10),
  make_option("--logfc", type = "double", default = 0.25)
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.na(opt$microglia_rds) || opt$microglia_rds == "") stop("Please provide --microglia_rds or set MA_MICROGLIA_RDS.")

cfg <- list(
  microglia_rds = opt$microglia_rds,
  outdir        = opt$outdir,
  assay         = opt$assay,
  cluster_col   = if (is.na(opt$cluster_col) || opt$cluster_col == "") NULL else opt$cluster_col,
  donor_col     = if (is.na(opt$donor_col) || opt$donor_col == "") NULL else opt$donor_col,
  dx_col        = if (is.na(opt$dx_col) || opt$dx_col == "") NULL else opt$dx_col,
  top_n         = opt$top_n,
  min_pct       = opt$min_pct,
  logfc         = opt$logfc
)

dir.create(cfg$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cfg$outdir, "Cluster_Signatures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cfg$outdir, "FeaturePlots"), recursive = TRUE, showWarnings = FALSE)

msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ..., "\n", sep = "")
}

save_session_info <- function(path) {
  writeLines(capture.output(sessionInfo()), con = path)
}

find_first_col <- function(df, patterns) {
  nm <- names(df)
  for (p in patterns) {
    hit <- nm[grepl(p, nm, ignore.case = TRUE)]
    if (length(hit) > 0) return(hit[1])
  }
  NULL
}

normalize_dx <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  y <- x
  y[grepl("asd|autism|case", x, ignore.case = TRUE)] <- "ASD"
  y[grepl("control|ctrl|ctl|healthy|non[- ]?psy", x, ignore.case = TRUE)] <- "Control"
  y
}

choose_cluster_col <- function(obj, user_col = NULL) {
  meta <- obj@meta.data
  if (!is.null(user_col) && user_col %in% names(meta)) return(user_col)
  if ("seurat_clusters" %in% names(meta)) return("seurat_clusters")
  cand <- find_first_col(meta, c("^cluster$", "cluster", "seurat", "res\\.", "RNA_snn", "wsnn"))
  if (!is.null(cand)) return(cand)
  obj$cluster_auto_from_idents <- as.character(Idents(obj))
  return("cluster_auto_from_idents")
}

choose_donor_col <- function(meta, user_col = NULL) {
  if (!is.null(user_col) && user_col %in% names(meta)) return(user_col)
  find_first_col(
    meta,
    c("^donor$", "donor_id", "individual", "subject", "patient", "sample", "orig.ident", "library", "case_id")
  )
}

choose_dx_col <- function(meta, user_col = NULL) {
  if (!is.null(user_col) && user_col %in% names(meta)) return(user_col)
  find_first_col(
    meta,
    c("^dx$", "^group$", "diagnosis", "condition", "status", "phenotype", "disease")
  )
}

get_fc_col <- function(markers) {
  fc_candidates <- c("avg_log2FC", "avg_logFC", "log2FC", "logFC")
  hit <- fc_candidates[fc_candidates %in% colnames(markers)]
  if (length(hit) == 0) return(NULL)
  hit[1]
}

calc_wilcox_p <- function(df, value_col, group_col) {
  df <- as.data.frame(df)
  g <- unique(df[[group_col]])
  g <- g[!is.na(g)]
  if (length(g) != 2) return(NA_real_)

  x1 <- df[[value_col]][df[[group_col]] == g[1]]
  x2 <- df[[value_col]][df[[group_col]] == g[2]]

  x1 <- x1[!is.na(x1)]
  x2 <- x2[!is.na(x2)]

  if (length(x1) < 1 || length(x2) < 1) return(NA_real_)
  suppressWarnings(wilcox.test(x1, x2, exact = FALSE)$p.value)
}

mean_diff_asd_minus_ctrl <- function(df, value_col, group_col) {
  df <- as.data.frame(df)
  grp <- as.character(df[[group_col]])
  val <- df[[value_col]]
  if (!all(c("ASD", "Control") %in% unique(grp))) return(NA_real_)
  mean(val[grp == "ASD"], na.rm = TRUE) - mean(val[grp == "Control"], na.rm = TRUE)
}

safe_scale <- function(x) {
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  if (sd(x, na.rm = TRUE) == 0) return(rep(0, length(x)))
  as.numeric(scale(x))
}

row_zscore_by_gene <- function(mat) {
  if (is.null(dim(mat))) return(mat)
  z <- t(apply(mat, 1, function(v) {
    if (all(is.na(v))) return(rep(NA_real_, length(v)))
    if (sd(v, na.rm = TRUE) == 0) return(rep(0, length(v)))
    as.numeric(scale(v))
  }))
  rownames(z) <- rownames(mat)
  colnames(z) <- colnames(mat)
  z
}

# curated gene sets
homeostatic_genes <- c(
  "P2RY12", "CX3CR1", "TMEM119", "CSF1R", "AIF1",
  "C1QA", "C1QB", "C1QC", "MERTK", "TREM2", "TYROBP"
)

activation_phago_genes <- c(
  "SPP1", "GAS6", "APOE", "LPL", "GPNMB", "FABP5",
  "CD68", "CTSB", "CTSD", "C3", "CD74", "LGALS3",
  "ABCA1", "MSR1", "LIPA", "LILRB4", "FCER1G"
)

anchor_genes <- c(
  "SPP1", "GAS6", "APOE", "MERTK", "CD74", "C3",
  "P2RY12", "CX3CR1", "TMEM119", "LPL", "GPNMB", "FABP5"
)

# --------------------------------------------------
# 1. read object
# --------------------------------------------------
msg("Reading microglia object ...")
obj <- readRDS(cfg$microglia_rds)
if (!inherits(obj, "Seurat")) stop("Input RDS is not a Seurat object.")

if (cfg$assay %in% Assays(obj)) {
  DefaultAssay(obj) <- cfg$assay
}

cluster_col <- choose_cluster_col(obj, cfg$cluster_col)
obj$cluster_use <- as.character(obj@meta.data[[cluster_col]])
Idents(obj) <- obj$cluster_use

meta <- obj@meta.data
donor_col <- choose_donor_col(meta, cfg$donor_col)
dx_col    <- choose_dx_col(meta, cfg$dx_col)

if (is.null(donor_col)) stop("Could not identify donor column.")
if (is.null(dx_col)) stop("Could not identify diagnosis column.")

cluster_levels <- unique(as.character(obj$cluster_use))
num_try <- suppressWarnings(as.numeric(cluster_levels))
if (all(!is.na(num_try))) {
  cluster_levels <- as.character(sort(num_try))
} else {
  cluster_levels <- sort(cluster_levels)
}

audit <- data.frame(
  item  = c("n_cells", "n_features", "default_assay", "cluster_col", "donor_col", "dx_col"),
  value = c(ncol(obj), nrow(obj), DefaultAssay(obj), cluster_col, donor_col, dx_col)
)
fwrite(audit, file.path(cfg$outdir, "Package1_Audit_Metadata.tsv"), sep = "\t")

# --------------------------------------------------
# 2. UMAP
# --------------------------------------------------
if ("umap" %in% names(obj@reductions)) {
  msg("Saving all-cluster UMAP ...")
  obj$cluster_use <- factor(obj$cluster_use, levels = cluster_levels)

  pal <- setNames(grDevices::hcl.colors(length(cluster_levels), "Dynamic"), cluster_levels)

  p_umap <- DimPlot(
    obj,
    reduction = "umap",
    group.by = "cluster_use",
    label = TRUE,
    repel = TRUE,
    raster = FALSE
  ) +
    scale_color_manual(values = pal) +
    theme_classic(base_size = 12) +
    labs(title = "PsychENCODE microglia: all clusters (discovery audit)") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  ggsave(file.path(cfg$outdir, "01_microglia_allcluster_umap.pdf"), p_umap, width = 9, height = 7)
  ggsave(file.path(cfg$outdir, "01_microglia_allcluster_umap.png"), p_umap, width = 9, height = 7, dpi = 300)

  if (donor_col %in% colnames(obj@meta.data)) {
    p_umap_donor <- DimPlot(
      obj,
      reduction = "umap",
      group.by = donor_col,
      raster = FALSE
    ) +
      theme_classic(base_size = 10) +
      labs(title = "PsychENCODE microglia UMAP colored by donor") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))

    ggsave(file.path(cfg$outdir, "02_microglia_umap_by_donor.pdf"), p_umap_donor, width = 12, height = 9)
    ggsave(file.path(cfg$outdir, "02_microglia_umap_by_donor.png"), p_umap_donor, width = 12, height = 9, dpi = 300)
  }
}

# --------------------------------------------------
# 3. all-cluster markers
# --------------------------------------------------
msg("Running FindAllMarkers ...")
markers <- FindAllMarkers(
  object = obj,
  only.pos = TRUE,
  min.pct = cfg$min_pct,
  logfc.threshold = cfg$logfc
)

if (!"gene" %in% colnames(markers)) {
  markers$gene <- rownames(markers)
}

fc_col <- get_fc_col(markers)
if (is.null(fc_col)) stop("Could not identify fold-change column in marker table.")

markers <- markers %>%
  arrange(cluster, desc(.data[[fc_col]]), p_val_adj)

fwrite(markers, file.path(cfg$outdir, "03_allcluster_markers.tsv.gz"), sep = "\t")

top10_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = .data[[fc_col]], n = 10, with_ties = FALSE) %>%
  ungroup()

fwrite(top10_markers, file.path(cfg$outdir, "04_allcluster_markers_top10.tsv.gz"), sep = "\t")

# signature files
sig_audit <- list()
for (cl in cluster_levels) {
  subm <- markers %>%
    filter(cluster == cl) %>%
    arrange(desc(.data[[fc_col]]), p_val_adj)

  genes <- head(unique(subm$gene), cfg$top_n)
  writeLines(
    genes,
    file.path(cfg$outdir, "Cluster_Signatures", paste0("Cluster", cl, "_top", cfg$top_n, ".txt"))
  )

  sig_audit[[length(sig_audit) + 1]] <- data.frame(
    cluster = cl,
    n_signature_genes = length(genes),
    top_marker_preview = paste(head(genes, 12), collapse = ", "),
    stringsAsFactors = FALSE
  )
}
sig_audit_df <- rbindlist(sig_audit, fill = TRUE)
fwrite(sig_audit_df, file.path(cfg$outdir, "05_cluster_signature_preview.tsv"), sep = "\t")

# --------------------------------------------------
# 4. donor-level proportions
# --------------------------------------------------
msg("Calculating donor-level proportions ...")
meta2 <- obj@meta.data
meta2$cell_id <- colnames(obj)
meta2$donor   <- as.character(meta2[[donor_col]])
meta2$dx      <- normalize_dx(meta2[[dx_col]])
meta2$cluster <- as.character(meta2$cluster_use)

prop_df <- meta2 %>%
  count(donor, dx, cluster, name = "n_cluster") %>%
  left_join(
    meta2 %>% count(donor, dx, name = "n_total"),
    by = c("donor", "dx")
  ) %>%
  mutate(prop = n_cluster / n_total) %>%
  arrange(cluster, dx, donor)

fwrite(prop_df, file.path(cfg$outdir, "06_allcluster_donor_proportions.tsv.gz"), sep = "\t")

abund_stats <- lapply(cluster_levels, function(cl) {
  sub <- prop_df %>% filter(cluster == cl)
  data.frame(
    cluster = cl,
    n_donors = length(unique(sub$donor)),
    mean_prop_ASD = mean(sub$prop[sub$dx == "ASD"], na.rm = TRUE),
    mean_prop_Control = mean(sub$prop[sub$dx == "Control"], na.rm = TRUE),
    effect_ASD_minus_Control = mean_diff_asd_minus_ctrl(sub, "prop", "dx"),
    wilcox_p = calc_wilcox_p(sub, "prop", "dx"),
    stringsAsFactors = FALSE
  )
})
abund_stats_df <- rbindlist(abund_stats, fill = TRUE)
abund_stats_df$wilcox_fdr <- p.adjust(abund_stats_df$wilcox_p, method = "BH")
fwrite(abund_stats_df, file.path(cfg$outdir, "07_allcluster_abundance_stats.tsv"), sep = "\t")

# donor-level proportion plot
prop_plot_df <- prop_df
prop_plot_df$cluster <- factor(prop_plot_df$cluster, levels = cluster_levels)

p_prop_all <- ggplot(prop_plot_df, aes(x = dx, y = prop, fill = dx)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.75) +
  geom_jitter(width = 0.12, size = 1.3, alpha = 0.8) +
  facet_wrap(~ cluster, scales = "free_y") +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 20, hjust = 1)
  ) +
  labs(
    title = "Discovery donor-level proportions across all microglial clusters",
    x = NULL,
    y = "Proportion within donor microglia"
  )

ggsave(file.path(cfg$outdir, "08_allcluster_donor_proportions_boxplot.pdf"), p_prop_all, width = 12, height = 8)
ggsave(file.path(cfg$outdir, "08_allcluster_donor_proportions_boxplot.png"), p_prop_all, width = 12, height = 8, dpi = 300)

# --------------------------------------------------
# 5. curated program scores
# --------------------------------------------------
msg("Computing curated program scores across clusters ...")
avg_list <- AverageExpression(
  obj,
  assays = DefaultAssay(obj),
  group.by = "cluster_use",
  slot = "data",
  verbose = FALSE
)

avg_mat <- avg_list[[DefaultAssay(obj)]]
avg_mat <- avg_mat[, cluster_levels, drop = FALSE]

get_program_score <- function(mat, genes) {
  genes_use <- intersect(genes, rownames(mat))
  if (length(genes_use) == 0) {
    return(list(score = setNames(rep(NA_real_, ncol(mat)), colnames(mat)), n_genes = 0, genes_used = character(0)))
  }
  sub <- mat[genes_use, , drop = FALSE]
  zsub <- row_zscore_by_gene(sub)
  sc <- colMeans(zsub, na.rm = TRUE)
  list(score = sc, n_genes = length(genes_use), genes_used = genes_use)
}

homeo_sc <- get_program_score(avg_mat, homeostatic_genes)
act_sc   <- get_program_score(avg_mat, activation_phago_genes)

expr_prog_df <- data.frame(
  cluster = cluster_levels,
  homeostatic_score_z = as.numeric(homeo_sc$score[cluster_levels]),
  activation_phago_score_z = as.numeric(act_sc$score[cluster_levels]),
  stringsAsFactors = FALSE
)
fwrite(expr_prog_df, file.path(cfg$outdir, "09_curated_program_scores.tsv"), sep = "\t")

writeLines(homeo_sc$genes_used, file.path(cfg$outdir, "GenesUsed_homeostatic.txt"))
writeLines(act_sc$genes_used, file.path(cfg$outdir, "GenesUsed_activation_phago.txt"))

# --------------------------------------------------
# 6. marker-based audit
# --------------------------------------------------
msg("Summarizing marker support per cluster ...")
marker_audit <- lapply(cluster_levels, function(cl) {
  subm <- markers %>%
    filter(cluster == cl) %>%
    arrange(desc(.data[[fc_col]]), p_val_adj)

  top50 <- head(unique(subm$gene), cfg$top_n)

  data.frame(
    cluster = cl,
    n_homeostatic_markers_in_top50 = sum(top50 %in% homeostatic_genes),
    n_activation_markers_in_top50 = sum(top50 %in% activation_phago_genes),
    flag_SPP1_top50 = as.integer("SPP1" %in% top50),
    flag_GAS6_top50 = as.integer("GAS6" %in% top50),
    flag_APOE_top50 = as.integer("APOE" %in% top50),
    flag_P2RY12_top50 = as.integer("P2RY12" %in% top50),
    flag_CX3CR1_top50 = as.integer("CX3CR1" %in% top50),
    flag_TMEM119_top50 = as.integer("TMEM119" %in% top50),
    top50_preview = paste(head(top50, 12), collapse = ", "),
    stringsAsFactors = FALSE
  )
})
marker_audit_df <- rbindlist(marker_audit, fill = TRUE)
fwrite(marker_audit_df, file.path(cfg$outdir, "10_marker_support.tsv"), sep = "\t")

# --------------------------------------------------
# 7. prioritization table
# --------------------------------------------------
msg("Building unbiased prioritization table ...")
prior_df <- abund_stats_df %>%
  left_join(expr_prog_df, by = "cluster") %>%
  left_join(marker_audit_df, by = "cluster") %>%
  left_join(sig_audit_df, by = "cluster")

prior_df$abundance_effect_pos <- pmax(prior_df$effect_ASD_minus_Control, 0)
prior_df$abundance_effect_pos_z <- safe_scale(prior_df$abundance_effect_pos)
prior_df$abundance_p_support <- ifelse(prior_df$effect_ASD_minus_Control > 0, -log10(pmax(prior_df$wilcox_p, 1e-300)), 0)
prior_df$abundance_p_support_z <- safe_scale(prior_df$abundance_p_support)
prior_df$homeostatic_score_z_scaled <- safe_scale(prior_df$homeostatic_score_z)
prior_df$activation_phago_score_z_scaled <- safe_scale(prior_df$activation_phago_score_z)

prior_df$marker_flag_sum <- with(
  prior_df,
  n_activation_markers_in_top50 +
    flag_SPP1_top50 + flag_GAS6_top50 + flag_APOE_top50 +
    flag_P2RY12_top50 + flag_CX3CR1_top50 + flag_TMEM119_top50
)
prior_df$marker_flag_sum_z <- safe_scale(prior_df$marker_flag_sum)

prior_df$heuristic_priority_score <- with(
  prior_df,
  1.0 * abundance_effect_pos_z +
    0.8 * abundance_p_support_z +
    1.0 * activation_phago_score_z_scaled +
    0.5 * homeostatic_score_z_scaled +
    0.6 * marker_flag_sum_z
)

prior_df <- prior_df %>%
  arrange(desc(heuristic_priority_score), desc(effect_ASD_minus_Control), wilcox_p)

prior_df$heuristic_rank <- seq_len(nrow(prior_df))
fwrite(prior_df, file.path(cfg$outdir, "11_cluster_prioritization.tsv"), sep = "\t")

top_hit <- prior_df %>% slice(1)
write.table(top_hit, file.path(cfg$outdir, "12_top_candidate.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cluster2_row <- prior_df %>% filter(cluster == "2")
if (nrow(cluster2_row) > 0) {
  write.table(cluster2_row, file.path(cfg$outdir, "13_cluster2_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}

# --------------------------------------------------
# 8. plots
# --------------------------------------------------
msg("Saving audit plots ...")

plot_df <- prior_df
plot_df$cluster_label <- paste0("C", plot_df$cluster)
plot_df$is_cluster2 <- ifelse(plot_df$cluster == "2", "Cluster2", "Others")
plot_df$neglog10p <- -log10(pmax(plot_df$wilcox_p, 1e-300))

p_scatter <- ggplot(plot_df, aes(x = effect_ASD_minus_Control, y = activation_phago_score_z, size = neglog10p, shape = is_cluster2)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
  geom_hline(yintercept = 0, linetype = 2, color = "grey60") +
  geom_point() +
  geom_text(aes(label = cluster_label), vjust = -0.8, size = 3, check_overlap = TRUE) +
  theme_classic(base_size = 12) +
  labs(
    title = "Unbiased audit of all microglial clusters",
    x = "Donor-level abundance effect (ASD - Control)",
    y = "Activation/phagocytic program score (z)",
    size = "-log10(p)",
    shape = NULL
  )

ggsave(file.path(cfg$outdir, "14_priority_scatter.pdf"), p_scatter, width = 7.2, height = 6.0)
ggsave(file.path(cfg$outdir, "14_priority_scatter.png"), p_scatter, width = 7.2, height = 6.0, dpi = 300)

rank_plot_df <- prior_df %>%
  mutate(cluster_label = factor(paste0("C", cluster), levels = rev(paste0("C", cluster))))

p_rank <- ggplot(rank_plot_df, aes(x = heuristic_priority_score, y = cluster_label)) +
  geom_col() +
  theme_classic(base_size = 12) +
  labs(
    title = "Heuristic cluster prioritization (audit only)",
    x = "Heuristic priority score",
    y = "Cluster"
  )

ggsave(file.path(cfg$outdir, "15_heuristic_ranking.pdf"), p_rank, width = 6.5, height = 6.5)
ggsave(file.path(cfg$outdir, "15_heuristic_ranking.png"), p_rank, width = 6.5, height = 6.5, dpi = 300)

# --------------------------------------------------
# 9. dotplots / featureplots
# --------------------------------------------------
obj$cluster_use <- factor(obj$cluster_use, levels = cluster_levels)

homeo_dot_genes <- intersect(c("P2RY12", "CX3CR1", "TMEM119", "CSF1R", "AIF1"), rownames(obj))
if (length(homeo_dot_genes) > 0) {
  p_homeo_dot <- DotPlot(obj, features = homeo_dot_genes, group.by = "cluster_use") +
    RotatedAxis() +
    theme_classic(base_size = 12) +
    labs(title = "Homeostatic microglial markers across clusters")
  ggsave(file.path(cfg$outdir, "16_homeostatic_marker_dotplot.pdf"), p_homeo_dot, width = 7.5, height = 4.2)
  ggsave(file.path(cfg$outdir, "16_homeostatic_marker_dotplot.png"), p_homeo_dot, width = 7.5, height = 4.2, dpi = 300)
}

anchor_genes_use <- intersect(anchor_genes, rownames(obj))
if (length(anchor_genes_use) > 0) {
  p_dot <- DotPlot(obj, features = anchor_genes_use, group.by = "cluster_use") +
    RotatedAxis() +
    theme_classic(base_size = 12) +
    labs(title = "Curated candidate/homeostatic genes across microglial clusters")

  ggsave(file.path(cfg$outdir, "17_curated_gene_dotplot.pdf"), p_dot, width = 10, height = 5.2)
  ggsave(file.path(cfg$outdir, "17_curated_gene_dotplot.png"), p_dot, width = 10, height = 5.2, dpi = 300)

  if ("umap" %in% names(obj@reductions)) {
    fp_list <- lapply(anchor_genes_use, function(g) {
      FeaturePlot(obj, features = g, reduction = "umap", raster = FALSE) +
        theme_classic(base_size = 10) +
        ggtitle(g)
    })
    p_fp <- wrap_plots(fp_list, ncol = 3)
    ggsave(file.path(cfg$outdir, "18_curated_gene_featureplots_umap.pdf"), p_fp, width = 12, height = 14)
    ggsave(file.path(cfg$outdir, "18_curated_gene_featureplots_umap.png"), p_fp, width = 12, height = 14, dpi = 300)
  }
}

saveRDS(obj, file.path(cfg$outdir, "19_microglia_object_with_audit.rds"))
save_session_info(file.path(cfg$outdir, "sessionInfo_package1.txt"))

manifest <- data.frame(
  output_file = c(
    "01_microglia_allcluster_umap",
    "02_microglia_umap_by_donor",
    "03_allcluster_markers",
    "04_allcluster_markers_top10",
    "06_allcluster_donor_proportions",
    "07_allcluster_abundance_stats",
    "08_allcluster_donor_proportions_boxplot",
    "09_curated_program_scores",
    "10_marker_support",
    "11_cluster_prioritization",
    "13_cluster2_summary",
    "16_homeostatic_marker_dotplot",
    "17_curated_gene_dotplot",
    "18_curated_gene_featureplots_umap"
  ),
  description = c(
    "all-cluster UMAP with unique colors",
    "UMAP colored by donor",
    "all microglia cluster markers",
    "top10 markers per cluster",
    "donor-level cluster proportions",
    "donor-level abundance statistics",
    "boxplot of donor-level cluster proportions",
    "homeostatic and activation/phagocytic program scores",
    "marker support summary",
    "heuristic prioritization of clusters",
    "Cluster 2 summary table",
    "homeostatic marker dotplot",
    "candidate/homeostatic curated gene dotplot",
    "UMAP feature plots for curated genes"
  ),
  stringsAsFactors = FALSE
)
fwrite(manifest, file.path(cfg$outdir, "PACKAGE1_manifest.tsv"), sep = "\t")

msg("Package 1 discovery audit finished.")
