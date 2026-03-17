#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(Matrix)
})

options(stringsAsFactors = FALSE)

cfg <- list(
  microglia_rds = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/PsychENCODE_Science_2024/Step3_Microglia/Microglia_Clustered.rds",
  outdir        = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Figure1_Unbiased_ClusterAudit",
  assay         = "RNA",
  cluster_col   = NULL,
  donor_col     = NULL,
  dx_col        = NULL,
  top_n         = 50,
  min_pct       = 0.10,
  logfc         = 0.25
)

dir.create(cfg$outdir, recursive = TRUE, showWarnings = FALSE)

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

make_signature_name <- function(cl) paste0("Cluster", cl)

# curated gene programs
core_microglia_genes <- c(
  "P2RY12", "CX3CR1", "TMEM119", "CSF1R", "AIF1",
  "C1QA", "C1QB", "C1QC", "MERTK", "TREM2", "TYROBP"
)

risk_lipid_phago_genes <- c(
  "SPP1", "GAS6", "APOE", "LPL", "GPNMB", "FABP5",
  "CD68", "CTSB", "CTSD", "C3", "CD74", "LGALS3",
  "ABCA1", "MSR1", "LIPA", "LILRB4", "FCER1G"
)

anchor_genes <- c("SPP1", "GAS6", "APOE", "P2RY12", "CX3CR1", "TMEM119", "MERTK", "CD74", "C3", "LPL", "GPNMB", "FABP5")

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
fwrite(audit, file.path(cfg$outdir, "Figure1_Unbiased_Audit.tsv"), sep = "\t")

# --------------------------------------------------
# 2. optional all-cluster UMAP
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
    labs(title = "PsychENCODE microglia: all clusters (unbiased audit)") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))

  ggsave(file.path(cfg$outdir, "Fig1_Audit_AllCluster_UMAP.pdf"), p_umap, width = 9, height = 7)
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

fwrite(markers, file.path(cfg$outdir, "Fig1_Audit_AllCluster_Markers.csv"))

top10_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = .data[[fc_col]], n = 10, with_ties = FALSE) %>%
  ungroup()

fwrite(top10_markers, file.path(cfg$outdir, "Fig1_Audit_AllCluster_Markers_Top10.csv"))

# signature files for each cluster
dir.create(file.path(cfg$outdir, "Cluster_Signatures"), showWarnings = FALSE)
sig_audit <- list()

for (cl in cluster_levels) {
  subm <- markers %>%
    filter(cluster == cl) %>%
    arrange(desc(.data[[fc_col]]), p_val_adj)

  genes <- head(unique(subm$gene), cfg$top_n)
  writeLines(genes, file.path(cfg$outdir, "Cluster_Signatures", paste0("Cluster", cl, "_top", cfg$top_n, ".txt")))

  sig_audit[[length(sig_audit) + 1]] <- data.frame(
    cluster = cl,
    n_signature_genes = length(genes),
    top_marker_preview = paste(head(genes, 12), collapse = ", "),
    stringsAsFactors = FALSE
  )
}
sig_audit_df <- rbindlist(sig_audit, fill = TRUE)
fwrite(sig_audit_df, file.path(cfg$outdir, "Fig1_Audit_ClusterSignature_Preview.tsv"), sep = "\t")

# --------------------------------------------------
# 4. donor-level proportions for all clusters
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

fwrite(prop_df, file.path(cfg$outdir, "Fig1_Audit_AllCluster_Proportions.csv"))

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
fwrite(abund_stats_df, file.path(cfg$outdir, "Fig1_Audit_AllCluster_AbundanceStats.tsv"), sep = "\t")

# --------------------------------------------------
# 5. program scores from AverageExpression
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

core_sc <- get_program_score(avg_mat, core_microglia_genes)
risk_sc <- get_program_score(avg_mat, risk_lipid_phago_genes)

expr_prog_df <- data.frame(
  cluster = cluster_levels,
  core_microglia_score_z = as.numeric(core_sc$score[cluster_levels]),
  risk_lipid_phago_score_z = as.numeric(risk_sc$score[cluster_levels]),
  stringsAsFactors = FALSE
)
fwrite(expr_prog_df, file.path(cfg$outdir, "Fig1_Audit_CuratedProgramScores.tsv"), sep = "\t")

writeLines(core_sc$genes_used, file.path(cfg$outdir, "GenesUsed_core_microglia.txt"))
writeLines(risk_sc$genes_used, file.path(cfg$outdir, "GenesUsed_risk_lipid_phago.txt"))

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
    n_core_microglia_markers_in_top50 = sum(top50 %in% core_microglia_genes),
    n_risk_lipid_phago_markers_in_top50 = sum(top50 %in% risk_lipid_phago_genes),
    flag_SPP1_top50 = as.integer("SPP1" %in% top50),
    flag_GAS6_top50 = as.integer("GAS6" %in% top50),
    flag_APOE_top50 = as.integer("APOE" %in% top50),
    flag_P2RY12_top50 = as.integer("P2RY12" %in% top50),
    flag_CX3CR1_top50 = as.integer("CX3CR1" %in% top50),
    flag_TMEm119_top50 = as.integer("TMEM119" %in% top50),
    top50_preview = paste(head(top50, 12), collapse = ", "),
    stringsAsFactors = FALSE
  )
})
marker_audit_df <- rbindlist(marker_audit, fill = TRUE)
fwrite(marker_audit_df, file.path(cfg$outdir, "Fig1_Audit_MarkerSupport.tsv"), sep = "\t")

# --------------------------------------------------
# 7. merge into prioritization table
# --------------------------------------------------
msg("Building unbiased prioritization table ...")
prior_df <- abund_stats_df %>%
  left_join(expr_prog_df, by = "cluster") %>%
  left_join(marker_audit_df, by = "cluster") %>%
  left_join(sig_audit_df, by = "cluster")

# heuristic score for audit only
prior_df$abundance_effect_pos <- pmax(prior_df$effect_ASD_minus_Control, 0)
prior_df$abundance_effect_pos_z <- safe_scale(prior_df$abundance_effect_pos)
prior_df$abundance_p_support <- ifelse(prior_df$effect_ASD_minus_Control > 0, -log10(pmax(prior_df$wilcox_p, 1e-300)), 0)
prior_df$abundance_p_support_z <- safe_scale(prior_df$abundance_p_support)
prior_df$core_microglia_score_z_scaled <- safe_scale(prior_df$core_microglia_score_z)
prior_df$risk_lipid_phago_score_z_scaled <- safe_scale(prior_df$risk_lipid_phago_score_z)

prior_df$marker_flag_sum <- with(
  prior_df,
  n_risk_lipid_phago_markers_in_top50 +
    flag_SPP1_top50 + flag_GAS6_top50 + flag_APOE_top50 +
    flag_P2RY12_top50 + flag_CX3CR1_top50 + flag_TMEm119_top50
)
prior_df$marker_flag_sum_z <- safe_scale(prior_df$marker_flag_sum)

prior_df$heuristic_priority_score <- with(
  prior_df,
  1.0 * abundance_effect_pos_z +
    0.8 * abundance_p_support_z +
    1.0 * risk_lipid_phago_score_z_scaled +
    0.5 * core_microglia_score_z_scaled +
    0.6 * marker_flag_sum_z
)

prior_df <- prior_df %>%
  arrange(desc(heuristic_priority_score), desc(effect_ASD_minus_Control), wilcox_p)

prior_df$heuristic_rank <- seq_len(nrow(prior_df))
fwrite(prior_df, file.path(cfg$outdir, "Fig1_Audit_ClusterPrioritization.tsv"), sep = "\t")

top_hit <- prior_df %>% slice(1)
write.table(top_hit, file.path(cfg$outdir, "Fig1_Audit_TopCandidate.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

cluster2_row <- prior_df %>% filter(cluster == "2")
if (nrow(cluster2_row) > 0) {
  write.table(cluster2_row, file.path(cfg$outdir, "Fig1_Audit_Cluster2_Summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}

# --------------------------------------------------
# 8. plotting
# --------------------------------------------------
msg("Saving audit plots ...")

plot_df <- prior_df
plot_df$cluster_label <- paste0("C", plot_df$cluster)
plot_df$is_cluster2 <- ifelse(plot_df$cluster == "2", "Cluster2", "Others")
plot_df$neglog10p <- -log10(pmax(plot_df$wilcox_p, 1e-300))

# abundance vs program scatter
p_scatter <- ggplot(plot_df, aes(x = effect_ASD_minus_Control, y = risk_lipid_phago_score_z, size = neglog10p, shape = is_cluster2)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
  geom_hline(yintercept = 0, linetype = 2, color = "grey60") +
  geom_point() +
  geom_text(aes(label = cluster_label), vjust = -0.8, size = 3, check_overlap = TRUE) +
  theme_classic(base_size = 12) +
  labs(
    title = "Unbiased audit of all microglial clusters",
    x = "Donor-level abundance effect (ASD - Control)",
    y = "Risk/lipid/phagocytic program score (z)",
    size = "-log10(p)",
    shape = NULL
  )

ggsave(file.path(cfg$outdir, "Fig1_Audit_PriorityScatter.pdf"), p_scatter, width = 7.2, height = 6.0)

# heuristic ranking barplot
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

ggsave(file.path(cfg$outdir, "Fig1_Audit_HeuristicRanking.pdf"), p_rank, width = 6.5, height = 6.5)

# all-cluster donor proportion boxplots
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

ggsave(file.path(cfg$outdir, "Fig1_Audit_AllCluster_Proportions_Boxplot.pdf"), p_prop_all, width = 12, height = 8)

# curated gene dotplot
genes_dot <- intersect(anchor_genes, rownames(obj))
if (length(genes_dot) > 0) {
  obj$cluster_use <- factor(obj$cluster_use, levels = cluster_levels)
  p_dot <- DotPlot(obj, features = genes_dot, group.by = "cluster_use") +
    RotatedAxis() +
    theme_classic(base_size = 12) +
    labs(title = "Curated marker/program genes across microglial clusters")

  ggsave(file.path(cfg$outdir, "Fig1_Audit_CuratedGenes_DotPlot.pdf"), p_dot, width = 9.5, height = 5.2)
}

save_session_info(file.path(cfg$outdir, "sessionInfo_Figure1_Unbiased_ClusterAudit.txt"))
msg("Figure 1 unbiased cluster audit finished.")
