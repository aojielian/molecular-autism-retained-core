#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(Matrix)
})

options(stringsAsFactors = FALSE)

# =========================================================
# Figure 1 enhanced script for Molecular Autism revision
# 目标：
# 1) 输出 all-cluster unique-color UMAP
# 2) 输出 all-cluster markers
# 3) 输出 all-cluster donor-level proportions
# 4) 保留一个 cluster 2 donor-level proportion 主图
# =========================================================

cfg <- list(
  microglia_rds = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/PsychENCODE_Science_2024/Step3_Microglia/Microglia_Clustered.rds",
  outdir       = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Figure1_Enhanced",
  assay        = "RNA",
  cluster_col  = NULL,   # 为空时自动识别
  donor_col    = NULL,   # 为空时自动识别
  dx_col       = NULL,   # 为空时自动识别
  target_cluster = "2",
  min_pct      = 0.10,
  logfc        = 0.25
)

dir.create(cfg$outdir, recursive = TRUE, showWarnings = FALSE)

msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ..., "\n", sep = "")
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x

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
  if (!is.null(user_col) && user_col %in% names(meta)) {
    return(user_col)
  }
  if ("seurat_clusters" %in% names(meta)) return("seurat_clusters")

  cand <- find_first_col(
    meta,
    c("^cluster$", "cluster", "seurat", "res\\.", "RNA_snn", "wsnn", "integrated_snn")
  )
  if (!is.null(cand)) return(cand)

  # 回退到当前 Idents
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

ensure_reduction <- function(obj) {
  reds <- names(obj@reductions)
  if ("umap" %in% reds) return(list(obj = obj, reduction = "umap"))
  if ("tsne" %in% reds) return(list(obj = obj, reduction = "tsne"))

  msg("No UMAP/tSNE found. Recomputing PCA + UMAP ...")
  DefaultAssay(obj) <- if ("RNA" %in% Assays(obj)) "RNA" else DefaultAssay(obj)

  if (ncol(GetAssayData(obj, assay = DefaultAssay(obj), slot = "data")) == 0) {
    obj <- NormalizeData(obj, verbose = FALSE)
  }
  obj <- FindVariableFeatures(obj, nfeatures = min(3000, nrow(obj)), verbose = FALSE)
  obj <- ScaleData(obj, features = VariableFeatures(obj), verbose = FALSE)
  obj <- RunPCA(obj, features = VariableFeatures(obj), verbose = FALSE)
  npc <- min(30, ncol(Embeddings(obj, "pca")))
  obj <- RunUMAP(obj, dims = 1:npc, verbose = FALSE)
  list(obj = obj, reduction = "umap")
}

get_fc_col <- function(markers) {
  fc_candidates <- c("avg_log2FC", "avg_logFC", "log2FC", "logFC")
  fc_candidates[fc_candidates %in% colnames(markers)][1]
}

plot_p_from_wilcox <- function(df, value_col, group_col) {
  g <- unique(df[[group_col]])
  g <- g[!is.na(g)]
  if (length(g) != 2) return(NA_real_)
  x1 <- df[df[[group_col]] == g[1], value_col, drop = TRUE]
  x2 <- df[df[[group_col]] == g[2], value_col, drop = TRUE]
  if (length(x1) < 1 || length(x2) < 1) return(NA_real_)
  wilcox.test(x1, x2, exact = FALSE)$p.value
}

# -----------------------------
# 1. 读取对象
# -----------------------------
msg("Reading object: ", cfg$microglia_rds)
obj <- readRDS(cfg$microglia_rds)
if (!inherits(obj, "Seurat")) stop("Input RDS is not a Seurat object.")

if (cfg$assay %in% Assays(obj)) {
  DefaultAssay(obj) <- cfg$assay
} else {
  DefaultAssay(obj) <- Assays(obj)[1]
}

cluster_col <- choose_cluster_col(obj, cfg$cluster_col)
obj$cluster_use <- as.character(obj@meta.data[[cluster_col]])
Idents(obj) <- obj$cluster_use

meta <- obj@meta.data
donor_col <- choose_donor_col(meta, cfg$donor_col)
dx_col    <- choose_dx_col(meta, cfg$dx_col)

if (is.null(donor_col)) stop("Could not identify donor column. Please set cfg$donor_col manually.")
if (is.null(dx_col))    stop("Could not identify diagnosis column. Please set cfg$dx_col manually.")

red_res <- ensure_reduction(obj)
obj <- red_res$obj
reduction_use <- red_res$reduction

audit <- data.frame(
  item  = c("n_cells", "n_features", "default_assay", "cluster_col", "donor_col", "dx_col", "reduction"),
  value = c(ncol(obj), nrow(obj), DefaultAssay(obj), cluster_col, donor_col, dx_col, reduction_use)
)
fwrite(audit, file.path(cfg$outdir, "Figure1_Audit.tsv"), sep = "\t")

# -----------------------------
# 2. all-cluster UMAP
# -----------------------------
msg("Plotting all-cluster UMAP ...")
cluster_levels <- sort(unique(as.character(obj$cluster_use)))
obj$cluster_use <- factor(as.character(obj$cluster_use), levels = cluster_levels)
Idents(obj) <- obj$cluster_use

pal <- setNames(
  grDevices::hcl.colors(length(cluster_levels), palette = "Dynamic"),
  cluster_levels
)

p_umap <- DimPlot(
  obj,
  reduction = reduction_use,
  group.by  = "cluster_use",
  label     = TRUE,
  repel     = TRUE,
  raster    = FALSE
) +
  scale_color_manual(values = pal) +
  labs(title = "PsychENCODE microglia: all clusters") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "right"
  )

ggsave(
  filename = file.path(cfg$outdir, "Supp_MG_UMAP_AllClusters.pdf"),
  plot = p_umap, width = 9, height = 7
)

# -----------------------------
# 3. all-cluster markers
# -----------------------------
msg("Running FindAllMarkers ...")
markers <- FindAllMarkers(
  object = obj,
  only.pos = TRUE,
  min.pct = cfg$min_pct,
  logfc.threshold = cfg$logfc
)

if (nrow(markers) == 0) {
  warning("FindAllMarkers returned 0 rows.")
  markers <- data.frame()
} else {
  fc_col <- get_fc_col(markers)
  if (!is.null(fc_col)) {
    markers <- markers %>%
      arrange(cluster, desc(.data[[fc_col]]), p_val_adj)
  }
}

fwrite(markers, file.path(cfg$outdir, "Supp_AllCluster_Markers.csv"))

if (nrow(markers) > 0) {
  fc_col <- get_fc_col(markers)
  top10 <- markers %>%
    group_by(cluster) %>%
    slice_max(order_by = .data[[fc_col]], n = 10, with_ties = FALSE) %>%
    ungroup()
  fwrite(top10, file.path(cfg$outdir, "Supp_AllCluster_Markers_Top10.csv"))
}

# -----------------------------
# 4. donor-level proportions
# -----------------------------
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

fwrite(prop_df, file.path(cfg$outdir, "Supp_AllCluster_Proportions.csv"))

p_prop_all <- ggplot(prop_df, aes(x = dx, y = prop, fill = dx)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.7) +
  geom_jitter(width = 0.12, size = 1.6, alpha = 0.9) +
  facet_wrap(~ cluster, scales = "free_y") +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "none",
    axis.text.x = element_text(angle = 20, hjust = 1)
  ) +
  labs(
    title = "Donor-level proportions across all microglial clusters",
    x = NULL, y = "Proportion within donor microglia"
  )

ggsave(
  filename = file.path(cfg$outdir, "Supp_AllCluster_Proportions_Boxplot.pdf"),
  plot = p_prop_all, width = 12, height = 8
)

# -----------------------------
# 5. cluster 2 donor-level 主图
# -----------------------------
if (cfg$target_cluster %in% unique(prop_df$cluster)) {
  msg("Plotting target cluster donor-level proportion: cluster ", cfg$target_cluster)
  c2_df <- prop_df %>% filter(cluster == cfg$target_cluster)
  fwrite(c2_df, file.path(cfg$outdir, paste0("Fig1C_Cluster", cfg$target_cluster, "_Proportion_DonorLevel.csv")))

  pval_c2 <- plot_p_from_wilcox(c2_df, "prop", "dx")
  p_lab_c2 <- ifelse(is.na(pval_c2), "Wilcoxon p = NA", paste0("Wilcoxon p = ", signif(pval_c2, 3)))

  p_c2 <- ggplot(c2_df, aes(x = dx, y = prop, fill = dx)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.75) +
    geom_jitter(width = 0.12, size = 2, alpha = 0.9) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 20, hjust = 1)
    ) +
    labs(
      title = paste0("Cluster ", cfg$target_cluster, " donor-level proportion"),
      subtitle = p_lab_c2,
      x = NULL,
      y = "Proportion within donor microglia"
    )

  ggsave(
    filename = file.path(cfg$outdir, paste0("Fig1C_Cluster", cfg$target_cluster, "_Proportion_DonorLevel.pdf")),
    plot = p_c2, width = 5.5, height = 4.8
  )
} else {
  msg("Target cluster ", cfg$target_cluster, " not found in object. Skipped cluster-specific plot.")
}

save_session_info(file.path(cfg$outdir, "sessionInfo_Figure1.txt"))
msg("Figure 1 enhanced analysis finished.")
