#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(edgeR)
  library(limma)
  library(optparse)
})

options(stringsAsFactors = FALSE)

DEFAULT_PROJECT_DIR <- normalizePath(Sys.getenv("MA_PROJECT_DIR", "."), mustWork = FALSE)

option_list <- list(
  make_option("--strict_rds", type = "character", default = Sys.getenv("MA_STRICT_RDS", NA_character_),
              help = "Path to strict retained-core Seurat object (.rds) [required or set MA_STRICT_RDS]"),
  make_option("--outdir", type = "character",
              default = file.path(DEFAULT_PROJECT_DIR, "results", "Package2b_ResidualBackgroundAudit")),
  make_option("--assay", type = "character", default = "RNA"),
  make_option("--cluster_col", type = "character", default = "cluster_use"),
  make_option("--donor_col", type = "character", default = "individual_ID"),
  make_option("--dx_col", type = "character", default = "Diagnosis"),
  make_option("--ref_cluster", type = "character", default = "0"),
  make_option("--cand_cluster", type = "character", default = "2"),
  make_option("--min_cells_per_pb_sample", type = "integer", default = 20),
  make_option("--neuronal_high_quantile", type = "double", default = 0.90),
  make_option("--scrublet_high_quantile", type = "double", default = 0.95)
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.na(opt$strict_rds) || opt$strict_rds == "") stop("Please provide --strict_rds or set MA_STRICT_RDS.")

cfg <- list(
  strict_rds = opt$strict_rds,
  outdir = opt$outdir,
  assay = opt$assay,
  cluster_col = opt$cluster_col,
  donor_col = opt$donor_col,
  dx_col = opt$dx_col,
  ref_cluster = opt$ref_cluster,
  cand_cluster = opt$cand_cluster,
  min_cells_per_pb_sample = opt$min_cells_per_pb_sample,
  neuronal_high_quantile = opt$neuronal_high_quantile,
  scrublet_high_quantile = opt$scrublet_high_quantile,
  qc_cols_try = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent_mito", "single_scrublet_score"),
  microglia_homeostatic = c("P2RY12","CX3CR1","TMEM119","CSF1R","AIF1","C1QA","C1QB","C1QC","TREM2","TYROBP","MERTK","FCER1G"),
  microglia_activation = c("SPP1","GAS6","APOE","CD74","C3","LGALS3","GPNMB","LPL","FABP5","CTSB","CTSD","ABCA1","MSR1"),
  neuronal_synaptic = c("NRG3","RBFOX1","LRRTM4","KCNIP4","LRP1B","CNTNAP2","CSMD1","SYT1","ROBO2","CADM2","GPM6A","LSAMP","DPP10","NRXN3","DLGAP1","PCDH9","DLG2","SNAP25","STMN2","RBFOX3"),
  astrocyte = c("GFAP","AQP4","ALDH1L1","SLC1A2","SLC1A3","GJA1","ALDOC"),
  oligodendrocyte = c("MBP","MOG","PLP1","MOBP","MAG","CNP"),
  opc = c("PDGFRA","VCAN","PTPRZ1","CSPG4","OLIG1","OLIG2"),
  endothelial = c("CLDN5","FLT1","KDR","VWF","EPAS1","ESAM","RAMP2"),
  lymphocyte = c("CD3D","CD3E","CD247","TRBC1","TRBC2","THEMIS","MS4A1","CD79A","CD79B","IGHM","PAX5"),
  neuronal_background_genes = c("NRG3","RBFOX1","LRRTM4","KCNIP4","LRP1B","CNTNAP2","CSMD1","SYT1","ROBO2","CADM2","GPM6A","LSAMP","DPP10","NRXN3","DLGAP1","PCDH9","DLG2"),
  audit_genes = c("SPP1","GAS6","APOE","MERTK","CD74","C3","P2RY12","CX3CR1","TMEM119")
)

dir.create(cfg$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cfg$outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cfg$outdir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cfg$outdir, "objects"), recursive = TRUE, showWarnings = FALSE)

msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ..., "\n", sep = "")
}

save_session_info <- function(path) {
  writeLines(capture.output(sessionInfo()), con = path)
}

normalize_dx <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  y <- x
  y[grepl("asd|autism|case", x, ignore.case = TRUE)] <- "ASD"
  y[grepl("control|ctrl|ctl|healthy|non[- ]?psy", x, ignore.case = TRUE)] <- "Control"
  y
}

norm_cluster_label <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- sub("^Cluster", "", x, ignore.case = TRUE)
  x <- sub("^cl", "", x, ignore.case = TRUE)
  x <- sub("\\.0$", "", x)
  x
}

row_zscore_dense <- function(mat) {
  z <- t(apply(mat, 1, function(v) {
    if (all(is.na(v))) return(rep(NA_real_, length(v)))
    s <- sd(v, na.rm = TRUE)
    if (is.na(s) || s == 0) return(rep(0, length(v)))
    as.numeric(scale(v))
  }))
  rownames(z) <- rownames(mat)
  colnames(z) <- colnames(mat)
  z
}

compute_module_score <- function(expr_mat, genes) {
  genes_use <- intersect(genes, rownames(expr_mat))
  if (length(genes_use) == 0) {
    return(list(score = rep(NA_real_, ncol(expr_mat)), genes_used = character(0)))
  }
  sub <- as.matrix(expr_mat[genes_use, , drop = FALSE])
  zsub <- row_zscore_dense(sub)
  sc <- colMeans(zsub, na.rm = TRUE)
  list(score = sc, genes_used = genes_use)
}

calc_wilcox <- function(df, value_col, group_col) {
  dd <- as.data.frame(df)
  g <- unique(dd[[group_col]])
  g <- g[!is.na(g)]
  if (length(g) != 2) return(NA_real_)
  x1 <- dd[[value_col]][dd[[group_col]] == g[1]]
  x2 <- dd[[value_col]][dd[[group_col]] == g[2]]
  x1 <- x1[!is.na(x1)]
  x2 <- x2[!is.na(x2)]
  if (length(x1) < 1 || length(x2) < 1) return(NA_real_)
  suppressWarnings(wilcox.test(x1, x2, exact = FALSE)$p.value)
}

paired_wilcox_by_donor <- function(df, value_col) {
  dd <- as.data.frame(df)
  need <- dd %>% count(donor) %>% filter(n == 2) %>% pull(donor)
  dd <- dd %>% filter(donor %in% need)
  if (nrow(dd) == 0) return(NA_real_)
  wide <- reshape(dd[, c("donor","cluster", value_col)],
                  idvar = "donor", timevar = "cluster", direction = "wide")
  c0 <- paste0(value_col, ".", cfg$ref_cluster)
  c2 <- paste0(value_col, ".", cfg$cand_cluster)
  if (!all(c(c0, c2) %in% colnames(wide))) return(NA_real_)
  x0 <- wide[[c0]]
  x2 <- wide[[c2]]
  keep <- !is.na(x0) & !is.na(x2)
  if (sum(keep) < 2) return(NA_real_)
  suppressWarnings(wilcox.test(x2[keep], x0[keep], paired = TRUE, exact = FALSE)$p.value)
}

sparse_pseudobulk <- function(count_mat, groups) {
  groups <- factor(groups)
  mm <- sparse.model.matrix(~ 0 + groups)
  colnames(mm) <- sub("^groups", "", colnames(mm))
  pb <- count_mat %*% mm
  pb
}

run_paired_pseudobulk <- function(obj_sub, donor_col, dx_col, cluster_col, out_prefix, outdir, min_cells = 20) {
  md <- obj_sub@meta.data
  md$cell_id <- colnames(obj_sub)
  md$cluster <- norm_cluster_label(md[[cluster_col]])
  md$donor <- as.character(md[[donor_col]])
  md$dx <- normalize_dx(md[[dx_col]])

  count_per_sample <- md %>%
    count(donor, cluster, dx, name = "n_cells") %>%
    filter(cluster %in% c(cfg$ref_cluster, cfg$cand_cluster),
           n_cells >= min_cells)

  paired_donors <- count_per_sample %>%
    count(donor, name = "n_clusters_kept") %>%
    filter(n_clusters_kept == 2) %>%
    pull(donor) %>% as.character()

  pb_meta <- count_per_sample %>%
    filter(donor %in% paired_donors) %>%
    mutate(
      sample_id = paste0(donor, "__cl", cluster),
      cluster_simple = ifelse(cluster == cfg$ref_cluster, "ref", "cand")
    ) %>%
    arrange(donor, cluster)

  if (nrow(pb_meta) == 0) {
    writeLines("No paired pseudobulk samples available after filtering.",
               file.path(outdir, paste0(out_prefix, "__SKIPPED.txt")))
    return(NULL)
  }

  pb_cells <- md %>%
    filter(donor %in% paired_donors, cluster %in% c(cfg$ref_cluster, cfg$cand_cluster)) %>%
    semi_join(pb_meta %>% select(donor, cluster), by = c("donor","cluster")) %>%
    mutate(sample_id = paste0(donor, "__cl", cluster))

  counts <- GetAssayData(obj_sub, assay = DefaultAssay(obj_sub), slot = "counts")
  pb_counts <- sparse_pseudobulk(counts[, pb_cells$cell_id, drop = FALSE], pb_cells$sample_id)
  pb_counts <- pb_counts[, pb_meta$sample_id, drop = FALSE]

  dge <- DGEList(counts = pb_counts)
  keep_gene <- filterByExpr(dge, group = factor(pb_meta$cluster_simple))
  dge <- dge[keep_gene, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)

  pb_meta$donor_factor <- factor(pb_meta$donor)
  pb_meta$cluster_simple <- factor(pb_meta$cluster_simple, levels = c("ref","cand"))

  design <- model.matrix(~ donor_factor + cluster_simple, data = pb_meta)
  v <- voom(dge, design = design, plot = FALSE)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)

  coef_name <- grep("^cluster_simplecand$", colnames(design), value = TRUE)
  if (length(coef_name) != 1) stop("Could not identify pseudobulk cluster coefficient.")
  de <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
  de$gene <- rownames(de)
  de <- de %>% arrange(adj.P.Val, desc(abs(logFC)))

  fwrite(pb_meta, file.path(outdir, paste0(out_prefix, "__pb_meta.tsv")), sep = "\t")
  fwrite(as.data.table(de), file.path(outdir, paste0(out_prefix, "__pb_DE.tsv.gz")), sep = "\t")
  saveRDS(pb_counts, file.path(outdir, paste0(out_prefix, "__pb_counts.rds")))

  pca <- prcomp(t(cpm(dge, log = TRUE, prior.count = 2)), center = TRUE, scale. = TRUE)
  plot_df <- data.frame(
    sample_id = rownames(pca$x),
    PC1 = pca$x[,1],
    PC2 = pca$x[,2],
    cluster = pb_meta$cluster,
    donor = pb_meta$donor,
    dx = pb_meta$dx,
    stringsAsFactors = FALSE
  )
  p <- ggplot(plot_df, aes(x = PC1, y = PC2, color = cluster, shape = dx)) +
    geom_point(size = 3) +
    theme_classic(base_size = 12) +
    labs(title = paste0(out_prefix, ": pseudobulk PCA"))
  ggsave(file.path(outdir, paste0(out_prefix, "__pb_PCA.pdf")), p, width = 7, height = 5.5)
  ggsave(file.path(outdir, paste0(out_prefix, "__pb_PCA.png")), p, width = 7, height = 5.5, dpi = 300)

  de
}

compare_de_tables <- function(base_de, filt_de, label_base, label_filt) {
  if (is.null(base_de) || is.null(filt_de)) return(NULL)
  x <- base_de[, c("gene","logFC","adj.P.Val")]
  y <- filt_de[, c("gene","logFC","adj.P.Val")]
  colnames(x) <- c("gene", paste0("logFC_", label_base), paste0("FDR_", label_base))
  colnames(y) <- c("gene", paste0("logFC_", label_filt), paste0("FDR_", label_filt))
  z <- merge(x, y, by = "gene", all = TRUE)
  z$direction_same <- sign(z[[paste0("logFC_", label_base)]]) == sign(z[[paste0("logFC_", label_filt)]])
  z
}

# ============================================================
# 1. Read object
# ============================================================
msg("Reading strict object ...")
obj <- readRDS(cfg$strict_rds)
if (!inherits(obj, "Seurat")) stop("strict_rds is not a Seurat object.")
if (cfg$assay %in% Assays(obj)) DefaultAssay(obj) <- cfg$assay

if (!cfg$cluster_col %in% colnames(obj@meta.data)) stop("cluster_col not found in object metadata.")
if (!cfg$donor_col %in% colnames(obj@meta.data)) stop("donor_col not found in object metadata.")
if (!cfg$dx_col %in% colnames(obj@meta.data)) stop("dx_col not found in object metadata.")

obj$cluster_use <- norm_cluster_label(obj@meta.data[[cfg$cluster_col]])
obj$donor_use <- as.character(obj@meta.data[[cfg$donor_col]])
obj$dx_use <- normalize_dx(obj@meta.data[[cfg$dx_col]])

cells_keep <- colnames(obj)[obj$cluster_use %in% c(cfg$ref_cluster, cfg$cand_cluster)]
if (length(cells_keep) == 0) stop("No cells found for requested clusters in strict object.")
obj <- subset(obj, cells = cells_keep)
obj$cluster_use <- factor(obj$cluster_use, levels = c(cfg$ref_cluster, cfg$cand_cluster))
Idents(obj) <- obj$cluster_use

msg("Retained cells by cluster:")
print(table(obj$cluster_use, useNA = "ifany"))

audit_meta <- data.frame(
  item = c("n_cells","n_features","default_assay","cluster_col","donor_col","dx_col","ref_cluster","cand_cluster"),
  value = c(ncol(obj), nrow(obj), DefaultAssay(obj), cfg$cluster_col, cfg$donor_col, cfg$dx_col, cfg$ref_cluster, cfg$cand_cluster),
  stringsAsFactors = FALSE
)
fwrite(audit_meta, file.path(cfg$outdir, "tables", "00_package2b_metadata.tsv"), sep = "\t")
saveRDS(obj, file.path(cfg$outdir, "objects", "Package2b_input_object.rds"))

# ============================================================
# 2. Add module scores
# ============================================================
msg("Computing module scores ...")
expr_data <- GetAssayData(obj, assay = DefaultAssay(obj), slot = "data")

panel_list <- list(
  microglia_homeostatic = cfg$microglia_homeostatic,
  microglia_activation = cfg$microglia_activation,
  neuronal_synaptic = cfg$neuronal_synaptic,
  astrocyte = cfg$astrocyte,
  oligodendrocyte = cfg$oligodendrocyte,
  opc = cfg$opc,
  endothelial = cfg$endothelial,
  lymphocyte = cfg$lymphocyte
)

genes_used_df <- rbindlist(lapply(names(panel_list), function(nm) {
  gs <- intersect(panel_list[[nm]], rownames(expr_data))
  data.frame(panel = nm, gene = gs, stringsAsFactors = FALSE)
}), fill = TRUE)
fwrite(genes_used_df, file.path(cfg$outdir, "tables", "01_module_genes_used.tsv"), sep = "\t")

for (nm in names(panel_list)) {
  res <- compute_module_score(expr_data, panel_list[[nm]])
  obj[[paste0("ms_", nm)]] <- res$score
}

saveRDS(obj, file.path(cfg$outdir, "objects", "Package2b_object_with_module_scores.rds"))

# ============================================================
# 3. Cell-level and donor-level module summaries
# ============================================================
msg("Summarizing module scores ...")
qc_cols_present <- cfg$qc_cols_try[cfg$qc_cols_try %in% colnames(obj@meta.data)]

cell_df <- obj@meta.data %>%
  mutate(
    cell_id = colnames(obj),
    cluster = as.character(cluster_use),
    donor = donor_use,
    dx = dx_use
  ) %>%
  select(cell_id, cluster, donor, dx,
         starts_with("ms_"),
         all_of(qc_cols_present))

fwrite(as.data.table(cell_df), file.path(cfg$outdir, "tables", "02_cell_level_module_and_qc.tsv.gz"), sep = "\t")

donor_df <- cell_df %>%
  group_by(donor, dx, cluster) %>%
  summarise(across(c(starts_with("ms_"), all_of(qc_cols_present)), ~mean(.x, na.rm = TRUE)), .groups = "drop")

fwrite(as.data.table(donor_df), file.path(cfg$outdir, "tables", "03_donor_level_module_and_qc.tsv.gz"), sep = "\t")

module_cols <- grep("^ms_", colnames(cell_df), value = TRUE)

cell_module_stats <- lapply(module_cols, function(vv) {
  data.frame(
    variable = vv,
    level = "cell",
    mean_cluster0 = mean(cell_df[[vv]][cell_df$cluster == cfg$ref_cluster], na.rm = TRUE),
    mean_cluster2 = mean(cell_df[[vv]][cell_df$cluster == cfg$cand_cluster], na.rm = TRUE),
    effect_cluster2_minus_0 = mean(cell_df[[vv]][cell_df$cluster == cfg$cand_cluster], na.rm = TRUE) -
      mean(cell_df[[vv]][cell_df$cluster == cfg$ref_cluster], na.rm = TRUE),
    wilcox_p = calc_wilcox(cell_df, vv, "cluster"),
    stringsAsFactors = FALSE
  )
})
cell_module_stats <- rbindlist(cell_module_stats, fill = TRUE)

donor_module_stats <- lapply(module_cols, function(vv) {
  data.frame(
    variable = vv,
    level = "donor",
    mean_cluster0 = mean(donor_df[[vv]][donor_df$cluster == cfg$ref_cluster], na.rm = TRUE),
    mean_cluster2 = mean(donor_df[[vv]][donor_df$cluster == cfg$cand_cluster], na.rm = TRUE),
    effect_cluster2_minus_0 = mean(donor_df[[vv]][donor_df$cluster == cfg$cand_cluster], na.rm = TRUE) -
      mean(donor_df[[vv]][donor_df$cluster == cfg$ref_cluster], na.rm = TRUE),
    paired_wilcox_p = paired_wilcox_by_donor(donor_df, vv),
    stringsAsFactors = FALSE
  )
})
donor_module_stats <- rbindlist(donor_module_stats, fill = TRUE)

fwrite(cell_module_stats, file.path(cfg$outdir, "tables", "04_cell_level_module_stats_cluster2_vs_0.tsv"), sep = "\t")
fwrite(donor_module_stats, file.path(cfg$outdir, "tables", "05_donor_level_module_stats_cluster2_vs_0.tsv"), sep = "\t")

# plots: module violins and donor paired
plot_modules_main <- c("ms_microglia_homeostatic","ms_microglia_activation","ms_neuronal_synaptic","ms_astrocyte","ms_oligodendrocyte","ms_endothelial","ms_lymphocyte")
plot_modules_main <- plot_modules_main[plot_modules_main %in% module_cols]

if (length(plot_modules_main) > 0) {
  p1 <- cell_df %>%
    select(cluster, all_of(plot_modules_main)) %>%
    tidyr::pivot_longer(cols = all_of(plot_modules_main), names_to = "module", values_to = "score") %>%
    ggplot(aes(x = cluster, y = score, fill = cluster)) +
    geom_violin(scale = "width", trim = TRUE) +
    facet_wrap(~ module, scales = "free_y") +
    theme_classic(base_size = 12) +
    labs(title = "Cell-level module scores: Cluster 0 vs Cluster 2", x = NULL, y = "Module score")
  ggsave(file.path(cfg$outdir, "plots", "06_cell_level_module_violin.pdf"), p1, width = 12, height = 7)
  ggsave(file.path(cfg$outdir, "plots", "06_cell_level_module_violin.png"), p1, width = 12, height = 7, dpi = 300)

  p2 <- donor_df %>%
    select(donor, cluster, dx, all_of(plot_modules_main)) %>%
    tidyr::pivot_longer(cols = all_of(plot_modules_main), names_to = "module", values_to = "score") %>%
    ggplot(aes(x = cluster, y = score, group = donor, color = dx)) +
    geom_line(alpha = 0.45) +
    geom_point(size = 1.6) +
    facet_wrap(~ module, scales = "free_y") +
    theme_classic(base_size = 12) +
    labs(title = "Donor-level paired module scores", x = NULL, y = "Average module score")
  ggsave(file.path(cfg$outdir, "plots", "07_donor_level_module_paired.pdf"), p2, width = 12, height = 7)
  ggsave(file.path(cfg$outdir, "plots", "07_donor_level_module_paired.png"), p2, width = 12, height = 7, dpi = 300)
}

# ============================================================
# 4. Gene-level residual background audit
# ============================================================
msg("Running gene-level background audit ...")
gene_panel <- unique(c(cfg$neuronal_background_genes, cfg$audit_genes))
gene_panel_use <- intersect(gene_panel, rownames(expr_data))

gene_summary_list <- lapply(gene_panel_use, function(g) {
  vv <- as.numeric(expr_data[g, ])
  data.frame(
    gene = g,
    cluster = as.character(cell_df$cluster),
    donor = cell_df$donor,
    dx = cell_df$dx,
    expr = vv,
    detected = as.integer(vv > 0),
    stringsAsFactors = FALSE
  )
})
gene_long <- rbindlist(gene_summary_list, fill = TRUE)

gene_cluster_summary <- gene_long %>%
  group_by(gene, cluster) %>%
  summarise(
    avg_expr = mean(expr, na.rm = TRUE),
    pct_exp = mean(detected, na.rm = TRUE),
    .groups = "drop"
  )
fwrite(as.data.table(gene_cluster_summary), file.path(cfg$outdir, "tables", "08_gene_cluster_summary.tsv.gz"), sep = "\t")

gene_donor_summary <- gene_long %>%
  group_by(gene, donor, dx, cluster) %>%
  summarise(
    avg_expr = mean(expr, na.rm = TRUE),
    pct_exp = mean(detected, na.rm = TRUE),
    .groups = "drop"
  )
fwrite(as.data.table(gene_donor_summary), file.path(cfg$outdir, "tables", "09_gene_donor_summary.tsv.gz"), sep = "\t")

if (length(gene_panel_use) > 0) {
  p3 <- DotPlot(obj, features = gene_panel_use, group.by = "cluster_use") +
    RotatedAxis() +
    theme_classic(base_size = 11) +
    labs(title = "Residual background audit genes across retained clusters")
  ggsave(file.path(cfg$outdir, "plots", "08_residual_background_gene_dotplot.pdf"), p3, width = 14, height = 5.5)
  ggsave(file.path(cfg$outdir, "plots", "08_residual_background_gene_dotplot.png"), p3, width = 14, height = 5.5, dpi = 300)
}

# ============================================================
# 5. QC / doublet correlation audit
# ============================================================
msg("Auditing QC and scrublet associations ...")
qc_numeric <- qc_cols_present

if (length(qc_numeric) > 0) {
  qc_stats <- lapply(qc_numeric, function(vv) {
    data.frame(
      variable = vv,
      mean_cluster0 = mean(cell_df[[vv]][cell_df$cluster == cfg$ref_cluster], na.rm = TRUE),
      mean_cluster2 = mean(cell_df[[vv]][cell_df$cluster == cfg$cand_cluster], na.rm = TRUE),
      effect_cluster2_minus_0 = mean(cell_df[[vv]][cell_df$cluster == cfg$cand_cluster], na.rm = TRUE) -
        mean(cell_df[[vv]][cell_df$cluster == cfg$ref_cluster], na.rm = TRUE),
      wilcox_p = calc_wilcox(cell_df, vv, "cluster"),
      stringsAsFactors = FALSE
    )
  })
  qc_stats <- rbindlist(qc_stats, fill = TRUE)
  fwrite(qc_stats, file.path(cfg$outdir, "tables", "10_qc_stats_cluster2_vs_0.tsv"), sep = "\t")

  p4 <- cell_df %>%
    select(cluster, all_of(qc_numeric)) %>%
    tidyr::pivot_longer(cols = all_of(qc_numeric), names_to = "qc", values_to = "value") %>%
    ggplot(aes(x = cluster, y = value, fill = cluster)) +
    geom_violin(scale = "width", trim = TRUE) +
    facet_wrap(~ qc, scales = "free_y") +
    theme_classic(base_size = 12) +
    labs(title = "QC / scrublet metrics by cluster", x = NULL, y = "Value")
  ggsave(file.path(cfg$outdir, "plots", "09_qc_violin_by_cluster.pdf"), p4, width = 11, height = 6.5)
  ggsave(file.path(cfg$outdir, "plots", "09_qc_violin_by_cluster.png"), p4, width = 11, height = 6.5, dpi = 300)
}

# scatter if columns present
if ("single_scrublet_score" %in% qc_numeric && "ms_neuronal_synaptic" %in% module_cols) {
  p5 <- ggplot(cell_df, aes(x = single_scrublet_score, y = ms_neuronal_synaptic, color = cluster)) +
    geom_point(alpha = 0.45, size = 0.7) +
    theme_classic(base_size = 12) +
    labs(title = "Neuronal module vs scrublet score")
  ggsave(file.path(cfg$outdir, "plots", "10_neuronal_module_vs_scrublet.pdf"), p5, width = 6.5, height = 5.5)
  ggsave(file.path(cfg$outdir, "plots", "10_neuronal_module_vs_scrublet.png"), p5, width = 6.5, height = 5.5, dpi = 300)
}

if ("nFeature_RNA" %in% qc_numeric && "ms_neuronal_synaptic" %in% module_cols) {
  p6 <- ggplot(cell_df, aes(x = nFeature_RNA, y = ms_neuronal_synaptic, color = cluster)) +
    geom_point(alpha = 0.45, size = 0.7) +
    theme_classic(base_size = 12) +
    labs(title = "Neuronal module vs nFeature_RNA")
  ggsave(file.path(cfg$outdir, "plots", "11_neuronal_module_vs_nFeature.pdf"), p6, width = 6.5, height = 5.5)
  ggsave(file.path(cfg$outdir, "plots", "11_neuronal_module_vs_nFeature.png"), p6, width = 6.5, height = 5.5, dpi = 300)
}

# ============================================================
# 6. Sensitivity filters
# ============================================================
msg("Running sensitivity filters ...")
if (!"ms_neuronal_synaptic" %in% colnames(obj@meta.data)) stop("Neuronal module score missing.")
neuronal_thr <- quantile(obj$ms_neuronal_synaptic, probs = cfg$neuronal_high_quantile, na.rm = TRUE)
if ("single_scrublet_score" %in% colnames(obj@meta.data)) {
  scrublet_thr <- quantile(obj$single_scrublet_score, probs = cfg$scrublet_high_quantile, na.rm = TRUE)
} else {
  scrublet_thr <- NA_real_
}

thr_df <- data.frame(
  neuronal_high_quantile = cfg$neuronal_high_quantile,
  neuronal_threshold = neuronal_thr,
  scrublet_high_quantile = cfg$scrublet_high_quantile,
  scrublet_threshold = scrublet_thr,
  stringsAsFactors = FALSE
)
fwrite(thr_df, file.path(cfg$outdir, "tables", "11_sensitivity_thresholds.tsv"), sep = "\t")

obj$keep_neuronal_filtered <- obj$ms_neuronal_synaptic <= neuronal_thr
if (!is.na(scrublet_thr)) {
  obj$keep_scrublet_filtered <- obj$single_scrublet_score <= scrublet_thr
} else {
  obj$keep_scrublet_filtered <- TRUE
}
obj$keep_combined_filtered <- obj$keep_neuronal_filtered & obj$keep_scrublet_filtered

saveRDS(obj, file.path(cfg$outdir, "objects", "Package2b_object_with_filter_flags.rds"))

filter_summary <- data.frame(
  filter = c("baseline","neuronal_filtered","scrublet_filtered","combined_filtered"),
  n_cells = c(
    ncol(obj),
    sum(obj$keep_neuronal_filtered, na.rm = TRUE),
    sum(obj$keep_scrublet_filtered, na.rm = TRUE),
    sum(obj$keep_combined_filtered, na.rm = TRUE)
  ),
  stringsAsFactors = FALSE
)
fwrite(filter_summary, file.path(cfg$outdir, "tables", "12_filter_summary.tsv"), sep = "\t")

obj_neur <- subset(obj, cells = colnames(obj)[obj$keep_neuronal_filtered %in% TRUE])
obj_scru <- subset(obj, cells = colnames(obj)[obj$keep_scrublet_filtered %in% TRUE])
obj_comb <- subset(obj, cells = colnames(obj)[obj$keep_combined_filtered %in% TRUE])

saveRDS(obj_neur, file.path(cfg$outdir, "objects", "Package2b_neuronal_filtered.rds"))
saveRDS(obj_scru, file.path(cfg$outdir, "objects", "Package2b_scrublet_filtered.rds"))
saveRDS(obj_comb, file.path(cfg$outdir, "objects", "Package2b_combined_filtered.rds"))

# ============================================================
# 7. Baseline and filtered pseudobulk DE
# ============================================================
msg("Running baseline and filtered paired pseudobulk DE ...")
base_de <- run_paired_pseudobulk(
  obj_sub = obj,
  donor_col = "donor_use",
  dx_col = "dx_use",
  cluster_col = "cluster_use",
  out_prefix = "13_baseline",
  outdir = file.path(cfg$outdir, "tables"),
  min_cells = cfg$min_cells_per_pb_sample
)

neur_de <- run_paired_pseudobulk(
  obj_sub = obj_neur,
  donor_col = "donor_use",
  dx_col = "dx_use",
  cluster_col = "cluster_use",
  out_prefix = "14_neuronal_filtered",
  outdir = file.path(cfg$outdir, "tables"),
  min_cells = cfg$min_cells_per_pb_sample
)

scru_de <- run_paired_pseudobulk(
  obj_sub = obj_scru,
  donor_col = "donor_use",
  dx_col = "dx_use",
  cluster_col = "cluster_use",
  out_prefix = "15_scrublet_filtered",
  outdir = file.path(cfg$outdir, "tables"),
  min_cells = cfg$min_cells_per_pb_sample
)

comb_de <- run_paired_pseudobulk(
  obj_sub = obj_comb,
  donor_col = "donor_use",
  dx_col = "dx_use",
  cluster_col = "cluster_use",
  out_prefix = "16_combined_filtered",
  outdir = file.path(cfg$outdir, "tables"),
  min_cells = cfg$min_cells_per_pb_sample
)

cmp_neur <- compare_de_tables(base_de, neur_de, "baseline", "neuronal_filtered")
cmp_scru <- compare_de_tables(base_de, scru_de, "baseline", "scrublet_filtered")
cmp_comb <- compare_de_tables(base_de, comb_de, "baseline", "combined_filtered")

if (!is.null(cmp_neur)) fwrite(cmp_neur, file.path(cfg$outdir, "tables", "17_compare_baseline_vs_neuronal_filtered.tsv.gz"), sep = "\t")
if (!is.null(cmp_scru)) fwrite(cmp_scru, file.path(cfg$outdir, "tables", "18_compare_baseline_vs_scrublet_filtered.tsv.gz"), sep = "\t")
if (!is.null(cmp_comb)) fwrite(cmp_comb, file.path(cfg$outdir, "tables", "19_compare_baseline_vs_combined_filtered.tsv.gz"), sep = "\t")

# compact top-gene tables
extract_top <- function(de, prefix) {
  if (is.null(de)) return(NULL)
  up2 <- de %>% filter(adj.P.Val < 0.05, logFC > 0) %>% head(30)
  up0 <- de %>% filter(adj.P.Val < 0.05, logFC < 0) %>% head(30)
  fwrite(up2, file.path(cfg$outdir, "tables", paste0(prefix, "__top_up_in_cluster2.tsv")), sep = "\t")
  fwrite(up0, file.path(cfg$outdir, "tables", paste0(prefix, "__top_up_in_cluster0.tsv")), sep = "\t")
}
extract_top(base_de, "20_baseline")
extract_top(neur_de, "21_neuronal_filtered")
extract_top(scru_de, "22_scrublet_filtered")
extract_top(comb_de, "23_combined_filtered")

# ============================================================
# 8. Manuscript-style compact summary
# ============================================================
msg("Writing compact summary ...")
summary_main <- donor_module_stats %>%
  mutate(panel = sub("^ms_", "", variable)) %>%
  arrange(desc(effect_cluster2_minus_0))

fwrite(summary_main, file.path(cfg$outdir, "tables", "24_module_summary_for_manuscript.tsv"), sep = "\t")

resid_gene_summary <- gene_cluster_summary %>%
  mutate(cluster = as.character(cluster)) %>%
  tidyr::pivot_wider(
    names_from = cluster,
    values_from = c(avg_expr, pct_exp),
    names_glue = "{.value}_cluster_{cluster}"
  )

# 如果某个 cluster 列缺失，先补 NA，避免 mutate 再报错
for (nm in c("avg_expr_cluster_0", "avg_expr_cluster_2", "pct_exp_cluster_0", "pct_exp_cluster_2")) {
  if (!nm %in% colnames(resid_gene_summary)) {
    resid_gene_summary[[nm]] <- NA_real_
  }
}

resid_gene_summary <- resid_gene_summary %>%
  mutate(
    avg_expr_diff_2_minus_0 = avg_expr_cluster_2 - avg_expr_cluster_0,
    pct_exp_diff_2_minus_0 = pct_exp_cluster_2 - pct_exp_cluster_0
  ) %>%
  arrange(desc(abs(avg_expr_diff_2_minus_0)))

fwrite(
  resid_gene_summary,
  file.path(cfg$outdir, "tables", "25_residual_gene_summary_for_manuscript.tsv"),
  sep = "\t"
)
manifest <- data.frame(
  output_file = c(
    "04_cell_level_module_stats_cluster2_vs_0.tsv",
    "05_donor_level_module_stats_cluster2_vs_0.tsv",
    "08_gene_cluster_summary.tsv.gz",
    "09_gene_donor_summary.tsv.gz",
    "10_qc_stats_cluster2_vs_0.tsv",
    "13_baseline__pb_DE.tsv.gz",
    "17_compare_baseline_vs_neuronal_filtered.tsv.gz",
    "19_compare_baseline_vs_combined_filtered.tsv.gz",
    "24_module_summary_for_manuscript.tsv",
    "25_residual_gene_summary_for_manuscript.tsv"
  ),
  description = c(
    "cell-level module score comparison between clusters 2 and 0",
    "donor-level paired module score comparison between clusters 2 and 0",
    "gene-level average expression and detection summary by cluster",
    "gene-level donor summary by cluster",
    "QC and scrublet comparison by cluster",
    "baseline paired pseudobulk DE",
    "baseline vs neuronal-filtered DE comparison",
    "baseline vs combined-filtered DE comparison",
    "compact module summary for manuscript/rebuttal",
    "compact residual-background gene summary for manuscript/rebuttal"
  ),
  stringsAsFactors = FALSE
)
fwrite(manifest, file.path(cfg$outdir, "PACKAGE2b_manifest.tsv"), sep = "\t")

save_session_info(file.path(cfg$outdir, "sessionInfo_package2b.txt"))
msg("Package 2b finished.")
