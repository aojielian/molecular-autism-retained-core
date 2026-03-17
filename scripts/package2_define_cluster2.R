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
})

options(stringsAsFactors = FALSE)

cfg <- list(
  strict_rds = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Package2_strict_input_fixed_from_original.rds",
  global_rds = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/PsychENCODE_Science_2024/PsychENCODE_global_object.rds",   # 若要做 global background audit，请填全局对象路径；否则保持 NA
  outdir = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Package2_Define_Cluster2",
  assay = "RNA",

  strict_cluster_col = NULL,
  strict_donor_col = NULL,
  strict_dx_col = NULL,

  global_celltype_col = "annotation",  # 若 global_rds 提供，尽量自动识别；必要时手动填

  ref_cluster = "0",
  cand_cluster = "2",

  min_cells_per_pb_sample = 20,
  top_n_core = 50,
  top_n_global_audit = 15,

  audit_genes = c("SPP1","GAS6","APOE","MERTK","CD74","C3"),
  homeostatic_genes = c("P2RY12","CX3CR1","TMEM119","CSF1R","AIF1","C1QA","C1QB","C1QC"),
  optional_gmt_files = character(0)  # 可填多个 GMT 文件路径；若为空则跳过 fgsea
)

dir.create(cfg$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cfg$outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cfg$outdir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cfg$outdir, "objects"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cfg$outdir, "fgsea"), recursive = TRUE, showWarnings = FALSE)

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
  cand <- find_first_col(meta, c("^cluster_use$", "^cluster$", "cluster", "seurat", "res\\.", "RNA_snn", "wsnn"))
  if (!is.null(cand)) return(cand)
  obj$cluster_auto_from_idents <- as.character(Idents(obj))
  return("cluster_auto_from_idents")
}

choose_donor_col <- function(meta, user_col = NULL) {
  if (!is.null(user_col) && user_col %in% names(meta)) return(user_col)
  find_first_col(meta, c("^donor$", "donor_id", "individual", "subject", "patient", "sample", "orig.ident", "library", "case_id"))
}

choose_dx_col <- function(meta, user_col = NULL) {
  if (!is.null(user_col) && user_col %in% names(meta)) return(user_col)
  find_first_col(meta, c("^dx$", "^group$", "diagnosis", "condition", "status", "phenotype", "disease"))
}

choose_celltype_col <- function(meta, user_col = NULL) {
  if (!is.null(user_col) && user_col %in% names(meta)) return(user_col)
  find_first_col(meta, c("^celltype$", "cell_type", "annotation", "annot", "class", "major", "broad"))
}

get_fc_col <- function(df) {
  fc_candidates <- c("avg_log2FC", "avg_logFC", "log2FC", "logFC")
  hit <- fc_candidates[fc_candidates %in% colnames(df)]
  if (length(hit) == 0) return(NULL)
  hit[1]
}

safe_scale <- function(x) {
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  if (sd(x, na.rm = TRUE) == 0) return(rep(0, length(x)))
  as.numeric(scale(x))
}

parse_gmt <- function(gmt_file) {
  lines <- readLines(gmt_file)
  res <- lapply(lines, function(x) {
    sp <- strsplit(x, "\t")[[1]]
    if (length(sp) < 3) return(NULL)
    name <- sp[1]
    genes <- unique(sp[-c(1,2)])
    list(name = name, genes = genes)
  })
  res <- res[!sapply(res, is.null)]
  gs <- lapply(res, `[[`, "genes")
  names(gs) <- sapply(res, `[[`, "name")
  gs
}

sparse_pseudobulk <- function(count_mat, groups) {
  groups <- factor(groups)
  mm <- sparse.model.matrix(~ 0 + groups)
  colnames(mm) <- sub("^groups", "", colnames(mm))
  pb <- count_mat %*% mm
  pb
}

make_pb_pca_plot <- function(logcpm_mat, meta_df, out_pdf, out_png, title_txt) {
  pca <- prcomp(t(logcpm_mat), center = TRUE, scale. = TRUE)
  plot_df <- data.frame(
    sample_id = rownames(pca$x),
    PC1 = pca$x[,1],
    PC2 = pca$x[,2],
    cluster = meta_df$cluster,
    donor = meta_df$donor,
    dx = meta_df$dx,
    stringsAsFactors = FALSE
  )
  p <- ggplot(plot_df, aes(x = PC1, y = PC2, color = cluster, shape = dx)) +
    geom_point(size = 3) +
    theme_classic(base_size = 12) +
    labs(title = title_txt)
  ggsave(out_pdf, p, width = 7, height = 5.5)
  ggsave(out_png, p, width = 7, height = 5.5, dpi = 300)
}

run_fgsea_if_possible <- function(rank_vec, gmt_files, outdir) {
  if (length(gmt_files) == 0) return(invisible(NULL))
  if (!requireNamespace("fgsea", quietly = TRUE)) {
    writeLines("fgsea not installed; pathway analysis skipped.", file.path(outdir, "FGSEA_SKIPPED.txt"))
    return(invisible(NULL))
  }

  for (gmt in gmt_files) {
    if (!file.exists(gmt)) next
    gs <- parse_gmt(gmt)
    gs <- gs[sapply(gs, length) >= 10 & sapply(gs, length) <= 500]
    if (length(gs) == 0) next
    fg <- fgsea::fgseaMultilevel(pathways = gs, stats = rank_vec)
    fg <- as.data.frame(fg) %>% arrange(padj, desc(abs(NES)))
    base <- tools::file_path_sans_ext(basename(gmt))
    fwrite(fg, file.path(outdir, paste0(base, "__fgsea.tsv.gz")), sep = "\t")

    top_show <- fg %>% filter(!is.na(padj)) %>% arrange(padj) %>% head(20)
    if (nrow(top_show) > 0) {
      top_show$pathway <- factor(top_show$pathway, levels = rev(top_show$pathway))
      p <- ggplot(top_show, aes(x = NES, y = pathway, size = -log10(padj))) +
        geom_point() +
        theme_classic(base_size = 11) +
        labs(title = paste0("FGSEA: ", base), y = NULL)
      ggsave(file.path(outdir, paste0(base, "__fgsea_top20.pdf")), p, width = 7.5, height = 6)
      ggsave(file.path(outdir, paste0(base, "__fgsea_top20.png")), p, width = 7.5, height = 6, dpi = 300)
    }
  }
}

# ============================================================
# 1. Read strict object
# ============================================================
msg("Reading strict object ...")
obj <- readRDS(cfg$strict_rds)
if (!inherits(obj, "Seurat")) stop("strict_rds is not a Seurat object.")

if (cfg$assay %in% Assays(obj)) DefaultAssay(obj) <- cfg$assay

cluster_col <- choose_cluster_col(obj, cfg$strict_cluster_col)
donor_col <- choose_donor_col(obj@meta.data, cfg$strict_donor_col)
dx_col <- choose_dx_col(obj@meta.data, cfg$strict_dx_col)

if (is.null(cluster_col)) stop("Could not identify cluster column in strict object.")
if (is.null(donor_col)) stop("Could not identify donor column in strict object.")
if (is.null(dx_col)) stop("Could not identify diagnosis column in strict object.")

obj$cluster_use <- as.character(obj@meta.data[[cluster_col]])
obj$donor_use <- as.character(obj@meta.data[[donor_col]])
obj$dx_use <- normalize_dx(obj@meta.data[[dx_col]])
Idents(obj) <- obj$cluster_use

keep_clusters <- c(cfg$ref_cluster, cfg$cand_cluster)
obj <- subset(obj, cells = colnames(obj)[obj$cluster_use %in% keep_clusters])
obj$cluster_use <- factor(obj$cluster_use, levels = keep_clusters)
Idents(obj) <- obj$cluster_use

audit_meta <- data.frame(
  item = c("n_cells","n_features","default_assay","cluster_col","donor_col","dx_col","ref_cluster","cand_cluster"),
  value = c(ncol(obj), nrow(obj), DefaultAssay(obj), cluster_col, donor_col, dx_col, cfg$ref_cluster, cfg$cand_cluster),
  stringsAsFactors = FALSE
)
fwrite(audit_meta, file.path(cfg$outdir, "tables", "00_package2_metadata.tsv"), sep = "\t")

saveRDS(obj, file.path(cfg$outdir, "objects", "Package2_strict_input_0_vs_2.rds"))

# ============================================================
# 2. Basic summary and plots
# ============================================================
msg("Saving basic summaries ...")
cell_count_df <- obj@meta.data %>%
  mutate(cluster = as.character(cluster_use), donor = donor_use, dx = dx_use) %>%
  count(cluster, dx, donor, name = "n_cells") %>%
  arrange(cluster, dx, donor)
fwrite(cell_count_df, file.path(cfg$outdir, "tables", "01_cell_counts_by_cluster_dx_donor.tsv"), sep = "\t")

cluster_summary <- obj@meta.data %>%
  mutate(cluster = as.character(cluster_use), dx = dx_use) %>%
  count(cluster, dx, name = "n_cells") %>%
  arrange(cluster, dx)
fwrite(cluster_summary, file.path(cfg$outdir, "tables", "02_cluster_summary.tsv"), sep = "\t")

if ("umap" %in% names(obj@reductions)) {
  p_umap <- DimPlot(obj, reduction = "umap", group.by = "cluster_use", label = TRUE, repel = TRUE, raster = FALSE) +
    theme_classic(base_size = 12) +
    labs(title = "Strict retained microglia object: Cluster 0 vs Cluster 2")
  ggsave(file.path(cfg$outdir, "plots", "03_umap_cluster0_vs_2.pdf"), p_umap, width = 7, height = 6)
  ggsave(file.path(cfg$outdir, "plots", "03_umap_cluster0_vs_2.png"), p_umap, width = 7, height = 6, dpi = 300)
}

audit_panel <- unique(c(cfg$audit_genes, cfg$homeostatic_genes))
audit_panel_use <- intersect(audit_panel, rownames(obj))

if (length(audit_panel_use) > 0) {
  p_dot <- DotPlot(obj, features = audit_panel_use, group.by = "cluster_use") +
    RotatedAxis() +
    theme_classic(base_size = 12) +
    labs(title = "Audit panel in strict microglia object")
  ggsave(file.path(cfg$outdir, "plots", "04_audit_panel_dotplot_strict.pdf"), p_dot, width = 9, height = 4.5)
  ggsave(file.path(cfg$outdir, "plots", "04_audit_panel_dotplot_strict.png"), p_dot, width = 9, height = 4.5, dpi = 300)

  if ("umap" %in% names(obj@reductions)) {
    fp_list <- lapply(audit_panel_use, function(g) {
      FeaturePlot(obj, features = g, reduction = "umap", raster = FALSE) +
        theme_classic(base_size = 10) +
        ggtitle(g)
    })
    p_fp <- wrap_plots(fp_list, ncol = 3)
    ggsave(file.path(cfg$outdir, "plots", "05_audit_panel_featureplots_strict.pdf"), p_fp, width = 10, height = ceiling(length(fp_list)/3) * 3.2)
    ggsave(file.path(cfg$outdir, "plots", "05_audit_panel_featureplots_strict.png"), p_fp, width = 10, height = ceiling(length(fp_list)/3) * 3.2, dpi = 300)
  }
}

# ============================================================
# 3. Cell-level DE: Cluster 2 vs Cluster 0
# ============================================================
msg("Running cell-level DE ...")
cell_de <- FindMarkers(
  object = obj,
  ident.1 = cfg$cand_cluster,
  ident.2 = cfg$ref_cluster,
  logfc.threshold = 0,
  min.pct = 0.05,
  test.use = "wilcox"
)

cell_de$gene <- rownames(cell_de)
cell_fc_col <- get_fc_col(cell_de)
if (is.null(cell_fc_col)) stop("Could not identify cell-level FC column.")

cell_de <- cell_de %>%
  arrange(p_val_adj, desc(abs(.data[[cell_fc_col]])))
fwrite(as.data.table(cell_de), file.path(cfg$outdir, "tables", "06_celllevel_DE_cluster2_vs_0.tsv.gz"), sep = "\t")

# ============================================================
# 4. Donor-paired pseudobulk DE
# ============================================================
msg("Building donor-paired pseudobulk matrix ...")
md <- obj@meta.data
md$cell_id <- colnames(obj)
md$cluster <- as.character(md$cluster_use)
md$donor <- as.character(md$donor_use)
md$dx <- as.character(md$dx_use)

count_per_sample <- md %>%
  count(donor, cluster, dx, name = "n_cells") %>%
  filter(cluster %in% keep_clusters)

# keep donor-cluster samples with enough cells
count_per_sample <- count_per_sample %>%
  filter(n_cells >= cfg$min_cells_per_pb_sample)

# paired donors: must have both clusters after thresholding
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

fwrite(pb_meta, file.path(cfg$outdir, "tables", "07_pseudobulk_sample_metadata.tsv"), sep = "\t")

pb_cells <- md %>%
  filter(donor %in% paired_donors, cluster %in% keep_clusters) %>%
  semi_join(pb_meta %>% select(donor, cluster), by = c("donor","cluster")) %>%
  mutate(sample_id = paste0(donor, "__cl", cluster))

counts <- GetAssayData(obj, assay = DefaultAssay(obj), slot = "counts")
pb_counts <- sparse_pseudobulk(counts[, pb_cells$cell_id, drop = FALSE], pb_cells$sample_id)

# order columns to pb_meta sample order
pb_counts <- pb_counts[, pb_meta$sample_id, drop = FALSE]

fwrite(as.data.table(as.matrix(pb_counts[1:min(10, nrow(pb_counts)), , drop = FALSE]), keep.rownames = "gene"),
       file.path(cfg$outdir, "tables", "08_pseudobulk_counts_preview_top10genes.tsv"),
       sep = "\t")

saveRDS(pb_counts, file.path(cfg$outdir, "objects", "Package2_pseudobulk_counts.rds"))

msg("Running edgeR/limma paired pseudobulk DE ...")
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
pb_de <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P")
pb_de$gene <- rownames(pb_de)
pb_de <- pb_de %>%
  arrange(adj.P.Val, desc(abs(logFC)))
fwrite(as.data.table(pb_de), file.path(cfg$outdir, "tables", "09_pseudobulk_DE_cluster2_vs_0.tsv.gz"), sep = "\t")

# paired pseudobulk PCA
logcpm <- cpm(dge, log = TRUE, prior.count = 2)
make_pb_pca_plot(
  logcpm_mat = logcpm,
  meta_df = pb_meta,
  out_pdf = file.path(cfg$outdir, "plots", "10_pseudobulk_PCA.pdf"),
  out_png = file.path(cfg$outdir, "plots", "10_pseudobulk_PCA.png"),
  title_txt = "Pseudobulk PCA (paired donors; Cluster 2 vs Cluster 0)"
)

# volcano
pb_vol <- pb_de %>%
  mutate(
    neglog10FDR = -log10(pmax(adj.P.Val, 1e-300)),
    sig = case_when(
      adj.P.Val < 0.05 & logFC > 0 ~ "up_in_cluster2",
      adj.P.Val < 0.05 & logFC < 0 ~ "up_in_cluster0",
      TRUE ~ "ns"
    )
  )

p_vol <- ggplot(pb_vol, aes(x = logFC, y = neglog10FDR, color = sig)) +
  geom_point(size = 1.0, alpha = 0.7) +
  theme_classic(base_size = 12) +
  labs(title = "Pseudobulk DE: Cluster 2 vs Cluster 0", y = "-log10(FDR)")
ggsave(file.path(cfg$outdir, "plots", "11_pseudobulk_volcano.pdf"), p_vol, width = 7, height = 5.5)
ggsave(file.path(cfg$outdir, "plots", "11_pseudobulk_volcano.png"), p_vol, width = 7, height = 5.5, dpi = 300)

# ============================================================
# 5. Cell-level vs pseudobulk concordance
# ============================================================
msg("Building concordance table ...")
cell_tab <- cell_de[, c("gene", cell_fc_col, "p_val_adj", "pct.1", "pct.2")]
colnames(cell_tab) <- c("gene", "cell_logFC", "cell_FDR", "pct_cluster2", "pct_cluster0")

pb_tab <- pb_de[, c("gene", "logFC", "adj.P.Val", "AveExpr", "t", "P.Value")]
colnames(pb_tab) <- c("gene", "pb_logFC", "pb_FDR", "pb_AveExpr", "pb_t", "pb_P")

concord <- merge(cell_tab, pb_tab, by = "gene", all = TRUE)
concord$direction_concordant <- with(concord, sign(cell_logFC) == sign(pb_logFC))
concord$cell_sig <- !is.na(concord$cell_FDR) & concord$cell_FDR < 0.05
concord$pb_sig <- !is.na(concord$pb_FDR) & concord$pb_FDR < 0.05
concord$both_sig <- concord$cell_sig & concord$pb_sig
concord$both_sig_same_direction <- concord$both_sig & concord$direction_concordant

# data-driven core genes
core_up <- concord %>%
  filter(
    both_sig_same_direction,
    pb_logFC > 0,
    cell_logFC > 0
  ) %>%
  arrange(pb_FDR, desc(abs(pb_logFC)), cell_FDR) %>%
  head(cfg$top_n_core)

core_down <- concord %>%
  filter(
    both_sig_same_direction,
    pb_logFC < 0,
    cell_logFC < 0
  ) %>%
  arrange(pb_FDR, desc(abs(pb_logFC)), cell_FDR) %>%
  head(cfg$top_n_core)

fwrite(as.data.table(concord), file.path(cfg$outdir, "tables", "12_cell_vs_pseudobulk_concordance.tsv.gz"), sep = "\t")
fwrite(as.data.table(core_up), file.path(cfg$outdir, "tables", "13_data_driven_core_genes_up_in_cluster2.tsv"), sep = "\t")
fwrite(as.data.table(core_down), file.path(cfg$outdir, "tables", "14_data_driven_core_genes_up_in_cluster0.tsv"), sep = "\t")

# top concordant plot
top_concord <- concord %>%
  filter(!is.na(pb_logFC), !is.na(cell_logFC)) %>%
  mutate(cat = ifelse(both_sig_same_direction, "both_sig_same_direction", "other"))

p_concord <- ggplot(top_concord, aes(x = cell_logFC, y = pb_logFC, color = cat)) +
  geom_point(alpha = 0.7, size = 1.1) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
  geom_hline(yintercept = 0, linetype = 2, color = "grey60") +
  theme_classic(base_size = 12) +
  labs(title = "Cell-level vs pseudobulk concordance", x = "Cell-level logFC (2 vs 0)", y = "Pseudobulk logFC (2 vs 0)")
ggsave(file.path(cfg$outdir, "plots", "15_concordance_scatter.pdf"), p_concord, width = 6.5, height = 5.5)
ggsave(file.path(cfg$outdir, "plots", "15_concordance_scatter.png"), p_concord, width = 6.5, height = 5.5, dpi = 300)

# ============================================================
# 6. Per-donor average expression for audit genes
# ============================================================
msg("Computing per-donor average expression for audit genes ...")
audit_plus_core <- unique(c(
  cfg$audit_genes,
  cfg$homeostatic_genes,
  head(core_up$gene, cfg$top_n_global_audit),
  head(core_down$gene, cfg$top_n_global_audit)
))
audit_plus_core_use <- intersect(audit_plus_core, rownames(obj))

expr_data <- GetAssayData(obj, assay = DefaultAssay(obj), slot = "data")
audit_cell_meta <- obj@meta.data
audit_cell_meta$cell_id <- colnames(obj)
audit_cell_meta$cluster <- as.character(audit_cell_meta$cluster_use)
audit_cell_meta$donor <- as.character(audit_cell_meta$donor_use)
audit_cell_meta$dx <- as.character(audit_cell_meta$dx_use)

if (length(audit_plus_core_use) > 0) {
  audit_long_list <- lapply(audit_plus_core_use, function(g) {
    v <- as.numeric(expr_data[g, ])
    data.frame(
      gene = g,
      expr = v,
      donor = audit_cell_meta$donor,
      cluster = audit_cell_meta$cluster,
      dx = audit_cell_meta$dx,
      stringsAsFactors = FALSE
    )
  })
  audit_long <- rbindlist(audit_long_list, fill = TRUE)

  donor_gene_avg <- audit_long %>%
    group_by(gene, donor, cluster, dx) %>%
    summarise(avg_expr = mean(expr, na.rm = TRUE), .groups = "drop")

  fwrite(donor_gene_avg, file.path(cfg$outdir, "tables", "16_donor_gene_average_expression_strict.tsv.gz"), sep = "\t")

  plot_genes <- intersect(c(cfg$audit_genes, cfg$homeostatic_genes), unique(donor_gene_avg$gene))
  if (length(plot_genes) > 0) {
    p_gene <- donor_gene_avg %>%
      filter(gene %in% plot_genes) %>%
      mutate(cluster = factor(cluster, levels = keep_clusters)) %>%
      ggplot(aes(x = cluster, y = avg_expr, color = dx)) +
      geom_boxplot(outlier.shape = NA) +
      geom_jitter(width = 0.12, size = 1.1) +
      facet_wrap(~ gene, scales = "free_y") +
      theme_classic(base_size = 12) +
      labs(title = "Per-donor average expression in strict object", x = NULL, y = "Average expression")
    ggsave(file.path(cfg$outdir, "plots", "17_donor_average_expression_audit_panel.pdf"), p_gene, width = 11, height = 7)
    ggsave(file.path(cfg$outdir, "plots", "17_donor_average_expression_audit_panel.png"), p_gene, width = 11, height = 7, dpi = 300)
  }
}

# ============================================================
# 7. Optional FGSEA on pseudobulk ranked list
# ============================================================
msg("Running optional FGSEA if possible ...")
rank_vec <- pb_de$pb_t <- pb_de$t
names(rank_vec) <- pb_de$gene
rank_vec <- rank_vec[!is.na(rank_vec)]
rank_vec <- sort(rank_vec, decreasing = TRUE)

run_fgsea_if_possible(
  rank_vec = rank_vec,
  gmt_files = cfg$optional_gmt_files,
  outdir = file.path(cfg$outdir, "fgsea")
)

# ============================================================
# 8. Optional global background audit
# ============================================================
if (!is.na(cfg$global_rds) && file.exists(cfg$global_rds)) {
  msg("Reading global object for all-cell-type audit ...")
  gobj <- readRDS(cfg$global_rds)
  if (!inherits(gobj, "Seurat")) stop("global_rds is not a Seurat object.")
  if (cfg$assay %in% Assays(gobj)) DefaultAssay(gobj) <- cfg$assay

  global_celltype_col <- choose_celltype_col(gobj@meta.data, cfg$global_celltype_col)
  if (is.null(global_celltype_col)) {
    writeLines("Could not identify global cell type column; global audit skipped.",
               file.path(cfg$outdir, "tables", "18_GLOBAL_AUDIT_SKIPPED.txt"))
  } else {
    gobj$global_celltype_use <- as.character(gobj@meta.data[[global_celltype_col]])
    global_genes_use <- intersect(audit_plus_core, rownames(gobj))

    fwrite(
      data.frame(global_celltype_col = global_celltype_col, stringsAsFactors = FALSE),
      file.path(cfg$outdir, "tables", "18_global_celltype_column_used.tsv"),
      sep = "\t"
    )

    if (length(global_genes_use) > 0) {
      avg_global <- AverageExpression(
        gobj,
        assays = DefaultAssay(gobj),
        group.by = "global_celltype_use",
        slot = "data",
        verbose = FALSE
      )[[DefaultAssay(gobj)]]

      avg_global <- avg_global[intersect(rownames(avg_global), global_genes_use), , drop = FALSE]
      fwrite(as.data.table(as.matrix(avg_global), keep.rownames = "gene"),
             file.path(cfg$outdir, "tables", "19_global_average_expression_selected_genes.tsv.gz"),
             sep = "\t")

      # audit panel only
      global_audit_panel <- intersect(c(cfg$audit_genes, cfg$homeostatic_genes), global_genes_use)
      if (length(global_audit_panel) > 0) {
        p_gdot1 <- DotPlot(gobj, features = global_audit_panel, group.by = "global_celltype_use") +
          RotatedAxis() +
          theme_classic(base_size = 11) +
          labs(title = "Global background audit: audit panel across major cell types")
        ggsave(file.path(cfg$outdir, "plots", "20_global_dotplot_audit_panel.pdf"), p_gdot1, width = 11, height = 5.5)
        ggsave(file.path(cfg$outdir, "plots", "20_global_dotplot_audit_panel.png"), p_gdot1, width = 11, height = 5.5, dpi = 300)
      }

      # data-driven top genes
      global_core_panel <- unique(c(head(core_up$gene, cfg$top_n_global_audit), head(core_down$gene, cfg$top_n_global_audit)))
      global_core_panel <- intersect(global_core_panel, global_genes_use)
      if (length(global_core_panel) > 0) {
        p_gdot2 <- DotPlot(gobj, features = global_core_panel, group.by = "global_celltype_use") +
          RotatedAxis() +
          theme_classic(base_size = 10) +
          labs(title = "Global background audit: data-driven core genes")
        ggsave(file.path(cfg$outdir, "plots", "21_global_dotplot_core_genes.pdf"), p_gdot2, width = 14, height = 6)
        ggsave(file.path(cfg$outdir, "plots", "21_global_dotplot_core_genes.png"), p_gdot2, width = 14, height = 6, dpi = 300)
      }
    }
  }
} else {
  writeLines("global_rds not provided or not found; global background audit skipped.",
             file.path(cfg$outdir, "tables", "18_GLOBAL_AUDIT_SKIPPED.txt"))
}

# ============================================================
# 9. Compact manuscript-oriented summary
# ============================================================
msg("Writing compact summary tables ...")
summary_up <- core_up %>%
  select(gene, cell_logFC, cell_FDR, pb_logFC, pb_FDR, pct_cluster2, pct_cluster0)

summary_down <- core_down %>%
  select(gene, cell_logFC, cell_FDR, pb_logFC, pb_FDR, pct_cluster2, pct_cluster0)

fwrite(summary_up, file.path(cfg$outdir, "tables", "22_cluster2_core_summary_up.tsv"), sep = "\t")
fwrite(summary_down, file.path(cfg$outdir, "tables", "23_cluster0_core_summary_up.tsv"), sep = "\t")

manifest <- data.frame(
  output_file = c(
    "06_celllevel_DE_cluster2_vs_0.tsv.gz",
    "09_pseudobulk_DE_cluster2_vs_0.tsv.gz",
    "12_cell_vs_pseudobulk_concordance.tsv.gz",
    "13_data_driven_core_genes_up_in_cluster2.tsv",
    "14_data_driven_core_genes_up_in_cluster0.tsv",
    "17_donor_average_expression_audit_panel",
    "20_global_dotplot_audit_panel",
    "21_global_dotplot_core_genes"
  ),
  description = c(
    "cell-level differential expression, Cluster 2 vs Cluster 0",
    "donor-paired pseudobulk differential expression, Cluster 2 vs Cluster 0",
    "concordance between cell-level and pseudobulk results",
    "data-driven core genes up in Cluster 2",
    "data-driven core genes up in Cluster 0",
    "per-donor audit panel expression in strict object",
    "all-major-cell-type background audit for audit panel",
    "all-major-cell-type background audit for data-driven core genes"
  ),
  stringsAsFactors = FALSE
)
fwrite(manifest, file.path(cfg$outdir, "PACKAGE2_manifest.tsv"), sep = "\t")

save_session_info(file.path(cfg$outdir, "sessionInfo_package2.txt"))
msg("Package 2 finished.")
