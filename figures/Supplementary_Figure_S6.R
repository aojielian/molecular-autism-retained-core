suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(optparse)
})

option_list <- list(
  make_option("--atlas_rds", type = "character",
              default = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/PsychENCODE_Science_2024/PsychENCODE_global_object.rds"),
  make_option("--retained_rds", type = "character",
              default = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Package2_strict_input_fixed_from_original.rds"),
  make_option("--outdir", type = "character",
              default = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/FigureS6_SPP1_context"),
  make_option("--gene", type = "character", default = "SPP1"),
  make_option("--atlas_celltype_col", type = "character", default = "annotation"),
  make_option("--atlas_dx_col", type = "character", default = "Diagnosis"),
  make_option("--atlas_donor_col", type = "character", default = "individual_ID"),
  make_option("--ret_dx_col", type = "character", default = NULL),
  make_option("--ret_donor_col", type = "character", default = NULL),
  make_option("--ret_cluster_col", type = "character", default = NULL),
  make_option("--min_cells_per_donor_celltype", type = "integer", default = 10),
  make_option("--min_cells_per_donor_cluster", type = "integer", default = 5)
)
opt <- parse_args(OptionParser(option_list = option_list))

log_msg <- function(...) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] ", ts), sprintf(...), "\n", sep = "")
}

find_first_col <- function(x, patterns) {
  for (p in patterns) {
    hit <- grep(p, x, ignore.case = TRUE, value = TRUE)
    if (length(hit) > 0) return(hit[1])
  }
  NULL
}

choose_dx_col <- function(meta) {
  find_first_col(colnames(meta), c("^dx$", "^group$", "diagnosis", "condition", "status", "phenotype", "disease"))
}

choose_donor_col <- function(meta) {
  find_first_col(colnames(meta), c("^donor$", "donor_id", "individual", "subject", "patient", "sample", "orig.ident", "library", "case_id"))
}

choose_cluster_col <- function(meta) {
  find_first_col(colnames(meta), c("^cluster_use$", "^cluster$", "seurat_clusters", "res\\.", "RNA_snn", "wsnn", "louvain"))
}

normalize_dx <- function(x) {
  x <- trimws(as.character(x))
  y <- x
  y[grepl("^(asd|autism|case)$", x, ignore.case = TRUE) | grepl("asd|autism", x, ignore.case = TRUE)] <- "ASD"
  y[grepl("^(ctl|ctrl|control|healthy)$", x, ignore.case = TRUE) | grepl("control|ctl|ctrl|healthy", x, ignore.case = TRUE)] <- "Control"
  y
}

get_first_available_reduction <- function(obj, preferred = c("umap", "wnn.umap", "tsne", "pca")) {
  reds <- names(obj@reductions)
  for (r in preferred) {
    if (r %in% reds) return(r)
  }
  if (length(reds) == 0) stop("No dimensional reduction found in object.")
  reds[1]
}

save_plot_multi <- function(p, out_base, width, height, dpi = 300) {
  ggsave(paste0(out_base, ".pdf"), p, width = width, height = height, useDingbats = FALSE)
  ggsave(paste0(out_base, ".png"), p, width = width, height = height, dpi = dpi)
}

ensure_gene_present <- function(obj, gene) {
  assays <- Assays(obj)
  found <- vapply(assays, function(a) gene %in% rownames(obj[[a]]), logical(1))
  if (!any(found)) {
    stop(sprintf("Gene %s not found in any assay of object.", gene))
  }
  assay <- assays[which(found)[1]]
  DefaultAssay(obj) <- assay
  obj
}

extract_expr <- function(obj, gene) {
  vals <- FetchData(obj, vars = gene)
  vals[[gene]]
}

# -----------------------------
# setup
# -----------------------------
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$outdir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$outdir, "plots", "standalone"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$outdir, "tables"), recursive = TRUE, showWarnings = FALSE)

log_msg("Reading atlas object: %s", opt$atlas_rds)
atlas <- readRDS(opt$atlas_rds)
log_msg("Reading retained-core object: %s", opt$retained_rds)
retained <- readRDS(opt$retained_rds)

if (!inherits(atlas, "Seurat")) stop("atlas_rds is not a Seurat object")
if (!inherits(retained, "Seurat")) stop("retained_rds is not a Seurat object")

atlas <- ensure_gene_present(atlas, opt$gene)
retained <- ensure_gene_present(retained, opt$gene)

atlas_meta <- atlas@meta.data
ret_meta <- retained@meta.data

# atlas columns are now fixed by the inspected object, but still checked for safety
for (nm in c(opt$atlas_celltype_col, opt$atlas_dx_col, opt$atlas_donor_col)) {
  if (!nm %in% colnames(atlas_meta)) {
    stop(sprintf("Atlas metadata column not found: %s", nm))
  }
}

ret_dx_col <- if (!is.null(opt$ret_dx_col)) opt$ret_dx_col else choose_dx_col(ret_meta)
ret_donor_col <- if (!is.null(opt$ret_donor_col)) opt$ret_donor_col else choose_donor_col(ret_meta)
ret_cluster_col <- if (!is.null(opt$ret_cluster_col)) opt$ret_cluster_col else choose_cluster_col(ret_meta)

if (is.null(ret_dx_col) || !ret_dx_col %in% colnames(ret_meta)) {
  stop("Could not identify retained diagnosis column. Please set --ret_dx_col.")
}
if (is.null(ret_donor_col) || !ret_donor_col %in% colnames(ret_meta)) {
  stop("Could not identify retained donor column. Please set --ret_donor_col.")
}
if (is.null(ret_cluster_col) || !ret_cluster_col %in% colnames(ret_meta)) {
  stop("Could not identify retained cluster column. Please set --ret_cluster_col.")
}

log_msg("Atlas columns: celltype=%s, dx=%s, donor=%s", opt$atlas_celltype_col, opt$atlas_dx_col, opt$atlas_donor_col)
log_msg("Retained columns: dx=%s, donor=%s, cluster=%s", ret_dx_col, ret_donor_col, ret_cluster_col)

atlas$celltype_use <- as.character(atlas_meta[[opt$atlas_celltype_col]])
atlas$dx_use <- normalize_dx(atlas_meta[[opt$atlas_dx_col]])
atlas$donor_use <- as.character(atlas_meta[[opt$atlas_donor_col]])

retained$dx_use <- normalize_dx(ret_meta[[ret_dx_col]])
retained$donor_use <- as.character(ret_meta[[ret_donor_col]])
retained$cluster_use <- as.character(ret_meta[[ret_cluster_col]])

# save metadata snapshots
fwrite(data.table(column = colnames(atlas_meta)), file.path(opt$outdir, "tables", "atlas_metadata_columns.tsv"), sep = "\t")
fwrite(data.table(column = colnames(ret_meta)), file.path(opt$outdir, "tables", "retained_metadata_columns.tsv"), sep = "\t")

# -----------------------------
# atlas summaries
# -----------------------------
atlas_reduction <- get_first_available_reduction(atlas)
ret_reduction <- get_first_available_reduction(retained)
log_msg("Atlas reduction: %s", atlas_reduction)
log_msg("Retained reduction: %s", ret_reduction)

atlas_expr <- extract_expr(atlas, opt$gene)
atlas_df <- data.table(
  donor = atlas$donor_use,
  dx = atlas$dx_use,
  celltype = atlas$celltype_use,
  expr = atlas_expr
)

atlas_df <- atlas_df[!is.na(donor) & !is.na(dx) & !is.na(celltype)]

donor_ct_df <- atlas_df[, .(
  mean_expr = mean(expr, na.rm = TRUE),
  pct_expr = mean(expr > 0, na.rm = TRUE),
  n_cells = .N
), by = .(donor, dx, celltype)]

donor_ct_df <- donor_ct_df[n_cells >= opt$min_cells_per_donor_celltype]
fwrite(donor_ct_df, file.path(opt$outdir, "tables", "atlas_donor_celltype_summary.tsv"), sep = "\t")

celltype_order <- donor_ct_df[, .N, by = celltype][order(-N), celltype]
if (length(celltype_order) > 0) donor_ct_df[, celltype := factor(celltype, levels = celltype_order)]

# Panel A1: atlas UMAP by major class
pA1 <- DimPlot(
  atlas,
  reduction = atlas_reduction,
  group.by = "celltype_use",
  label = TRUE,
  repel = TRUE,
  raster = FALSE
) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0),
    legend.position = "right"
  ) +
  labs(title = "A. Major cortical cell classes")

# Panel A2: FeaturePlot for SPP1
pA2 <- FeaturePlot(
  atlas,
  features = opt$gene,
  reduction = atlas_reduction,
  raster = FALSE,
  order = TRUE,
  max.cutoff = "q99.5"
) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0)) +
  labs(title = sprintf("B. %s expression across the full cortical atlas", opt$gene))

# Panel B: dotplot across classes
pB <- DotPlot(atlas, features = opt$gene, group.by = "celltype_use") +
  RotatedAxis() +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0),
    axis.title = element_blank()
  ) +
  labs(title = sprintf("C. %s across major cell classes", opt$gene))

# Panel C: donor-level mean expression by dx across cell classes
pC <- ggplot(donor_ct_df, aes(x = dx, y = mean_expr, fill = dx)) +
  geom_boxplot(outlier.shape = NA, width = 0.62, alpha = 0.9) +
  geom_jitter(width = 0.12, size = 0.8, alpha = 0.75) +
  facet_wrap(~celltype, scales = "free_y", ncol = 4) +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 8),
    axis.text.x = element_text(angle = 20, hjust = 1),
    plot.title = element_text(face = "bold", hjust = 0)
  ) +
  labs(
    title = sprintf("D. Donor-level %s expression by diagnosis across major cell classes", opt$gene),
    x = NULL,
    y = "Donor-level mean expression"
  )

# -----------------------------
# retained summaries
# -----------------------------
ret_expr <- extract_expr(retained, opt$gene)
ret_df <- data.table(
  donor = retained$donor_use,
  dx = retained$dx_use,
  cluster = retained$cluster_use,
  expr = ret_expr
)
ret_df <- ret_df[!is.na(donor) & !is.na(dx) & !is.na(cluster)]

ret_donor_df <- ret_df[, .(
  mean_expr = mean(expr, na.rm = TRUE),
  pct_expr = mean(expr > 0, na.rm = TRUE),
  n_cells = .N
), by = .(donor, dx, cluster)]
ret_donor_df <- ret_donor_df[n_cells >= opt$min_cells_per_donor_cluster]

# try to order cluster 0 and 2 first if present
clus_vals <- unique(as.character(ret_donor_df$cluster))
preferred <- c("0", "2")
clus_levels <- c(preferred[preferred %in% clus_vals], setdiff(sort(clus_vals), preferred))
ret_donor_df[, cluster := factor(cluster, levels = clus_levels)]
fwrite(ret_donor_df, file.path(opt$outdir, "tables", "retained_donor_cluster_summary.tsv"), sep = "\t")

# retained UMAP optional standalone
pD1 <- DimPlot(
  retained,
  reduction = ret_reduction,
  group.by = "cluster_use",
  label = TRUE,
  repel = TRUE,
  raster = FALSE
) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0)) +
  labs(title = "E. Retained microglial core clusters")

pD2 <- ggplot(ret_donor_df, aes(x = cluster, y = mean_expr, fill = dx)) +
  geom_boxplot(outlier.shape = NA, position = position_dodge(width = 0.75), width = 0.65) +
  geom_point(position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75), size = 1.2, alpha = 0.8) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0)) +
  labs(
    title = sprintf("F. Donor-level %s expression within the retained microglial core", opt$gene),
    x = "Retained cluster",
    y = "Donor-level mean expression"
  )

# -----------------------------
# save standalone panels
# -----------------------------
standalone_dir <- file.path(opt$outdir, "plots", "standalone")
save_plot_multi(pA1, file.path(standalone_dir, "FigureS6_panel_A_major_cell_classes"), width = 7.5, height = 6)
save_plot_multi(pA2, file.path(standalone_dir, "FigureS6_panel_B_gene_featureplot"), width = 7.5, height = 6)
save_plot_multi(pB, file.path(standalone_dir, "FigureS6_panel_C_dotplot_cellclasses"), width = 5.5, height = 4.5)
save_plot_multi(pC, file.path(standalone_dir, "FigureS6_panel_D_donor_celltype_dx"), width = 12, height = 8)
save_plot_multi(pD1, file.path(standalone_dir, "FigureS6_panel_E_retained_umap"), width = 6, height = 5)
save_plot_multi(pD2, file.path(standalone_dir, "FigureS6_panel_F_retained_donor_cluster"), width = 7, height = 5)

# combined figure
fig <- (pA1 | pA2) / (pB | pD1) / (pC / pD2 + plot_layout(heights = c(1.6, 1))) +
  plot_annotation(
    title = sprintf("Supplementary Figure S6. %s expression across major cortical cell classes and within the retained microglial core", opt$gene),
    theme = theme(plot.title = element_text(face = "bold", size = 13, hjust = 0))
  )

save_plot_multi(fig, file.path(opt$outdir, "plots", "Supplementary_Figure_S6"), width = 16, height = 18)

# write quick run summary
summary_lines <- c(
  sprintf("gene\t%s", opt$gene),
  sprintf("atlas_rds\t%s", opt$atlas_rds),
  sprintf("retained_rds\t%s", opt$retained_rds),
  sprintf("atlas_celltype_col\t%s", opt$atlas_celltype_col),
  sprintf("atlas_dx_col\t%s", opt$atlas_dx_col),
  sprintf("atlas_donor_col\t%s", opt$atlas_donor_col),
  sprintf("ret_dx_col\t%s", ret_dx_col),
  sprintf("ret_donor_col\t%s", ret_donor_col),
  sprintf("ret_cluster_col\t%s", ret_cluster_col),
  sprintf("atlas_reduction\t%s", atlas_reduction),
  sprintf("ret_reduction\t%s", ret_reduction),
  sprintf("atlas_donor_celltype_rows\t%s", nrow(donor_ct_df)),
  sprintf("retained_donor_cluster_rows\t%s", nrow(ret_donor_df))
)
writeLines(summary_lines, file.path(opt$outdir, "tables", "FigureS6_run_summary.tsv"))

log_msg("Done. Outputs written to: %s", opt$outdir)
