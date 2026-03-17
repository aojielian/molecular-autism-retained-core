#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(Matrix)
  library(patchwork)
})

options(stringsAsFactors = FALSE)

cfg <- list(
  discovery_rds = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/PsychENCODE_Science_2024/Step3_Microglia/Microglia_Clustered.rds",
  strict_rds    = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Package2_strict_input_fixed_from_original.rds",
  audit_tsv     = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Package1b_ClusterAnnotation_ContaminantAudit/tables/22_cluster_annotation_summary_for_manuscript.tsv",
  outdir        = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Figure1_Final_NoWorkflow_v5",
  discovery_assay = "RNA",
  strict_assay    = "RNA",
  cluster_col   = NULL,
  donor_col     = NULL,
  dx_col        = NULL,
  strict_cluster_col = NULL,
  strict_donor_col   = NULL,
  strict_dx_col      = NULL,
  cand_cluster  = "2",
  ref_cluster   = "0"
)

dir.create(cfg$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cfg$outdir, "panels"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cfg$outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
logfile <- file.path(cfg$outdir, "00_Figure1_build.log")

msg <- function(...) {
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", sprintf(...))
  cat(line, "\n")
  cat(line, "\n", file = logfile, append = TRUE)
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
  x <- trimws(as.character(x))
  out <- x
  out[grepl("asd|autism|case", x, ignore.case = TRUE)] <- "ASD"
  out[grepl("control|ctrl|ctl|healthy|non[- ]?psy", x, ignore.case = TRUE)] <- "Control"
  out
}

choose_cluster_col <- function(obj, user_col = NULL) {
  meta <- obj@meta.data
  if (!is.null(user_col) && user_col %in% names(meta)) return(user_col)
  if ("cluster_use" %in% names(meta)) return("cluster_use")
  if ("seurat_clusters" %in% names(meta)) return("seurat_clusters")
  cand <- find_first_col(meta, c("^cluster$", "cluster", "seurat", "res\\.", "RNA_snn", "wsnn"))
  if (!is.null(cand)) return(cand)
  obj$cluster_auto_from_idents <- as.character(Idents(obj))
  "cluster_auto_from_idents"
}

choose_donor_col <- function(meta, user_col = NULL) {
  if (!is.null(user_col) && user_col %in% names(meta)) return(user_col)
  find_first_col(meta, c("^donor$", "donor_id", "individual", "subject", "patient", "sample", "orig.ident", "library", "case_id"))
}

choose_dx_col <- function(meta, user_col = NULL) {
  if (!is.null(user_col) && user_col %in% names(meta)) return(user_col)
  find_first_col(meta, c("^dx$", "^group$", "diagnosis", "condition", "status", "phenotype", "disease"))
}

safe_write <- function(x, file) fwrite(as.data.table(x), file, sep = "\t")

safe_scale01 <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  rng <- range(x, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0) return(rep(0.5, length(x)))
  (x - rng[1]) / (rng[2] - rng[1])
}

row_zscore_by_gene <- function(mat) {
  z <- t(apply(mat, 1, function(v) {
    if (all(is.na(v))) return(rep(NA_real_, length(v)))
    s <- sd(v, na.rm = TRUE)
    if (!is.finite(s) || s == 0) return(rep(0, length(v)))
    as.numeric(scale(v))
  }))
  rownames(z) <- rownames(mat)
  colnames(z) <- colnames(mat)
  z
}

make_cluster_levels <- function(v) {
  lev <- unique(as.character(v))
  num_try <- suppressWarnings(as.numeric(lev))
  if (all(!is.na(num_try))) as.character(sort(num_try)) else sort(lev)
}

# marker panel for panel B
marker_groups <- list(
  microglial = c("P2RY12","CX3CR1","TMEM119","MERTK","C3"),
  neuronal   = c("SYT1","NRXN3","RBFOX1","GPM6A","GRM5"),
  astrocytic = c("GFAP","AQP4","SLC1A3"),
  oligodendroglial = c("MBP","MOG","PLP1"),
  opc = c("PDGFRA","CSPG4"),
  endothelial_lymphoid = c("CLDN5","FLT1","PECAM1","NKG7","CD3D")
)
lineage_marker_panel <- unique(unlist(marker_groups))

# used for fallback audit scoring if audit_tsv missing
panel_def <- list(
  microglia = list(
    anchor  = c("P2RY12","CX3CR1","TMEM119","CSF1R","AIF1","C1QA","C1QB","C1QC","TREM2","TYROBP"),
    support = c("MERTK","FCER1G","LAPTM5","SPP1","GAS6","APOE","CD74","C3","LGALS3","GPNMB","LPL","FABP5","CTSB","CTSD","LIPA","ABCA1","MSR1")
  ),
  neuronal = list(anchor = c("RBFOX3","SNAP25","SYT1","STMN2"), support = c("SLC17A7","GAD1","GAD2","DLG4","GRIN2B","CAMK2A")),
  astrocyte = list(anchor = c("GFAP","AQP4","ALDH1L1"), support = c("SLC1A2","SLC1A3","GJA1","ALDOC")),
  oligodendrocyte = list(anchor = c("MBP","MOG","PLP1"), support = c("MOBP","MAG","CNP")),
  opc = list(anchor = c("PDGFRA","VCAN","PTPRZ1"), support = c("CSPG4","OLIG1","OLIG2")),
  endothelial = list(anchor = c("CLDN5","FLT1","KDR","VWF","EPAS1"), support = c("ESAM","RAMP2","PLVAP","EMCN")),
  t_cell = list(anchor = c("CD3D","CD3E","CD247","TRBC1","TRBC2","THEMIS"), support = c("IL7R","LTB","ETS1")),
  b_cell = list(anchor = c("MS4A1","CD79A","CD79B","IGHM","PAX5"), support = c("BLK","FCRL1","CD22","BANK1")),
  pericyte = list(anchor = c("PDGFRB","RGS5"), support = c("MCAM","COL1A1","COL1A2","TAGLN")),
  erythroid = list(anchor = c("HBB","HBA1","HBA2"), support = c("ALAS2","AHSP"))
)
nonmg_lineages <- setdiff(names(panel_def), "microglia")

get_panel_scores <- function(mat, genes_anchor, genes_support) {
  anchor_use <- intersect(genes_anchor, rownames(mat))
  support_use <- intersect(genes_support, rownames(mat))
  anchor_score <- if (length(anchor_use) > 0) colMeans(row_zscore_by_gene(mat[anchor_use, , drop = FALSE]), na.rm = TRUE) else setNames(rep(NA_real_, ncol(mat)), colnames(mat))
  support_score <- if (length(support_use) > 0) colMeans(row_zscore_by_gene(mat[support_use, , drop = FALSE]), na.rm = TRUE) else setNames(rep(NA_real_, ncol(mat)), colnames(mat))
  total_score <- 0.75 * anchor_score + 0.25 * support_score
  list(anchor_score = anchor_score, support_score = support_score, total_score = total_score)
}

build_fallback_audit <- function(obj, cluster_levels, donor_col) {
  msg("Building fallback audit summary from discovery object ...")
  avg_list <- AverageExpression(obj, assays = DefaultAssay(obj), group.by = "cluster_use", slot = "data", return.seurat = FALSE, verbose = FALSE)
  avg_mat <- avg_list[[1]]
  avg_mat <- as.matrix(avg_mat[, cluster_levels, drop = FALSE])

  score_dt <- data.table(cluster = cluster_levels)
  for (lin in names(panel_def)) {
    sc <- get_panel_scores(avg_mat, panel_def[[lin]]$anchor, panel_def[[lin]]$support)
    score_dt[[paste0(lin, "_total_score")]] <- as.numeric(sc$total_score[cluster_levels])
  }

  donor_rep <- as.data.table(obj@meta.data, keep.rownames = "cell")
  donor_rep[, cluster_use := as.character(obj$cluster_use)]
  donor_rep[, donor_use := as.character(get(donor_col))]
  donor_dt <- donor_rep[, .(n_donors = uniqueN(donor_use)), by = cluster_use]
  setnames(donor_dt, "cluster_use", "cluster")

  out <- merge(score_dt, donor_dt, by = "cluster", all.x = TRUE)
  out
}

# -------------------------------------------------
# Read discovery object
# -------------------------------------------------
msg("Reading discovery object ...")
disc_obj <- readRDS(cfg$discovery_rds)
stopifnot(inherits(disc_obj, "Seurat"))
if (cfg$discovery_assay %in% Assays(disc_obj)) DefaultAssay(disc_obj) <- cfg$discovery_assay

disc_cluster_col <- choose_cluster_col(disc_obj, cfg$cluster_col)
disc_obj$cluster_use <- as.character(disc_obj@meta.data[[disc_cluster_col]])
Idents(disc_obj) <- disc_obj$cluster_use

disc_donor_col <- choose_donor_col(disc_obj@meta.data, cfg$donor_col)
disc_dx_col    <- choose_dx_col(disc_obj@meta.data, cfg$dx_col)
if (is.null(disc_donor_col)) stop("Could not identify donor column in discovery object.")
if (is.null(disc_dx_col)) stop("Could not identify diagnosis column in discovery object.")

disc_obj$dx_use <- normalize_dx(disc_obj@meta.data[[disc_dx_col]])
cluster_levels <- make_cluster_levels(disc_obj$cluster_use)
disc_obj$cluster_use <- factor(disc_obj$cluster_use, levels = cluster_levels)
Idents(disc_obj) <- disc_obj$cluster_use
safe_write(data.frame(cluster = cluster_levels), file.path(cfg$outdir, "tables", "01_discovery_clusters.tsv"))

# -------------------------------------------------
# Read strict retained object and derive retained clusters from REAL data
# -------------------------------------------------
msg("Reading strict retained object ...")
strict_obj <- readRDS(cfg$strict_rds)
stopifnot(inherits(strict_obj, "Seurat"))
if (cfg$strict_assay %in% Assays(strict_obj)) DefaultAssay(strict_obj) <- cfg$strict_assay
strict_cluster_col <- choose_cluster_col(strict_obj, cfg$strict_cluster_col)
strict_obj$cluster_use <- as.character(strict_obj@meta.data[[strict_cluster_col]])
Idents(strict_obj) <- strict_obj$cluster_use
strict_retained_clusters <- sort(unique(as.character(strict_obj$cluster_use)))
msg("Retained clusters inferred from strict object: %s", paste(strict_retained_clusters, collapse = ", "))
safe_write(data.frame(retained_cluster = strict_retained_clusters), file.path(cfg$outdir, "tables", "02_retained_clusters_from_strict_object.tsv"))

strict_donor_col <- choose_donor_col(strict_obj@meta.data, cfg$strict_donor_col)
strict_dx_col    <- choose_dx_col(strict_obj@meta.data, cfg$strict_dx_col)
if (is.null(strict_donor_col)) strict_donor_col <- disc_donor_col
if (is.null(strict_dx_col)) strict_dx_col <- disc_dx_col
strict_obj$dx_use <- normalize_dx(strict_obj@meta.data[[strict_dx_col]])

# -------------------------------------------------
# Panel A: all-cluster UMAP
# -------------------------------------------------
msg("Building Panel A ...")
if (!("umap" %in% names(disc_obj@reductions))) stop("Discovery object has no UMAP reduction.")
pal_all <- setNames(grDevices::hcl.colors(length(cluster_levels), "Dynamic"), cluster_levels)
p_A <- DimPlot(
  disc_obj,
  reduction = "umap",
  group.by = "cluster_use",
  label = TRUE,
  repel = TRUE,
  raster = FALSE
) +
  scale_color_manual(values = pal_all, drop = FALSE) +
  theme_classic(base_size = 12) +
  labs(title = "A  Reclustered discovery microglial object") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0),
    legend.position = "right"
  )
ggsave(file.path(cfg$outdir, "panels", "Figure1A_allcluster_umap.pdf"), p_A, width = 8.4, height = 6.8)
ggsave(file.path(cfg$outdir, "panels", "Figure1A_allcluster_umap.png"), p_A, width = 8.4, height = 6.8, dpi = 300)

# -------------------------------------------------
# Panel B: lineage/background dotplot
# -------------------------------------------------
msg("Building Panel B ...")
marker_use <- intersect(lineage_marker_panel, rownames(disc_obj))
if (length(marker_use) == 0) stop("No marker genes from the lineage/background panel were found in the discovery object.")
p_B <- DotPlot(disc_obj, features = marker_use, group.by = "cluster_use") +
  RotatedAxis() +
  theme_classic(base_size = 10) +
  labs(title = "B  Transparent lineage/background marker panel") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )
ggsave(file.path(cfg$outdir, "panels", "Figure1B_lineage_dotplot.pdf"), p_B, width = 13.5, height = 5.8)
ggsave(file.path(cfg$outdir, "panels", "Figure1B_lineage_dotplot.png"), p_B, width = 13.5, height = 5.8, dpi = 300)

# -------------------------------------------------
# Panel C: audit summary using real audit data + retained derived from strict object
# -------------------------------------------------
msg("Building Panel C ...")
if (!is.na(cfg$audit_tsv) && file.exists(cfg$audit_tsv)) {
  audit_dt <- fread(cfg$audit_tsv)
} else {
  audit_dt <- build_fallback_audit(disc_obj, cluster_levels, disc_donor_col)
}
audit_dt <- as.data.table(audit_dt)

# harmonize cluster column
if (!("cluster" %in% names(audit_dt))) {
  cand_cluster_col <- find_first_col(audit_dt, c("^cluster$", "cluster"))
  if (is.null(cand_cluster_col)) stop("Could not identify cluster column in audit summary table.")
  setnames(audit_dt, cand_cluster_col, "cluster")
}
audit_dt[, cluster := as.character(cluster)]

# ensure all discovery clusters present
base_dt <- data.table(cluster = cluster_levels)
audit_dt <- merge(base_dt, audit_dt, by = "cluster", all.x = TRUE)

# donor representation from real discovery object if absent
if (!("n_donors" %in% names(audit_dt))) {
  donor_dt <- as.data.table(disc_obj@meta.data, keep.rownames = "cell")
  donor_dt[, cluster_use := as.character(disc_obj$cluster_use)]
  donor_dt[, donor_use := as.character(get(disc_donor_col))]
  donor_dt <- donor_dt[, .(n_donors = uniqueN(donor_use)), by = cluster_use]
  setnames(donor_dt, "cluster_use", "cluster")
  audit_dt <- merge(audit_dt, donor_dt, by = "cluster", all.x = TRUE)
}

# derive microglial / non-microglial scores from real audit table if available; otherwise fallback from discovery object
need_fallback <- !("microglia_total_score" %in% names(audit_dt))
if (need_fallback) {
  fb <- build_fallback_audit(disc_obj, cluster_levels, disc_donor_col)
  audit_dt <- merge(audit_dt, fb, by = c("cluster", "n_donors"), all.x = TRUE, suffixes = c("", ".fb"))
}

# real retained decision MUST come from strict object, not hard-coded panel logic
retained_real <- cluster_levels %in% strict_retained_clusters
audit_dt[, retained_final := retained_real]

audit_dt[, microglial_identity := safe_scale01(microglia_total_score)]

# compute top non-microglial background from available real audit columns
nonmg_cols_present <- intersect(paste0(nonmg_lineages, "_total_score"), names(audit_dt))
if (length(nonmg_cols_present) > 0) {
  audit_dt[, top_nonmg_score := do.call(pmax, c(.SD, list(na.rm = TRUE))), .SDcols = nonmg_cols_present]
} else {
  audit_dt[, top_nonmg_score := NA_real_]
}
audit_dt[, non_microglial_background := safe_scale01(top_nonmg_score)]
audit_dt[, donor_representation := safe_scale01(as.numeric(n_donors))]
audit_dt[, retained_numeric := ifelse(retained_final, 1, 0)]

panelC_long <- rbindlist(list(
  audit_dt[, .(cluster, metric = "Microglial\nidentity", value = microglial_identity)],
  audit_dt[, .(cluster, metric = "Non-microglial\nbackground", value = non_microglial_background)],
  audit_dt[, .(cluster, metric = "Donor\nrepresentation", value = donor_representation)],
  audit_dt[, .(cluster, metric = "Retained", value = retained_numeric)]
), fill = TRUE)
panelC_long[, cluster := factor(cluster, levels = rev(cluster_levels))]
panelC_long[, metric := factor(metric, levels = c("Microglial\nidentity", "Non-microglial\nbackground", "Donor\nrepresentation", "Retained"))]

label_dt <- unique(panelC_long[, .(cluster)])
label_dt[, x_ret := "Retained"]
label_dt[, retained_label := ifelse(as.character(cluster) %in% strict_retained_clusters, "Yes", "No")]

panelC_cont <- panelC_long[metric != "Retained"]
panelC_ret_yes <- merge(panelC_long[metric == "Retained"], label_dt[retained_label == "Yes", .(cluster)], by = "cluster")
panelC_ret_no  <- merge(panelC_long[metric == "Retained"], label_dt[retained_label == "No",  .(cluster)], by = "cluster")

p_C <- ggplot() +
  geom_tile(data = panelC_cont, aes(x = metric, y = cluster, fill = value), color = "white", linewidth = 0.5) +
  scale_fill_gradient(low = "#f7f7f7", high = "#2166ac", limits = c(0, 1), na.value = "#f0f0f0") +
  geom_tile(data = panelC_ret_no, aes(x = metric, y = cluster), fill = "#E9E9E9", color = "white", linewidth = 0.5) +
  geom_tile(data = panelC_ret_yes, aes(x = metric, y = cluster), fill = "#2166ac", color = "white", linewidth = 0.5) +
  geom_text(data = label_dt, aes(x = x_ret, y = cluster, label = retained_label, color = retained_label), inherit.aes = FALSE, size = 3.4, fontface = "bold") +
  scale_color_manual(values = c("Yes" = "white", "No" = "black"), guide = "none") +
  geom_vline(xintercept = 3.5, linewidth = 0.8, color = "#636363") +
  theme_classic(base_size = 11) +
  labs(title = "C  Audit summary and retention decision", x = NULL, y = NULL) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0),
    axis.text.x = element_text(face = "bold"),
    legend.position = "none"
  )
ggsave(file.path(cfg$outdir, "panels", "Figure1C_audit_summary.pdf"), p_C, width = 6.6, height = 6.8)
ggsave(file.path(cfg$outdir, "panels", "Figure1C_audit_summary.png"), p_C, width = 6.6, height = 6.8, dpi = 300)
safe_write(audit_dt[, .(cluster, microglia_total_score, top_nonmg_score, n_donors, retained_final)], file.path(cfg$outdir, "tables", "03_panelC_underlying_values.tsv"))

# -------------------------------------------------
# Panel D: recomputed UMAP on strict retained object (clean, only clusters 0 and 2)
# -------------------------------------------------
msg("Building Panel D ...")
strict_keep <- unique(as.character(strict_obj$cluster_use))
strict_obj <- subset(strict_obj, cells = colnames(strict_obj)[as.character(strict_obj$cluster_use) %in% strict_keep])
strict_obj$cluster_use <- factor(as.character(strict_obj$cluster_use), levels = make_cluster_levels(strict_keep))
Idents(strict_obj) <- strict_obj$cluster_use

if (!inherits(tryCatch(GetAssayData(strict_obj, assay = DefaultAssay(strict_obj), slot = "data"), error = function(e) NULL), "Matrix")) {
  strict_obj <- NormalizeData(strict_obj, assay = DefaultAssay(strict_obj), verbose = FALSE)
}
# recompute a clean retained-object UMAP
strict_obj <- FindVariableFeatures(strict_obj, assay = DefaultAssay(strict_obj), selection.method = "vst", nfeatures = min(2000, nrow(strict_obj)), verbose = FALSE)
strict_obj <- ScaleData(strict_obj, assay = DefaultAssay(strict_obj), features = VariableFeatures(strict_obj), verbose = FALSE)
strict_obj <- RunPCA(strict_obj, assay = DefaultAssay(strict_obj), features = VariableFeatures(strict_obj), npcs = min(20, max(5, length(VariableFeatures(strict_obj)) - 1)), verbose = FALSE)
use_dims <- 1:min(15, ncol(Embeddings(strict_obj, "pca")))
strict_obj <- FindNeighbors(strict_obj, dims = use_dims, verbose = FALSE)
strict_obj <- RunUMAP(strict_obj, dims = use_dims, reduction = "pca", reduction.name = "umap_clean", umap.method = "uwot", metric = "cosine", verbose = FALSE)

pal_strict <- setNames(c("#4C86C6", "#E87547"), c(cfg$ref_cluster, cfg$cand_cluster))
pal_strict <- pal_strict[names(pal_strict) %in% levels(strict_obj$cluster_use)]
p_D <- DimPlot(
  strict_obj,
  reduction = "umap_clean",
  group.by = "cluster_use",
  label = TRUE,
  repel = TRUE,
  raster = FALSE
) +
  scale_color_manual(values = pal_strict, drop = FALSE) +
  theme_classic(base_size = 12) +
  labs(title = "D  Strict retained core microglia") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0),
    legend.position = "right"
  )
ggsave(file.path(cfg$outdir, "panels", "Figure1D_retainedcore_umap.pdf"), p_D, width = 7.0, height = 5.8)
ggsave(file.path(cfg$outdir, "panels", "Figure1D_retainedcore_umap.png"), p_D, width = 7.0, height = 5.8, dpi = 300)

# -------------------------------------------------
# Panel E: donor representation heatmap from strict object
# -------------------------------------------------
msg("Building Panel E ...")
strict_meta <- as.data.table(strict_obj@meta.data, keep.rownames = "cell")
strict_meta[, donor_use := as.character(get(strict_donor_col))]
strict_meta[, dx_use := normalize_dx(get(strict_dx_col))]
strict_meta[, cluster_use := as.character(cluster_use)]
strict_meta <- strict_meta[dx_use %in% c("Control", "ASD") & !is.na(donor_use) & donor_use != ""]

# within-donor proportions across strict retained cells only
count_dt <- strict_meta[, .N, by = .(donor_use, dx_use, cluster_use)]
count_dt[, donor_total := sum(N), by = donor_use]
count_dt[, within_donor_prop := N / donor_total]

# donor ordering by dx then candidate proportion then reference proportion
wide_dt <- dcast(count_dt, donor_use + dx_use ~ cluster_use, value.var = "within_donor_prop", fill = 0)
if (!(cfg$ref_cluster %in% names(wide_dt))) wide_dt[[cfg$ref_cluster]] <- 0
if (!(cfg$cand_cluster %in% names(wide_dt))) wide_dt[[cfg$cand_cluster]] <- 0
setnames(wide_dt, old = cfg$ref_cluster, new = "ref_prop")
setnames(wide_dt, old = cfg$cand_cluster, new = "cand_prop")
wide_dt[, dx_order := ifelse(dx_use == "Control", 0L, 1L)]
setorder(wide_dt, dx_order, -cand_prop, -ref_prop, donor_use)
donor_levels <- wide_dt$donor_use

count_dt[, donor_use := factor(donor_use, levels = donor_levels)]
count_dt[, cluster_use := factor(cluster_use, levels = c(cfg$ref_cluster, cfg$cand_cluster))]
count_dt[, dx_use := factor(dx_use, levels = c("Control", "ASD"))]

p_E <- ggplot(count_dt, aes(x = donor_use, y = cluster_use, fill = within_donor_prop)) +
  geom_tile(color = "white", linewidth = 0.2) +
  facet_grid(. ~ dx_use, scales = "free_x", space = "free_x") +
  scale_fill_gradient(low = "#f7fbff", high = "#08519c", limits = c(0, 1)) +
  theme_classic(base_size = 11) +
  labs(title = "E  Donor-level representation of retained clusters", x = NULL, y = NULL, fill = "Within-donor\nproportion") +
  theme(
    plot.title = element_text(face = "bold", hjust = 0),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )
ggsave(file.path(cfg$outdir, "panels", "Figure1E_donor_representation.pdf"), p_E, width = 10.8, height = 3.6)
ggsave(file.path(cfg$outdir, "panels", "Figure1E_donor_representation.png"), p_E, width = 10.8, height = 3.6, dpi = 300)
safe_write(count_dt, file.path(cfg$outdir, "tables", "04_panelE_donor_representation.tsv"))

# -------------------------------------------------
# Combine A-E
# -------------------------------------------------
msg("Combining panels ...")
combined <- (p_A | p_B) / (p_C | p_D | p_E) + plot_layout(heights = c(1, 1), widths = c(1, 1, 1))
ggsave(file.path(cfg$outdir, "Figure1_A_to_E_combined.pdf"), combined, width = 18, height = 11)
ggsave(file.path(cfg$outdir, "Figure1_A_to_E_combined.png"), combined, width = 18, height = 11, dpi = 300)

msg("Done.")
