#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(Matrix)
  library(patchwork)
  library(scales)
})

options(stringsAsFactors = FALSE)

cfg <- list(
  # -------- inputs --------
  microglia_rds = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/PsychENCODE_Science_2024/Step3_Microglia/Microglia_Clustered.rds",
  package1_markers = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Package1_Figure1_DiscoveryAudit/03_allcluster_markers.tsv.gz",
  package1b_audit = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Package1b_ClusterAnnotation_ContaminantAudit/tables/10_cluster_annotation_contaminant_audit.tsv",
  package1b_strict_clusters = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Package1b_ClusterAnnotation_ContaminantAudit/tables/11_strict_microglia_clusters.txt",
  strict_rds = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Package2_strict_input_fixed_from_original.rds",
  assay = "RNA",

  # -------- figure config --------
  outdir = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Figure1_Final",
  retained_clusters = c("0", "2"),
  ref_cluster = "0",
  cand_cluster = "2",
  top_n_markers = 50,
  panelA_png = NA_character_,   # optional external workflow schematic PNG/PDF exported as image

  # label/size tuning
  seed = 1,
  point_size_umap = 0.5,
  donor_heatmap_fill_limits = c(0, 1),
  width_final = 16,
  height_final = 12
)

# ------------------------------------------------------------
# helpers
# ------------------------------------------------------------
msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ..., "\n", sep = "")
}

find_first_col <- function(df, patterns) {
  nm <- names(df)
  for (p in patterns) {
    hit <- nm[grepl(p, nm, ignore.case = TRUE)]
    if (length(hit) > 0) return(hit[1])
  }
  NULL
}

choose_cluster_col <- function(obj) {
  meta <- obj@meta.data
  if ("cluster_use" %in% names(meta)) return("cluster_use")
  if ("seurat_clusters" %in% names(meta)) return("seurat_clusters")
  cand <- find_first_col(meta, c("^cluster$", "cluster", "seurat", "res\\.", "RNA_snn", "wsnn"))
  if (!is.null(cand)) return(cand)
  obj$cluster_auto_from_idents <- as.character(Idents(obj))
  return("cluster_auto_from_idents")
}

choose_donor_col <- function(meta) {
  find_first_col(meta, c("^donor$", "donor_id", "individual", "subject", "patient", "sample", "orig.ident", "library", "case_id"))
}

choose_dx_col <- function(meta) {
  find_first_col(meta, c("^dx$", "^group$", "diagnosis", "condition", "status", "phenotype", "disease"))
}

normalize_dx <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  y <- x
  y[grepl("asd|autism|case", x, ignore.case = TRUE)] <- "ASD"
  y[grepl("control|ctrl|ctl|healthy|non[- ]?psy", x, ignore.case = TRUE)] <- "Control"
  y
}

safe_read <- function(path, fun) {
  if (!file.exists(path)) return(NULL)
  tryCatch(fun(path), error = function(e) NULL)
}

get_fc_col <- function(markers) {
  fc_candidates <- c("avg_log2FC", "avg_logFC", "log2FC", "logFC")
  hit <- fc_candidates[fc_candidates %in% colnames(markers)]
  if (length(hit) == 0) return(NULL)
  hit[1]
}

read_any_markers <- function(path) {
  x <- fread(path)
  if (!"gene" %in% names(x)) x[, gene := rownames(x)]
  fc_col <- get_fc_col(x)
  if (is.null(fc_col)) stop("Could not identify FC column in marker table.")
  x[, cluster := as.character(cluster)]
  x[]
}

row_zscore_by_gene <- function(mat) {
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

# Lineage panels adapted from package1b
panel_def <- list(
  microglia = list(
    anchor  = c("P2RY12","CX3CR1","TMEM119","CSF1R","AIF1","C1QA","C1QB","C1QC","TREM2","TYROBP"),
    support = c("MERTK","FCER1G","LAPTM5","SPP1","GAS6","APOE","CD74","C3","LGALS3","GPNMB","LPL","FABP5","CTSB","CTSD","LIPA","ABCA1","MSR1")
  ),
  t_cell = list(
    anchor  = c("CD3D","CD3E","CD247","TRBC1","TRBC2","THEMIS"),
    support = c("IL7R","LTB","ETS1")
  ),
  b_cell = list(
    anchor  = c("MS4A1","CD79A","CD79B","IGHM","PAX5"),
    support = c("BLK","FCRL1","CD22","BANK1")
  ),
  endothelial = list(
    anchor  = c("CLDN5","FLT1","KDR","VWF","EPAS1"),
    support = c("ESAM","RAMP2","PLVAP","EMCN")
  ),
  neuronal = list(
    anchor  = c("RBFOX3","SNAP25","SYT1","STMN2"),
    support = c("SLC17A7","GAD1","GAD2","DLG4","GRIN2B","CAMK2A")
  ),
  astrocyte = list(
    anchor  = c("GFAP","AQP4","ALDH1L1"),
    support = c("SLC1A2","SLC1A3","GJA1","ALDOC")
  ),
  oligodendrocyte = list(
    anchor  = c("MBP","MOG","PLP1"),
    support = c("MOBP","MAG","CNP")
  ),
  opc = list(
    anchor  = c("PDGFRA","VCAN","PTPRZ1"),
    support = c("CSPG4","OLIG1","OLIG2")
  ),
  pericyte = list(
    anchor  = c("PDGFRB","RGS5"),
    support = c("MCAM","COL1A1","COL1A2","TAGLN")
  ),
  erythroid = list(
    anchor  = c("HBB","HBA1","HBA2"),
    support = c("ALAS2","AHSP")
  )
)

lineage_names <- names(panel_def)
lineage_colors <- c(
  microglia = "#1b9e77",
  mixed_microglia = "#7570b3",
  uncertain = "#999999",
  contaminant = "#d95f02"
)

lineage_marker_panel <- unique(unlist(lapply(panel_def, function(x) c(x$anchor, x$support))))

# ------------------------------------------------------------
# read source objects and tables
# ------------------------------------------------------------
dir.create(cfg$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cfg$outdir, "panels"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cfg$outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cfg$outdir, "objects"), recursive = TRUE, showWarnings = FALSE)

msg("Reading discovery object ...")
obj <- readRDS(cfg$microglia_rds)
if (!inherits(obj, "Seurat")) stop("microglia_rds is not a Seurat object")
if (cfg$assay %in% Assays(obj)) DefaultAssay(obj) <- cfg$assay

cluster_col <- choose_cluster_col(obj)
donor_col <- choose_donor_col(obj@meta.data)
dx_col <- choose_dx_col(obj@meta.data)
if (is.null(donor_col)) stop("Could not identify donor column.")
if (is.null(dx_col)) stop("Could not identify dx column.")

obj$cluster_use <- as.character(obj@meta.data[[cluster_col]])
Idents(obj) <- obj$cluster_use
cluster_levels <- unique(as.character(obj$cluster_use))
num_try <- suppressWarnings(as.numeric(cluster_levels))
cluster_levels <- if (all(!is.na(num_try))) as.character(sort(num_try)) else sort(cluster_levels)
obj$cluster_use <- factor(obj$cluster_use, levels = cluster_levels)

markers <- safe_read(cfg$package1_markers, read_any_markers)
if (is.null(markers)) {
  msg("Package1 marker table not found; recomputing FindAllMarkers ...")
  markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.10, logfc.threshold = 0.25)
  if (!"gene" %in% names(markers)) markers$gene <- rownames(markers)
  markers[, cluster := as.character(cluster)]
}
fc_col <- get_fc_col(markers)
markers <- as.data.table(markers)

final_anno <- safe_read(cfg$package1b_audit, fread)
if (is.null(final_anno)) stop("package1b audit table not found: ", cfg$package1b_audit)
final_anno[, cluster := as.character(cluster)]

strict_clusters_file <- safe_read(cfg$package1b_strict_clusters, readLines)
retained_clusters <- if (!is.null(strict_clusters_file)) unique(trimws(strict_clusters_file)) else cfg$retained_clusters
if (!all(cfg$retained_clusters %in% retained_clusters)) retained_clusters <- cfg$retained_clusters

strict_obj <- safe_read(cfg$strict_rds, readRDS)
if (is.null(strict_obj)) {
  msg("strict_rds not found; creating retained-only object from discovery object.")
  strict_obj <- subset(obj, cells = colnames(obj)[as.character(obj$cluster_use) %in% retained_clusters])
}
if (cfg$assay %in% Assays(strict_obj)) DefaultAssay(strict_obj) <- cfg$assay
if (!"cluster_use" %in% colnames(strict_obj@meta.data)) {
  strict_obj$cluster_use <- as.character(strict_obj@meta.data[[choose_cluster_col(strict_obj)]])
}
strict_obj$cluster_use <- factor(as.character(strict_obj$cluster_use), levels = cfg$retained_clusters)
Idents(strict_obj) <- strict_obj$cluster_use

saveRDS(strict_obj, file.path(cfg$outdir, "objects", "Figure1_strict_object_used.rds"))

# ------------------------------------------------------------
# Panel A prompt only (saved as txt for provenance)
# ------------------------------------------------------------
prompt_A <- paste(
  "Create a clean scientific workflow schematic for a Molecular Autism research paper figure.",
  "Style: professional biomedical journal infographic, white background, minimal, vector-like, no 3D, no cartoon cells, no decorative icons.",
  "Use a left-to-right or top-to-bottom flow with subtle arrows.",
  "Main steps:",
  "1) Discovery cortical microglia object from human ASD and control cortex.",
  "2) All-cluster audit and contaminant-aware lineage assessment.",
  "3) Exclude contaminant, mixed, and background-driven clusters.",
  "4) Retain strict core microglial compartment: clusters 0 and 2 only.",
  "5) Donor-aware discovery analyses: state definition, abundance trend, paired state shift, pooled and cluster-specific pseudobulk DE.",
  "6) Cross-cohort single-cell validation in Velmeshev dataset using score-based candidate/reference metrics.",
  "7) External bulk-cortex validation in Gandal 2022 dataset using frozen candidate and reference programs.",
  "8) Leave-one-region-out sensitivity analysis.",
  "Visually emphasize retained clusters 0 and 2 as the central decision point.",
  "Use neutral gray for excluded clusters, blue for cluster 0, orange-red for cluster 2.",
  "Add short panel-style labels only if space allows; do not include paragraph text.",
  "Make the final figure compact, publication-ready, and suitable as Figure 1A in a high-quality autism transcriptomics paper."
)
writeLines(prompt_A, file.path(cfg$outdir, "Figure1A_nano_banana_prompt.txt"))

# ------------------------------------------------------------
# Panel B: all-cluster UMAP
# ------------------------------------------------------------
msg("Building Panel B ...")
set.seed(cfg$seed)
cluster_pal <- setNames(hue_pal()(length(cluster_levels)), cluster_levels)
cluster_pal[cfg$ref_cluster] <- "#3B6FB6"
cluster_pal[cfg$cand_cluster] <- "#D95F02"

pB <- DimPlot(
  obj,
  reduction = if ("umap" %in% names(obj@reductions)) "umap" else names(obj@reductions)[1],
  group.by = "cluster_use",
  label = TRUE,
  repel = TRUE,
  raster = FALSE,
  pt.size = cfg$point_size_umap
) +
  scale_color_manual(values = cluster_pal, drop = FALSE) +
  theme_classic(base_size = 12) +
  labs(title = "All reclustered discovery microglial clusters") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

ggsave(file.path(cfg$outdir, "panels", "Figure1B_allcluster_umap.pdf"), pB, width = 6.8, height = 5.8)
ggsave(file.path(cfg$outdir, "panels", "Figure1B_allcluster_umap.png"), pB, width = 6.8, height = 5.8, dpi = 300)

# ------------------------------------------------------------
# Panel C: lineage marker dotplot
# ------------------------------------------------------------
msg("Building Panel C ...")
features_C <- intersect(lineage_marker_panel, rownames(obj))
# keep genes in a stable grouped order and not too many
ordered_C <- unique(c(
  panel_def$microglia$anchor[1:6], panel_def$microglia$support[1:6],
  panel_def$neuronal$anchor[1:4], panel_def$neuronal$support[1:4],
  panel_def$astrocyte$anchor[1:3], panel_def$astrocyte$support[1:2],
  panel_def$oligodendrocyte$anchor[1:3], panel_def$opc$anchor[1:3],
  panel_def$endothelial$anchor[1:4], panel_def$t_cell$anchor[1:4], panel_def$b_cell$anchor[1:4]
))
features_C <- intersect(ordered_C, rownames(obj))

pC <- DotPlot(obj, features = features_C, group.by = "cluster_use") +
  RotatedAxis() +
  theme_classic(base_size = 10) +
  labs(title = "Transparent lineage/background marker panel across clusters") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 8)
  )

ggsave(file.path(cfg$outdir, "panels", "Figure1C_lineage_marker_dotplot.pdf"), pC, width = 13.5, height = 5.6)
ggsave(file.path(cfg$outdir, "panels", "Figure1C_lineage_marker_dotplot.png"), pC, width = 13.5, height = 5.6, dpi = 300)

# ------------------------------------------------------------
# Panel D: audit summary heatmap/table-like plot
# ------------------------------------------------------------
msg("Building Panel D ...")
# choose compact display columns
keep_cols <- c(
  "cluster", "n_cells", "frac_cells",
  "microglia_total_score", "microglia_anchor_score", "microglia_support_score",
  "neuronal_total_score", "astrocyte_total_score", "oligodendrocyte_total_score",
  "opc_total_score", "endothelial_total_score", "t_cell_total_score", "b_cell_total_score",
  "suggested_superclass", "suggested_lineage", "suggested_label", "keep_strict"
)
keep_cols <- intersect(keep_cols, colnames(final_anno))
audit_compact <- copy(final_anno)[, ..keep_cols]
# compute top non-mg score if missing
score_cols <- grep("_total_score$", names(final_anno), value = TRUE)
nonmg_score_cols <- setdiff(score_cols, "microglia_total_score")
if (length(nonmg_score_cols) > 0) {
  final_anno[, top_nonmg_score := do.call(pmax, c(.SD, na.rm = TRUE)), .SDcols = nonmg_score_cols]
  final_anno[, top_nonmg_lineage := nonmg_score_cols[max.col(as.matrix(.SD), ties.method = "first")], .SDcols = nonmg_score_cols]
  final_anno[, top_nonmg_lineage := sub("_total_score$", "", top_nonmg_lineage)]
} else {
  final_anno[, top_nonmg_score := NA_real_]
  final_anno[, top_nonmg_lineage := NA_character_]
}

D_long <- rbindlist(list(
  data.table(cluster = final_anno$cluster, metric = "Microglia score", value = final_anno$microglia_total_score),
  data.table(cluster = final_anno$cluster, metric = "Top non-microglial score", value = final_anno$top_nonmg_score),
  data.table(cluster = final_anno$cluster, metric = "Retained", value = as.numeric(final_anno$keep_strict %in% TRUE))
), fill = TRUE)
D_long[, cluster := factor(cluster, levels = rev(final_anno$cluster))]
D_long[, metric := factor(metric, levels = c("Microglia score", "Top non-microglial score", "Retained"))]

pD_heat <- ggplot(D_long, aes(x = metric, y = cluster, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#f7f7f7", mid = "#9ecae1", high = "#08519c", midpoint = 0.5, na.value = "grey90") +
  theme_minimal(base_size = 11) +
  labs(title = "Audit summary and retention decision", x = NULL, y = "Cluster") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid = element_blank(),
    axis.text.x = element_text(face = "bold")
  )

D_text <- final_anno[, .(
  cluster,
  superclass = suggested_superclass,
  lineage = suggested_lineage,
  keep = ifelse(keep_strict %in% TRUE, "Retained", "Excluded")
)]
D_text$cluster <- factor(D_text$cluster, levels = rev(final_anno$cluster))

pD_text <- ggplot(D_text, aes(x = 1, y = cluster)) +
  geom_text(aes(label = paste0(superclass, " | ", lineage, " | ", keep)), hjust = 0, size = 3.2) +
  theme_void(base_size = 11) +
  xlim(1, 8) +
  labs(title = "Cluster class | dominant background | decision") +
  theme(plot.title = element_text(hjust = 0, face = "bold"))

pD <- pD_heat + pD_text + plot_layout(widths = c(1.2, 1.6))

ggsave(file.path(cfg$outdir, "panels", "Figure1D_audit_summary.pdf"), pD, width = 11.5, height = 6.5)
ggsave(file.path(cfg$outdir, "panels", "Figure1D_audit_summary.png"), pD, width = 11.5, height = 6.5, dpi = 300)

fwrite(final_anno, file.path(cfg$outdir, "tables", "Figure1D_full_audit_table.tsv"), sep = "\t")

# ------------------------------------------------------------
# Panel E: strict retained core UMAP
# ------------------------------------------------------------
msg("Building Panel E ...")
strict_obj$cluster_use <- factor(as.character(strict_obj$cluster_use), levels = cfg$retained_clusters)
pE <- DimPlot(
  strict_obj,
  reduction = if ("umap" %in% names(strict_obj@reductions)) "umap" else names(strict_obj@reductions)[1],
  group.by = "cluster_use",
  label = TRUE,
  repel = TRUE,
  raster = FALSE,
  pt.size = cfg$point_size_umap
) +
  scale_color_manual(values = c("0" = "#3B6FB6", "2" = "#D95F02"), drop = FALSE) +
  theme_classic(base_size = 12) +
  labs(title = "Strict retained core microglial compartment") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

ggsave(file.path(cfg$outdir, "panels", "Figure1E_retained_core_umap.pdf"), pE, width = 6.2, height = 5.6)
ggsave(file.path(cfg$outdir, "panels", "Figure1E_retained_core_umap.png"), pE, width = 6.2, height = 5.6, dpi = 300)

# ------------------------------------------------------------
# Panel F: donor representation heatmap for retained 0/2
# ------------------------------------------------------------
msg("Building Panel F ...")
mdF <- strict_obj@meta.data
mdF$cell_id <- colnames(strict_obj)
mdF$cluster <- as.character(mdF$cluster_use)
mdF$donor <- as.character(mdF[[choose_donor_col(mdF)]])
mdF$dx <- normalize_dx(mdF[[choose_dx_col(mdF)]])

propF <- as.data.table(mdF)[, .N, by = .(donor, dx, cluster)]
propF_total <- as.data.table(mdF)[, .(n_total = .N), by = .(donor, dx)]
propF <- merge(propF, propF_total, by = c("donor", "dx"), all.x = TRUE)
propF[, prop := N / n_total]

# donor ordering by diagnosis then cluster2 proportion
wide2 <- dcast(propF, donor + dx ~ cluster, value.var = "prop", fill = 0)
if (!cfg$cand_cluster %in% colnames(wide2)) wide2[, (cfg$cand_cluster) := 0]
setorder(wide2, dx, -get(cfg$cand_cluster), donor)
donor_levels <- wide2$donor

propF[, donor := factor(donor, levels = donor_levels)]
propF[, cluster := factor(cluster, levels = cfg$retained_clusters)]

pF <- ggplot(propF, aes(x = cluster, y = donor, fill = prop)) +
  geom_tile(color = "white") +
  facet_grid(dx ~ ., scales = "free_y", space = "free_y") +
  scale_fill_gradient(low = "#f7fbff", high = "#08306b", limits = cfg$donor_heatmap_fill_limits, oob = squish) +
  theme_minimal(base_size = 11) +
  labs(title = "Donor-level representation of retained clusters 0 and 2", x = "Retained cluster", y = "Donor") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold")
  )

ggsave(file.path(cfg$outdir, "panels", "Figure1F_donor_representation_heatmap.pdf"), pF, width = 7.8, height = 8.5)
ggsave(file.path(cfg$outdir, "panels", "Figure1F_donor_representation_heatmap.png"), pF, width = 7.8, height = 8.5, dpi = 300)

fwrite(propF, file.path(cfg$outdir, "tables", "Figure1F_donor_representation.tsv"), sep = "\t")

# ------------------------------------------------------------
# Optional combined figure assembly (A optional)
# ------------------------------------------------------------
msg("Assembling combined figure layouts ...")
# use patchwork objects directly for reproducibility
B <- pB
C <- pC
D <- pD
E <- pE
F <- pF

# B-F layout
fig_BF <- (
  (B | E) /
  C /
  D /
  F
) + plot_annotation(tag_levels = "A")

ggsave(file.path(cfg$outdir, "Figure1_B_to_F_combined.pdf"), fig_BF, width = cfg$width_final, height = cfg$height_final)
ggsave(file.path(cfg$outdir, "Figure1_B_to_F_combined.png"), fig_BF, width = cfg$width_final, height = cfg$height_final, dpi = 300)

# save notes for manual insertion of panel A
notes <- c(
  "Figure 1 build completed.",
  "Panel A is intentionally not code-generated.",
  "Use Figure1A_nano_banana_prompt.txt to generate the workflow schematic.",
  "Then combine Panel A externally with Figure1_B_to_F_combined.pdf/png, or use a final layout in Illustrator/PowerPoint/Affinity/patchwork after exporting Panel A to PNG."
)
writeLines(notes, file.path(cfg$outdir, "Figure1_build_notes.txt"))

# provenance summary
summary_dt <- data.table(
  panel = c("A","B","C","D","E","F"),
  description = c(
    "Workflow schematic (external generation using Nano Banana prompt)",
    "All-cluster UMAP",
    "Transparent lineage/background marker dotplot",
    "Audit summary and retention decision",
    "Strict retained-core UMAP (clusters 0 and 2)",
    "Donor-level representation heatmap for retained clusters"
  ),
  output = c(
    "Figure1A_nano_banana_prompt.txt",
    "panels/Figure1B_allcluster_umap.pdf",
    "panels/Figure1C_lineage_marker_dotplot.pdf",
    "panels/Figure1D_audit_summary.pdf",
    "panels/Figure1E_retained_core_umap.pdf",
    "panels/Figure1F_donor_representation_heatmap.pdf"
  )
)
fwrite(summary_dt, file.path(cfg$outdir, "Figure1_manifest.tsv"), sep = "\t")

writeLines(capture.output(sessionInfo()), file.path(cfg$outdir, "sessionInfo_Figure1_final.txt"))
msg("Figure 1 final build finished.")

