#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(optparse)
})

options(stringsAsFactors = FALSE)


project_dir_default <- Sys.getenv("MA_PROJECT_DIR", ".")
option_list <- list(
  make_option("--microglia_rds", type = "character", default = Sys.getenv("MA_MICROGLIA_RDS", file.path(project_dir_default, "inputs", "Microglia_Clustered.rds")), help = "Discovery microglia Seurat object [default: %default]"),
  make_option("--package1b_audit", type = "character", default = file.path(project_dir_default, "results", "Package1b_ClusterAnnotation_ContaminantAudit", "tables", "10_cluster_annotation_contaminant_audit.tsv"), help = "Package1b audit table [default: %default]"),
  make_option("--package1b_strict_clusters", type = "character", default = file.path(project_dir_default, "results", "Package1b_ClusterAnnotation_ContaminantAudit", "tables", "11_strict_microglia_clusters.txt"), help = "Strict cluster list [default: %default]"),
  make_option("--strict_rds", type = "character", default = Sys.getenv("MA_STRICT_RDS", file.path(project_dir_default, "results", "Package2_strict_input_fixed_from_original.rds")), help = "Strict retained-core Seurat object [default: %default]"),
  make_option("--assay", type = "character", default = "RNA", help = "Assay name [default: %default]"),
  make_option("--outdir", type = "character", default = file.path(project_dir_default, "results", "Figure1_Final_NoWorkflow"), help = "Output directory [default: %default]"),
  make_option("--ref_cluster", type = "character", default = "0", help = "Reference cluster ID [default: %default]"),
  make_option("--cand_cluster", type = "character", default = "2", help = "Candidate cluster ID [default: %default]"),
  make_option("--point_size_umap", type = "double", default = 0.45, help = "UMAP point size [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

cfg <- list(
  microglia_rds = opt$microglia_rds,
  package1b_audit = opt$package1b_audit,
  package1b_strict_clusters = opt$package1b_strict_clusters,
  strict_rds = opt$strict_rds,
  assay = opt$assay,
  outdir = opt$outdir,
  retained_clusters = c(opt$ref_cluster, opt$cand_cluster),
  ref_cluster = opt$ref_cluster,
  cand_cluster = opt$cand_cluster,
  point_size_umap = opt$point_size_umap,
  donor_fill_limits = c(0, 1)
)

msg <- function(...) {
  cat(sprintf("%s | %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = "")))
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
  "cluster_auto_from_idents"
}

choose_donor_col <- function(meta) {
  find_first_col(meta, c("^donor$", "donor_id", "individual", "subject", "patient", "sample", "orig.ident", "library", "case_id"))
}

choose_dx_col <- function(meta) {
  find_first_col(meta, c("^dx$", "^group$", "diagnosis", "condition", "status", "phenotype", "disease"))
}

normalize_dx <- function(x) {
  x <- trimws(as.character(x))
  y <- x
  y[grepl("asd|autism|case", x, ignore.case = TRUE)] <- "ASD"
  y[grepl("control|ctrl|ctl|healthy|non[- ]?psy", x, ignore.case = TRUE)] <- "Control"
  y
}

rescale01_safe <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  rng <- range(x, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0) return(rep(0.5, length(x)))
  (x - rng[1]) / diff(rng)
}

# transparent lineage/background marker panel
lineage_marker_panel <- c(
  "P2RY12","CX3CR1","TMEM119","MERTK","C3",
  "SYT1","NRXN3","RBFOX1","GPM6A","GRM5",
  "GFAP","AQP4","SLC1A3",
  "MBP","MOG","PLP1",
  "PDGFRA","CSPG4",
  "CLDN5","FLT1","PECAM1","NKG7","CD3D"
)

# output dirs

dir.create(cfg$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cfg$outdir, "panels"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cfg$outdir, "tables"), recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# read objects/tables
# ------------------------------------------------------------
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
obj$donor_use <- as.character(obj@meta.data[[donor_col]])
obj$dx_use <- normalize_dx(obj@meta.data[[dx_col]])
Idents(obj) <- obj$cluster_use

cluster_levels <- unique(as.character(obj$cluster_use))
num_try <- suppressWarnings(as.numeric(cluster_levels))
cluster_levels <- if (all(!is.na(num_try))) as.character(sort(num_try)) else sort(cluster_levels)
obj$cluster_use <- factor(as.character(obj$cluster_use), levels = cluster_levels)
Idents(obj) <- obj$cluster_use

final_anno <- fread(cfg$package1b_audit)
final_anno[, cluster := as.character(cluster)]

strict_clusters <- if (file.exists(cfg$package1b_strict_clusters)) unique(trimws(readLines(cfg$package1b_strict_clusters))) else cfg$retained_clusters
if (!all(cfg$retained_clusters %in% strict_clusters)) strict_clusters <- cfg$retained_clusters

strict_obj <- if (file.exists(cfg$strict_rds)) readRDS(cfg$strict_rds) else subset(obj, cells = colnames(obj)[as.character(obj$cluster_use) %in% strict_clusters])
if (cfg$assay %in% Assays(strict_obj)) DefaultAssay(strict_obj) <- cfg$assay
if (!"cluster_use" %in% colnames(strict_obj@meta.data)) {
  strict_obj$cluster_use <- as.character(strict_obj@meta.data[[choose_cluster_col(strict_obj)]])
}
if (!"donor_use" %in% colnames(strict_obj@meta.data)) {
  sc_donor_col <- choose_donor_col(strict_obj@meta.data)
  strict_obj$donor_use <- as.character(strict_obj@meta.data[[sc_donor_col]])
}
if (!"dx_use" %in% colnames(strict_obj@meta.data)) {
  sc_dx_col <- choose_dx_col(strict_obj@meta.data)
  strict_obj$dx_use <- normalize_dx(strict_obj@meta.data[[sc_dx_col]])
}
strict_obj$cluster_use <- factor(as.character(strict_obj$cluster_use), levels = cfg$retained_clusters)
Idents(strict_obj) <- strict_obj$cluster_use

# save metadata used
fwrite(
  data.table(
    item = c("microglia_rds", "audit_table", "strict_rds", "cluster_col", "donor_col", "dx_col"),
    value = c(cfg$microglia_rds, cfg$package1b_audit, cfg$strict_rds, cluster_col, donor_col, dx_col)
  ),
  file.path(cfg$outdir, "tables", "00_inputs_and_columns.tsv"),
  sep = "\t"
)

# ------------------------------------------------------------
# Panel A: all-cluster UMAP
# ------------------------------------------------------------
msg("Building Panel A ...")
other_clusters <- setdiff(cluster_levels, cfg$retained_clusters)
other_cols <- if (length(other_clusters) > 0) setNames(rep("#BFBFBF", length(other_clusters)), other_clusters) else c()
cluster_cols <- c(
  other_cols,
  setNames("#5B8CC0", cfg$ref_cluster),
  setNames("#E67E4E", cfg$cand_cluster)
)
cluster_cols <- cluster_cols[cluster_levels]

pA <- suppressWarnings(
  DimPlot(
    obj,
    reduction = if ("umap" %in% names(obj@reductions)) "umap" else names(obj@reductions)[1],
    group.by = "cluster_use",
    label = TRUE,
    repel = TRUE,
    raster = FALSE,
    pt.size = cfg$point_size_umap
  )
) +
  scale_color_manual(values = cluster_cols, drop = FALSE) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold")
  ) +
  labs(title = "A  Reclusted discovery microglial object")

ggsave(file.path(cfg$outdir, "panels", "Figure1A_allcluster_umap.pdf"), pA, width = 6.8, height = 5.7)
ggsave(file.path(cfg$outdir, "panels", "Figure1A_allcluster_umap.png"), pA, width = 6.8, height = 5.7, dpi = 300)

# ------------------------------------------------------------
# Panel B: lineage/background marker dotplot
# ------------------------------------------------------------
msg("Building Panel B ...")
features_B <- intersect(lineage_marker_panel, rownames(obj))
if (length(features_B) == 0) stop("No lineage/background marker genes were found in the object.")

pB <- suppressWarnings(
  DotPlot(obj, features = features_B, group.by = "cluster_use")
) +
  RotatedAxis() +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1),
    plot.title = element_text(face = "bold")
  ) +
  labs(title = "B  Transparent lineage/background marker panel", x = NULL, y = "Cluster")

ggsave(file.path(cfg$outdir, "panels", "Figure1B_lineage_dotplot.pdf"), pB, width = 13.8, height = 4.8)
ggsave(file.path(cfg$outdir, "panels", "Figure1B_lineage_dotplot.png"), pB, width = 13.8, height = 4.8, dpi = 300)

# ------------------------------------------------------------
# Panel C: audit summary + retention decision
# ------------------------------------------------------------
msg("Building Panel C ...")
score_cols <- grep("_total_score$", names(final_anno), value = TRUE)
mg_col <- if ("microglia_total_score" %in% score_cols) "microglia_total_score" else NA_character_
nonmg_cols <- setdiff(score_cols, mg_col)

final_anno[, top_nonmg_score := if (length(nonmg_cols) > 0) do.call(pmax, c(.SD, na.rm = TRUE)) else NA_real_, .SDcols = nonmg_cols]
final_anno[, donor_representation := fifelse("n_donors" %in% names(final_anno), as.numeric(n_donors), NA_real_)]
final_anno[, retained_numeric := fifelse("keep_strict" %in% names(final_anno) & keep_strict, 1, 0)]

plotC <- data.table(
  cluster = rep(final_anno$cluster, each = 4),
  metric = rep(c("Microglial\nidentity", "Non-microglial\nbackground", "Donor\nrepresentation", "Retained"), times = nrow(final_anno)),
  value = c(
    rescale01_safe(if (!is.na(mg_col)) final_anno[[mg_col]] else NA_real_),
    rescale01_safe(final_anno$top_nonmg_score),
    rescale01_safe(final_anno$donor_representation),
    final_anno$retained_numeric
  )
)
plotC[, cluster := factor(cluster, levels = rev(cluster_levels))]
plotC[, metric := factor(metric, levels = c("Microglial\nidentity", "Non-microglial\nbackground", "Donor\nrepresentation", "Retained"))]
plotC[, label := fifelse(metric == "Retained", ifelse(value >= 0.5, "Yes", "No"), sprintf("%.2f", value))]

pC <- ggplot(plotC, aes(x = metric, y = cluster, fill = value)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = label), size = 3) +
  scale_fill_gradient(low = "#F3F3F3", high = "#4C78A8", limits = c(0, 1), na.value = "#F3F3F3") +
  theme_classic(base_size = 11) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(face = "bold"),
    legend.position = "none",
    plot.title = element_text(face = "bold")
  ) +
  labs(title = "C  Audit summary and retention decision")

ggsave(file.path(cfg$outdir, "panels", "Figure1C_audit_summary.pdf"), pC, width = 6.6, height = 5.7)
ggsave(file.path(cfg$outdir, "panels", "Figure1C_audit_summary.png"), pC, width = 6.6, height = 5.7, dpi = 300)

# ------------------------------------------------------------
# Panel D: strict retained-core UMAP
# ------------------------------------------------------------
msg("Building Panel D ...")
ret_cols <- c(setNames("#5B8CC0", cfg$ref_cluster), setNames("#E67E4E", cfg$cand_cluster))
ret_cols <- ret_cols[cfg$retained_clusters]

pD <- suppressWarnings(
  DimPlot(
    strict_obj,
    reduction = if ("umap" %in% names(strict_obj@reductions)) "umap" else names(strict_obj@reductions)[1],
    group.by = "cluster_use",
    label = TRUE,
    repel = TRUE,
    raster = FALSE,
    pt.size = cfg$point_size_umap
  )
) +
  scale_color_manual(values = ret_cols, drop = FALSE) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold")
  ) +
  labs(title = "D  Strict retained core microglia")

ggsave(file.path(cfg$outdir, "panels", "Figure1D_retainedcore_umap.pdf"), pD, width = 5.8, height = 5.2)
ggsave(file.path(cfg$outdir, "panels", "Figure1D_retainedcore_umap.png"), pD, width = 5.8, height = 5.2, dpi = 300)

# ------------------------------------------------------------
# Panel E: donor representation heatmap for retained clusters
# ------------------------------------------------------------
msg("Building Panel E ...")
md <- as.data.table(strict_obj@meta.data, keep.rownames = "cell")
md[, donor := as.character(donor_use)]
md[, dx := normalize_dx(dx_use)]
md[, cluster := as.character(cluster_use)]
md <- md[!is.na(donor) & donor != "" & !is.na(cluster)]

prop_dt <- md[, .N, by = .(donor, dx, cluster)]
prop_dt[, donor_total := sum(N), by = donor]
prop_dt[, prop_within_donor := N / donor_total]

wide_dt <- dcast(prop_dt, donor + dx ~ cluster, value.var = "prop_within_donor", fill = 0)
if (!(cfg$ref_cluster %in% names(wide_dt))) wide_dt[, (cfg$ref_cluster) := 0]
if (!(cfg$cand_cluster %in% names(wide_dt))) wide_dt[, (cfg$cand_cluster) := 0]
wide_dt[, ref_prop := get(cfg$ref_cluster)]
wide_dt[, cand_prop := get(cfg$cand_cluster)]
wide_dt[, dx_order := fifelse(dx == "ASD", 1L, 0L)]
setorder(wide_dt, dx_order, -cand_prop, -ref_prop, donor)
donor_order <- wide_dt$donor

plotE <- copy(prop_dt)
plotE[, donor := factor(donor, levels = donor_order)]
plotE[, cluster := factor(cluster, levels = cfg$retained_clusters)]
plotE[, dx := factor(dx, levels = c("Control", "ASD"))]

pE <- ggplot(plotE, aes(x = donor, y = cluster, fill = prop_within_donor)) +
  geom_tile(color = "white", linewidth = 0.15) +
  facet_grid(. ~ dx, scales = "free_x", space = "free_x") +
  scale_fill_gradient(low = "#F3F3F3", high = "#3B6FB6", limits = cfg$donor_fill_limits, oob = squish) +
  theme_classic(base_size = 11) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  ) +
  labs(title = "E  Donor-level representation of retained clusters", y = NULL, fill = "Within-donor\nproportion")

ggsave(file.path(cfg$outdir, "panels", "Figure1E_donor_representation.pdf"), pE, width = 12.6, height = 2.8)
ggsave(file.path(cfg$outdir, "panels", "Figure1E_donor_representation.png"), pE, width = 12.6, height = 2.8, dpi = 300)

# export donor table
fwrite(wide_dt[, .(donor, dx, ref_prop, cand_prop)], file.path(cfg$outdir, "tables", "01_retained_donor_representation.tsv"), sep = "\t")

# ------------------------------------------------------------
# combined figure
# ------------------------------------------------------------
msg("Combining panels ...")
combined <- (pA | pB) / (pC | pD | pE) + plot_layout(heights = c(1.15, 1))

ggsave(file.path(cfg$outdir, "Figure1_A_to_E_combined.pdf"), combined, width = 16.5, height = 10.8)
ggsave(file.path(cfg$outdir, "Figure1_A_to_E_combined.png"), combined, width = 16.5, height = 10.8, dpi = 300)

msg("Done. Output dir: ", cfg$outdir)
