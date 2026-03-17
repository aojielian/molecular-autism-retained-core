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
  vel_raw_input = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/Velmeshev_scRNA/rawMatrix.zip",
  vel_meta_file = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/Velmeshev_scRNA/meta.tsv",
  outdir = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Figure2_Enhanced_v2",
  assay = "RNA",
  cluster_col = NULL,
  donor_col = NULL,
  dx_col = NULL,
  barcode_col = NULL,
  celltype_col = NULL,
  risk_cluster = "2",
  homeo_cluster = "4",
  top_n = 50,
  min_pct = 0.10,
  logfc = 0.25
)

dir.create(cfg$outdir, recursive = TRUE, showWarnings = FALSE)

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

normalize_dx <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  y <- x
  y[grepl("asd|autism|case", x, ignore.case = TRUE)] <- "ASD"
  y[grepl("control|ctrl|ctl|healthy|non[- ]?psy", x, ignore.case = TRUE)] <- "Control"
  y
}

strip_bc <- function(x) {
  x <- as.character(x)
  x <- sub("-1$", "", x)
  x
}

choose_cluster_col <- function(obj, user_col = NULL) {
  meta <- obj@meta.data
  if (!is.null(user_col) && user_col %in% names(meta)) return(user_col)
  if ("seurat_clusters" %in% names(meta)) return("seurat_clusters")
  cand <- find_first_col(meta, c("^cluster$", "cluster", "seurat", "res\\.", "RNA_snn"))
  if (!is.null(cand)) return(cand)
  obj$cluster_auto_from_idents <- as.character(Idents(obj))
  return("cluster_auto_from_idents")
}

get_fc_col <- function(markers) {
  fc_candidates <- c("avg_log2FC", "avg_logFC", "log2FC", "logFC")
  fc_candidates[fc_candidates %in% colnames(markers)][1]
}

harmonize_genes <- function(genes, features) {
  genes <- unique(as.character(genes))
  genes <- genes[!is.na(genes) & nzchar(genes)]

  exact <- intersect(genes, features)
  if (length(exact) >= 5) return(exact)

  fmap <- setNames(features, toupper(features))
  mapped <- unname(fmap[toupper(genes)])
  mapped <- unique(mapped[!is.na(mapped)])
  mapped
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

read_counts_auto <- function(path) {
  if (grepl("\\.zip$", path, ignore.case = TRUE)) {
    exdir <- file.path(tempdir(), paste0("Velmeshev_raw_", format(Sys.time(), "%Y%m%d_%H%M%S")))
    dir.create(exdir, recursive = TRUE, showWarnings = FALSE)
    unzip(path, exdir = exdir)
    path <- exdir
  }

  if (dir.exists(path)) {
    mtx_files <- list.files(path, pattern = "matrix\\.mtx(\\.gz)?$", recursive = TRUE, full.names = TRUE)
    if (length(mtx_files) == 0) stop("No matrix.mtx(.gz) found.")
    candidate_dirs <- unique(dirname(mtx_files))

    pick_dir <- NULL
    for (d in candidate_dirs) {
      ff <- list.files(d)
      has_mtx <- any(grepl("matrix\\.mtx(\\.gz)?$", ff))
      has_feat <- any(grepl("(features|genes)\\.tsv(\\.gz)?$", ff))
      has_bc <- any(grepl("barcodes\\.tsv(\\.gz)?$", ff))
      if (has_mtx && has_feat && has_bc) {
        pick_dir <- d
        break
      }
    }
    if (is.null(pick_dir)) pick_dir <- candidate_dirs[1]

    msg("Read10X directory: ", pick_dir)
    x <- Read10X(data.dir = pick_dir)
    if (is.list(x)) {
      if ("Gene Expression" %in% names(x)) {
        x <- x[["Gene Expression"]]
      } else {
        x <- x[[1]]
      }
    }
    return(x)
  }

  stop("Unsupported raw_input format: ", path)
}

detect_barcode_col <- function(meta, cell_ids, user_col = NULL) {
  if (!is.null(user_col) && user_col %in% names(meta)) return(user_col)

  cand <- names(meta)[sapply(meta, function(z) is.character(z) || is.factor(z))]
  if (length(cand) == 0) return(NULL)

  exact_hits <- sapply(cand, function(cc) sum(as.character(meta[[cc]]) %in% cell_ids))
  if (max(exact_hits) > 0) return(names(which.max(exact_hits)))

  stripped_ids <- strip_bc(cell_ids)
  stripped_hits <- sapply(cand, function(cc) sum(strip_bc(as.character(meta[[cc]])) %in% stripped_ids))
  if (max(stripped_hits) > 0) return(names(which.max(stripped_hits)))

  find_first_col(meta, c("barcode", "cell", "barcodes"))
}

align_meta_to_cells <- function(meta, barcode_col, cell_ids) {
  bc <- as.character(meta[[barcode_col]])
  idx <- match(cell_ids, bc)

  match_rate <- mean(!is.na(idx))
  if (match_rate < 0.5) {
    idx <- match(strip_bc(cell_ids), strip_bc(bc))
  }

  meta2 <- meta[idx, , drop = FALSE]
  rownames(meta2) <- cell_ids
  meta2
}

detect_celltype_col <- function(meta, user_col = NULL) {
  if (!is.null(user_col) && user_col %in% names(meta)) return(user_col)

  char_cols <- names(meta)[sapply(meta, function(z) is.character(z) || is.factor(z))]
  if (length(char_cols) == 0) return(NULL)

  priority <- char_cols[grepl("cell.?type|annotation|class|subclass|type|label|cluster", char_cols, ignore.case = TRUE)]
  cand <- unique(c(priority, char_cols))

  score <- sapply(cand, function(cc) {
    v <- as.character(meta[[cc]])
    sum(grepl("microgl", v, ignore.case = TRUE), na.rm = TRUE) +
      sum(grepl("^MG$", v, ignore.case = TRUE), na.rm = TRUE)
  })

  if (max(score) == 0) return(NULL)
  names(which.max(score))
}

detect_donor_col <- function(meta, user_col = NULL) {
  if (!is.null(user_col) && user_col %in% names(meta)) return(user_col)
  find_first_col(meta, c("^donor$", "donor_id", "individual", "subject", "patient", "sample", "orig.ident", "library", "case_id"))
}

detect_dx_col <- function(meta, user_col = NULL) {
  if (!is.null(user_col) && user_col %in% names(meta)) return(user_col)
  find_first_col(meta, c("^dx$", "^group$", "diagnosis", "condition", "status", "phenotype", "disease"))
}

# -----------------------------
# 1. 从 discovery microglia 对象提取 cluster marker signatures
# -----------------------------
msg("Reading discovery microglia object ...")
obj <- readRDS(cfg$microglia_rds)
DefaultAssay(obj) <- if (cfg$assay %in% Assays(obj)) cfg$assay else DefaultAssay(obj)

cluster_col <- choose_cluster_col(obj, cfg$cluster_col)
obj$cluster_use <- as.character(obj@meta.data[[cluster_col]])
Idents(obj) <- obj$cluster_use

msg("Running FindAllMarkers on discovery microglia ...")
markers <- FindAllMarkers(
  object = obj,
  only.pos = TRUE,
  min.pct = cfg$min_pct,
  logfc.threshold = cfg$logfc
)

fwrite(markers, file.path(cfg$outdir, "Discovery_Microglia_AllMarkers.csv"))

fc_col <- get_fc_col(markers)
if (is.na(fc_col) || is.null(fc_col)) stop("Could not find fold-change column in markers output.")

risk_markers <- markers %>%
  filter(cluster == cfg$risk_cluster) %>%
  arrange(desc(.data[[fc_col]]), p_val_adj)

homeo_markers <- markers %>%
  filter(cluster == cfg$homeo_cluster) %>%
  arrange(desc(.data[[fc_col]]), p_val_adj)

if (nrow(risk_markers) == 0) stop("No markers found for risk cluster: ", cfg$risk_cluster)
if (nrow(homeo_markers) == 0) stop("No markers found for homeo cluster: ", cfg$homeo_cluster)

risk_genes <- head(unique(risk_markers$gene), cfg$top_n)
homeo_genes <- head(unique(homeo_markers$gene), cfg$top_n)

writeLines(risk_genes, file.path(cfg$outdir, "RiskSignature_Cluster2_top50.txt"))
writeLines(homeo_genes, file.path(cfg$outdir, "HomeoSignature_Cluster4_top50.txt"))

sig_audit <- data.frame(
  signature = c("Risk", "Homeostatic"),
  cluster = c(cfg$risk_cluster, cfg$homeo_cluster),
  n_genes = c(length(risk_genes), length(homeo_genes)),
  preview = c(paste(head(risk_genes, 12), collapse = ", "),
              paste(head(homeo_genes, 12), collapse = ", "))
)
fwrite(sig_audit, file.path(cfg$outdir, "ClusterMarker_Signature_Audit.tsv"), sep = "\t")

# -----------------------------
# 2. 读取 Velmeshev counts + metadata
# -----------------------------
msg("Reading Velmeshev counts ...")
counts <- read_counts_auto(cfg$vel_raw_input)

msg("Reading Velmeshev metadata ...")
meta <- fread(cfg$vel_meta_file, data.table = FALSE)
meta <- as.data.frame(meta)

barcode_col <- detect_barcode_col(meta, colnames(counts), cfg$barcode_col)
if (is.null(barcode_col)) stop("Could not identify barcode column in metadata.")

meta2 <- align_meta_to_cells(meta, barcode_col, colnames(counts))
keep <- !is.na(meta2[[barcode_col]])
meta2 <- meta2[keep, , drop = FALSE]
counts <- counts[, rownames(meta2), drop = FALSE]

celltype_col <- detect_celltype_col(meta2, cfg$celltype_col)
if (is.null(celltype_col)) stop("Could not identify cell-type column containing microglia labels.")

celltype_values <- as.character(meta2[[celltype_col]])
mg_idx <- grepl("microgl", celltype_values, ignore.case = TRUE) | grepl("^MG$", celltype_values, ignore.case = TRUE)

if (sum(mg_idx) == 0) stop("No microglia cells matched in column: ", celltype_col)

meta_mg <- meta2[mg_idx, , drop = FALSE]
counts_mg <- counts[, rownames(meta_mg), drop = FALSE]

donor_col <- detect_donor_col(meta_mg, cfg$donor_col)
dx_col <- detect_dx_col(meta_mg, cfg$dx_col)

if (is.null(donor_col)) stop("Could not identify donor column.")
if (is.null(dx_col)) stop("Could not identify diagnosis column.")

audit_meta <- data.frame(
  item = c("barcode_col", "celltype_col", "donor_col", "dx_col", "n_all_cells", "n_mg_cells"),
  value = c(barcode_col, celltype_col, donor_col, dx_col, ncol(counts), ncol(counts_mg))
)
fwrite(audit_meta, file.path(cfg$outdir, "Figure2_Metadata_Audit_v2.tsv"), sep = "\t")

# -----------------------------
# 3. Velmeshev MG 打分
# -----------------------------
msg("Building Velmeshev microglia object ...")
mg <- CreateSeuratObject(
  counts = counts_mg,
  meta.data = meta_mg,
  assay = cfg$assay,
  project = "VelmeshevMG",
  min.cells = 0,
  min.features = 0
)
DefaultAssay(mg) <- cfg$assay
mg <- NormalizeData(mg, verbose = FALSE)

features_use <- rownames(mg)
risk_genes_use <- harmonize_genes(risk_genes, features_use)
homeo_genes_use <- harmonize_genes(homeo_genes, features_use)

mapped_tbl <- data.frame(
  signature = c("Risk", "Homeostatic"),
  n_input = c(length(risk_genes), length(homeo_genes)),
  n_mapped = c(length(risk_genes_use), length(homeo_genes_use)),
  preview = c(paste(head(risk_genes_use, 12), collapse = ", "),
              paste(head(homeo_genes_use, 12), collapse = ", "))
)
fwrite(mapped_tbl, file.path(cfg$outdir, "Mapped_Signature_Audit_v2.tsv"), sep = "\t")

mg <- AddModuleScore(mg, features = list(risk_genes_use), name = "RiskScore", assay = cfg$assay)
mg <- AddModuleScore(mg, features = list(homeo_genes_use), name = "HomeoScore", assay = cfg$assay)

md <- mg@meta.data
md$cell_id <- colnames(mg)
md$donor <- as.character(md[[donor_col]])
md$dx <- normalize_dx(md[[dx_col]])
md$Risk <- md$RiskScore1
md$Homeo <- md$HomeoScore1

# -----------------------------
# 4. donor-level 汇总
# -----------------------------
msg("Aggregating donor-level scores ...")
donor_df <- md %>%
  group_by(donor, dx) %>%
  summarise(
    n_cells = n(),
    Risk_mean = mean(Risk, na.rm = TRUE),
    Risk_median = median(Risk, na.rm = TRUE),
    Homeo_mean = mean(Homeo, na.rm = TRUE),
    Homeo_median = median(Homeo, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(dx, donor)

fwrite(donor_df, file.path(cfg$outdir, "Validation_DonorLevel_Scores_v2.csv"))
fwrite(
  md[, c("cell_id", "donor", "dx", "Risk", "Homeo"), drop = FALSE],
  file.path(cfg$outdir, "Validation_CellLevel_Scores_v2.tsv"),
  sep = "\t"
)

# -----------------------------
# 5. 出图
# -----------------------------
order_dx <- c("Control", "ASD")
donor_df$dx <- factor(donor_df$dx, levels = unique(c(order_dx, as.character(donor_df$dx))))

p_risk <- plot_p_from_wilcox(donor_df, "Risk_mean", "dx")
p_homeo <- plot_p_from_wilcox(donor_df, "Homeo_mean", "dx")

p_lab_risk <- ifelse(is.na(p_risk), "Wilcoxon p = NA", paste0("Wilcoxon p = ", signif(p_risk, 3)))
p_lab_homeo <- ifelse(is.na(p_homeo), "Wilcoxon p = NA", paste0("Wilcoxon p = ", signif(p_homeo, 3)))

g_risk <- ggplot(donor_df, aes(x = dx, y = Risk_mean, fill = dx)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.75) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.9) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 20, hjust = 1)
  ) +
  labs(
    title = "Velmeshev microglia donor-level risk score",
    subtitle = p_lab_risk,
    x = NULL, y = "Mean module score per donor"
  )

g_homeo <- ggplot(donor_df, aes(x = dx, y = Homeo_mean, fill = dx)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.75) +
  geom_jitter(width = 0.12, size = 2, alpha = 0.9) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 20, hjust = 1)
  ) +
  labs(
    title = "Velmeshev microglia donor-level homeostatic score",
    subtitle = p_lab_homeo,
    x = NULL, y = "Mean module score per donor"
  )

ggsave(file.path(cfg$outdir, "Fig2B_Risk_DonorLevel_v2.pdf"), g_risk, width = 5.5, height = 4.8)
ggsave(file.path(cfg$outdir, "Fig2C_Homeo_DonorLevel_v2.pdf"), g_homeo, width = 5.5, height = 4.8)

msg("Figure 2 v2 finished.")
