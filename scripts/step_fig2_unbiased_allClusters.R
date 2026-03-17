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
  # discovery
  microglia_rds = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/PsychENCODE_Science_2024/Step3_Microglia/Microglia_Clustered.rds",

  # validation
  vel_raw_input = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/Velmeshev_scRNA/rawMatrix.zip",
  vel_meta_file = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/Velmeshev_scRNA/meta.tsv",

  # output
  outdir = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Figure2_Unbiased_AllClusters",

  # general
  assay = "RNA",
  cluster_col = NULL,
  donor_col_discovery = NULL,
  dx_col_discovery = NULL,
  donor_col_validation = NULL,
  dx_col_validation = NULL,
  barcode_col_validation = NULL,
  celltype_col_validation = NULL,

  # marker/signature parameters
  top_n = 50,
  min_pct = 0.10,
  logfc = 0.25,

  # ranking preference
  prefer_positive_direction = TRUE,
  min_genes_mapped = 5
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

strip_bc <- function(x) {
  x <- as.character(x)
  x <- sub("-1$", "", x)
  x
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

read_counts_auto <- function(path) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    x <- readRDS(path)
    if (inherits(x, "Seurat")) return(GetAssayData(x, assay = DefaultAssay(x), slot = "counts"))
    if (inherits(x, "dgCMatrix")) return(x)
    stop("RDS loaded but object is neither Seurat nor dgCMatrix.")
  }

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
      has_mtx  <- any(grepl("matrix\\.mtx(\\.gz)?$", ff))
      has_feat <- any(grepl("(features|genes)\\.tsv(\\.gz)?$", ff))
      has_bc   <- any(grepl("barcodes\\.tsv(\\.gz)?$", ff))
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

  mean(val[grp == "ASD"], na.rm = TRUE) -
    mean(val[grp == "Control"], na.rm = TRUE)
}

make_signature_name <- function(cl) paste0("Cluster", cl)

# --------------------------------------------------
# 1. discovery: read microglia object
# --------------------------------------------------
msg("Reading discovery microglia object ...")
obj <- readRDS(cfg$microglia_rds)
if (!inherits(obj, "Seurat")) stop("Discovery RDS is not a Seurat object.")

if (cfg$assay %in% Assays(obj)) {
  DefaultAssay(obj) <- cfg$assay
}

cluster_col <- choose_cluster_col(obj, cfg$cluster_col)
obj$cluster_use <- as.character(obj@meta.data[[cluster_col]])
Idents(obj) <- obj$cluster_use

meta_disc <- obj@meta.data
donor_col_disc <- choose_donor_col(meta_disc, cfg$donor_col_discovery)
dx_col_disc    <- choose_dx_col(meta_disc, cfg$dx_col_discovery)

if (is.null(donor_col_disc)) stop("Could not identify donor column in discovery object.")
if (is.null(dx_col_disc)) stop("Could not identify diagnosis column in discovery object.")

audit_disc <- data.frame(
  item = c("n_cells", "n_features", "default_assay", "cluster_col", "donor_col_discovery", "dx_col_discovery"),
  value = c(ncol(obj), nrow(obj), DefaultAssay(obj), cluster_col, donor_col_disc, dx_col_disc)
)
fwrite(audit_disc, file.path(cfg$outdir, "Discovery_Audit.tsv"), sep = "\t")

# --------------------------------------------------
# 2. discovery: find markers for all clusters
# --------------------------------------------------
msg("Running FindAllMarkers on discovery microglia ...")
markers <- FindAllMarkers(
  object = obj,
  only.pos = TRUE,
  min.pct = cfg$min_pct,
  logfc.threshold = cfg$logfc
)
fwrite(markers, file.path(cfg$outdir, "Discovery_AllCluster_Markers.csv"))

fc_col <- get_fc_col(markers)
if (is.null(fc_col)) stop("Could not identify fold-change column in marker table.")

cluster_levels <- sort(unique(as.character(obj$cluster_use)))

top10_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = .data[[fc_col]], n = 10, with_ties = FALSE) %>%
  ungroup()
fwrite(top10_markers, file.path(cfg$outdir, "Discovery_AllCluster_Markers_Top10.csv"))

# generate signatures
dir.create(file.path(cfg$outdir, "Cluster_Signatures"), showWarnings = FALSE)
signature_list <- list()
signature_audit <- list()

for (cl in cluster_levels) {
  subm <- markers %>%
    filter(cluster == cl) %>%
    arrange(desc(.data[[fc_col]]), p_val_adj)

  genes <- head(unique(subm$gene), cfg$top_n)
  signature_list[[make_signature_name(cl)]] <- genes

  writeLines(
    genes,
    file.path(cfg$outdir, "Cluster_Signatures", paste0(make_signature_name(cl), "_top", cfg$top_n, ".txt"))
  )

  signature_audit[[length(signature_audit) + 1]] <- data.frame(
    cluster = cl,
    signature_name = make_signature_name(cl),
    n_genes = length(genes),
    top_marker_preview = paste(head(genes, 12), collapse = ", "),
    stringsAsFactors = FALSE
  )
}
signature_audit_df <- rbindlist(signature_audit, fill = TRUE)
fwrite(signature_audit_df, file.path(cfg$outdir, "Discovery_ClusterSignature_Audit.tsv"), sep = "\t")

# --------------------------------------------------
# 3. discovery: donor-level abundance for all clusters
# --------------------------------------------------
msg("Computing discovery donor-level cluster abundance ...")
meta_disc2 <- obj@meta.data
meta_disc2$cell_id <- colnames(obj)
meta_disc2$donor <- as.character(meta_disc2[[donor_col_disc]])
meta_disc2$dx <- normalize_dx(meta_disc2[[dx_col_disc]])
meta_disc2$cluster <- as.character(meta_disc2$cluster_use)

disc_prop_df <- meta_disc2 %>%
  count(donor, dx, cluster, name = "n_cluster") %>%
  left_join(
    meta_disc2 %>% count(donor, dx, name = "n_total"),
    by = c("donor", "dx")
  ) %>%
  mutate(prop = n_cluster / n_total) %>%
  arrange(cluster, dx, donor)

fwrite(disc_prop_df, file.path(cfg$outdir, "Discovery_AllCluster_Proportions.csv"))

disc_stats <- lapply(cluster_levels, function(cl) {
  sub <- disc_prop_df %>% filter(cluster == cl)
  data.frame(
    cluster = cl,
    signature_name = make_signature_name(cl),
    n_donors = length(unique(sub$donor)),
    mean_prop_ASD = mean(sub$prop[sub$dx == "ASD"], na.rm = TRUE),
    mean_prop_Control = mean(sub$prop[sub$dx == "Control"], na.rm = TRUE),
    effect_ASD_minus_Control = mean_diff_asd_minus_ctrl(sub, "prop", "dx"),
    wilcox_p = calc_wilcox_p(sub, "prop", "dx"),
    stringsAsFactors = FALSE
  )
})
disc_stats_df <- rbindlist(disc_stats, fill = TRUE)
disc_stats_df$wilcox_fdr <- p.adjust(disc_stats_df$wilcox_p, method = "BH")
fwrite(disc_stats_df, file.path(cfg$outdir, "Discovery_AllCluster_AbundanceStats.tsv"), sep = "\t")

# --------------------------------------------------
# 4. validation: read Velmeshev counts + metadata
# --------------------------------------------------
msg("Reading validation counts ...")
counts <- read_counts_auto(cfg$vel_raw_input)
msg("Validation counts dim: ", paste(dim(counts), collapse = " x "))

msg("Reading validation metadata ...")
meta_val <- fread(cfg$vel_meta_file, data.table = FALSE)
meta_val <- as.data.frame(meta_val)

barcode_col_val <- detect_barcode_col(meta_val, colnames(counts), cfg$barcode_col_validation)
if (is.null(barcode_col_val)) stop("Could not identify barcode column in validation metadata.")

meta_val2 <- align_meta_to_cells(meta_val, barcode_col_val, colnames(counts))
keep <- !is.na(meta_val2[[barcode_col_val]])
meta_val2 <- meta_val2[keep, , drop = FALSE]
counts <- counts[, rownames(meta_val2), drop = FALSE]

celltype_col_val <- detect_celltype_col(meta_val2, cfg$celltype_col_validation)
if (is.null(celltype_col_val)) stop("Could not identify validation cell-type column containing microglia labels.")

ctv <- as.character(meta_val2[[celltype_col_val]])
mg_idx <- grepl("microgl", ctv, ignore.case = TRUE) | grepl("^MG$", ctv, ignore.case = TRUE)
if (sum(mg_idx) == 0) stop("No microglia cells matched in validation metadata.")

meta_val_mg <- meta_val2[mg_idx, , drop = FALSE]
counts_val_mg <- counts[, rownames(meta_val_mg), drop = FALSE]

donor_col_val <- choose_donor_col(meta_val_mg, cfg$donor_col_validation)
dx_col_val    <- choose_dx_col(meta_val_mg, cfg$dx_col_validation)

if (is.null(donor_col_val)) stop("Could not identify donor column in validation metadata.")
if (is.null(dx_col_val)) stop("Could not identify diagnosis column in validation metadata.")

audit_val <- data.frame(
  item = c("barcode_col_validation", "celltype_col_validation", "donor_col_validation", "dx_col_validation", "n_all_cells_validation", "n_mg_cells_validation"),
  value = c(barcode_col_val, celltype_col_val, donor_col_val, dx_col_val, ncol(counts), ncol(counts_val_mg))
)
fwrite(audit_val, file.path(cfg$outdir, "Validation_Audit.tsv"), sep = "\t")

# --------------------------------------------------
# 5. validation: build MG object and score all signatures
# --------------------------------------------------
msg("Building validation microglia Seurat object ...")
mg <- CreateSeuratObject(
  counts = counts_val_mg,
  meta.data = meta_val_mg,
  assay = cfg$assay,
  project = "VelmeshevMG",
  min.cells = 0,
  min.features = 0
)
DefaultAssay(mg) <- cfg$assay
mg <- NormalizeData(mg, verbose = FALSE)

features_use <- rownames(mg)

mapped_signature_audit <- list()
mapped_sig_list <- list()

for (nm in names(signature_list)) {
  genes_in <- signature_list[[nm]]
  genes_map <- harmonize_genes(genes_in, features_use)
  mapped_sig_list[[nm]] <- genes_map

  mapped_signature_audit[[length(mapped_signature_audit) + 1]] <- data.frame(
    signature_name = nm,
    n_input = length(genes_in),
    n_mapped = length(genes_map),
    mapped_gene_preview = paste(head(genes_map, 12), collapse = ", "),
    stringsAsFactors = FALSE
  )
}
mapped_signature_audit_df <- rbindlist(mapped_signature_audit, fill = TRUE)
fwrite(mapped_signature_audit_df, file.path(cfg$outdir, "Validation_MappedSignature_Audit.tsv"), sep = "\t")

valid_sigs <- names(mapped_sig_list)[sapply(mapped_sig_list, length) >= cfg$min_genes_mapped]
if (length(valid_sigs) == 0) stop("No signatures passed min_genes_mapped threshold.")

msg("Scoring ", length(valid_sigs), " discovery-cluster signatures in validation microglia ...")

for (sig_name in valid_sigs) {
  feat_list <- list(mapped_sig_list[[sig_name]])
  prefix <- paste0(sig_name, "_Score")
  mg <- AddModuleScore(mg, features = feat_list, name = prefix, assay = cfg$assay)
}

md <- mg@meta.data
md$cell_id <- colnames(mg)
md$donor <- as.character(md[[donor_col_val]])
md$dx <- normalize_dx(md[[dx_col_val]])

score_cols <- grep("_Score1$", colnames(md), value = TRUE)
if (length(score_cols) == 0) stop("No module-score columns found after AddModuleScore.")

for (sc in score_cols) {
  sig_name <- sub("_Score1$", "", sc)
  md[[sig_name]] <- md[[sc]]
}

keep_cols_cell <- c("cell_id", "donor", "dx", valid_sigs)
keep_cols_cell <- keep_cols_cell[keep_cols_cell %in% colnames(md)]
fwrite(md[, keep_cols_cell, drop = FALSE], file.path(cfg$outdir, "Validation_CellLevel_AllSignatureScores.tsv"), sep = "\t")

# donor-level aggregation wide
donor_df <- md %>%
  group_by(donor, dx) %>%
  summarise(across(all_of(valid_sigs), ~ mean(.x, na.rm = TRUE)), n_cells = n(), .groups = "drop") %>%
  as.data.frame()

fwrite(donor_df, file.path(cfg$outdir, "Validation_DonorLevel_AllSignatureScores.csv"))

# long donor-level
donor_long <- as.data.frame(rbindlist(lapply(valid_sigs, function(sig_name) {
  data.frame(
    donor = donor_df$donor,
    dx = donor_df$dx,
    n_cells = donor_df$n_cells,
    signature_name = sig_name,
    score = donor_df[[sig_name]],
    stringsAsFactors = FALSE
  )
}), fill = TRUE))

fwrite(donor_long, file.path(cfg$outdir, "Validation_DonorLevel_AllSignatureScores_long.tsv"), sep = "\t")

# --------------------------------------------------
# 6. validation: cluster signature ranking
# --------------------------------------------------
msg("Ranking validation results across all discovery-cluster signatures ...")
val_stats <- lapply(valid_sigs, function(sig_name) {
  sub <- donor_long %>% filter(signature_name == sig_name)
  p <- calc_wilcox_p(sub, "score", "dx")
  eff <- mean_diff_asd_minus_ctrl(sub, "score", "dx")
  data.frame(
    signature_name = sig_name,
    cluster = sub("^Cluster", "", sig_name),
    n_donors = length(unique(sub$donor)),
    mean_score_ASD = mean(sub$score[sub$dx == "ASD"], na.rm = TRUE),
    mean_score_Control = mean(sub$score[sub$dx == "Control"], na.rm = TRUE),
    effect_ASD_minus_Control = eff,
    wilcox_p = p,
    stringsAsFactors = FALSE
  )
})
val_stats_df <- rbindlist(val_stats, fill = TRUE)
val_stats_df$wilcox_fdr <- p.adjust(val_stats_df$wilcox_p, method = "BH")

prior_df <- val_stats_df %>%
  left_join(
    disc_stats_df %>% select(cluster, effect_ASD_minus_Control, wilcox_p, wilcox_fdr),
    by = "cluster",
    suffix = c("_validation", "_discovery")
  ) %>%
  left_join(
    signature_audit_df %>% select(cluster, signature_name, top_marker_preview),
    by = c("cluster", "signature_name")
  ) %>%
  left_join(
    mapped_signature_audit_df %>% select(signature_name, n_mapped),
    by = "signature_name"
  )

fwrite(prior_df, file.path(cfg$outdir, "Validation_AllClusterSignature_Stats.tsv"), sep = "\t")

# --------------------------------------------------
# 7. ranking plot
# --------------------------------------------------
plot_df <- prior_df %>%
  mutate(
    cluster_label = paste0("C", cluster),
    neglog10p = -log10(pmax(wilcox_p_validation, 1e-300))
  ) %>%
  arrange(wilcox_p_validation, desc(effect_ASD_minus_Control_validation))

plot_df$cluster_label <- factor(plot_df$cluster_label, levels = plot_df$cluster_label)

p_rank <- ggplot(plot_df, aes(x = effect_ASD_minus_Control_validation, y = cluster_label)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
  geom_point(aes(size = neglog10p)) +
  theme_classic(base_size = 12) +
  labs(
    title = "Validation of discovery-derived microglial cluster signatures",
    x = "Donor-level mean score difference (ASD - Control)",
    y = "Discovery cluster signature",
    size = "-log10(p)"
  )

ggsave(file.path(cfg$outdir, "Fig2B_AllClusterSignature_Ranking.pdf"), p_rank, width = 7.2, height = 6.2)

# --------------------------------------------------
# 8. choose top hit and optional reference hit
# --------------------------------------------------
if (cfg$prefer_positive_direction) {
  pos_df <- plot_df %>% filter(effect_ASD_minus_Control_validation > 0)
  if (nrow(pos_df) > 0) {
    top_hit <- pos_df %>% arrange(wilcox_p_validation, desc(effect_ASD_minus_Control_validation)) %>% slice(1)
  } else {
    top_hit <- plot_df %>% arrange(wilcox_p_validation, desc(effect_ASD_minus_Control_validation)) %>% slice(1)
  }
} else {
  top_hit <- plot_df %>% arrange(wilcox_p_validation, desc(abs(effect_ASD_minus_Control_validation))) %>% slice(1)
}

neg_df <- plot_df %>% filter(effect_ASD_minus_Control_validation < 0)
ref_hit <- if (nrow(neg_df) > 0) neg_df %>% arrange(wilcox_p_validation, effect_ASD_minus_Control_validation) %>% slice(1) else NULL

write.table(top_hit, file.path(cfg$outdir, "TopHit_Signature.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
if (!is.null(ref_hit)) {
  write.table(ref_hit, file.path(cfg$outdir, "ReferenceHit_Signature.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}

# --------------------------------------------------
# 9. donor-level boxplot for top hit
# --------------------------------------------------
make_boxplot <- function(sig_name, title_text, out_pdf) {
  sub <- donor_long %>% filter(signature_name == sig_name)
  sub$dx <- factor(sub$dx, levels = c("Control", "ASD"))

  p <- calc_wilcox_p(sub, "score", "dx")
  p_lab <- ifelse(is.na(p), "Wilcoxon p = NA", paste0("Wilcoxon p = ", signif(p, 3)))

  g <- ggplot(sub, aes(x = dx, y = score, fill = dx)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.75) +
    geom_jitter(width = 0.12, size = 2, alpha = 0.9) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 20, hjust = 1)
    ) +
    labs(
      title = title_text,
      subtitle = p_lab,
      x = NULL,
      y = "Mean module score per donor"
    )

  ggsave(out_pdf, g, width = 5.6, height = 4.8)
}

top_sig <- as.character(top_hit$signature_name[1])
make_boxplot(
  sig_name = top_sig,
  title_text = paste0("Top validation hit: ", top_sig),
  out_pdf = file.path(cfg$outdir, "Fig2C_TopHit_DonorLevel.pdf")
)

if (!is.null(ref_hit) && nrow(ref_hit) > 0) {
  ref_sig <- as.character(ref_hit$signature_name[1])
  make_boxplot(
    sig_name = ref_sig,
    title_text = paste0("Reference/opposite-direction hit: ", ref_sig),
    out_pdf = file.path(cfg$outdir, "Fig2D_ReferenceHit_DonorLevel.pdf")
  )
}

save_session_info(file.path(cfg$outdir, "sessionInfo_Figure2_Unbiased_AllClusters.txt"))
msg("Figure 2 unbiased all-cluster validation finished.")
