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
# Figure 2 enhanced script for Molecular Autism revision
# 修正版
# 目标：
# 1) 读取 Velmeshev 原始矩阵 + meta
# 2) 自动识别 microglia
# 3) 从 PsychENCODE Dual_Signatures_List.rds 读取风险/稳态 signature
# 4) 输出 donor-level Risk/Homeo scores
# 5) 输出 Fig2B / Fig2C donor-level PDF
# 6) 额外输出 signature audit 文件，用于核对 homeostatic 定义
# =========================================================

cfg <- list(
  raw_input   = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/Velmeshev_scRNA/rawMatrix.zip",
  meta_file   = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/Velmeshev_scRNA/meta.tsv",
  sig_rds     = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/PsychENCODE_Science_2024/Step3_Microglia/Dual_Signatures_List.rds",
  outdir      = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Figure2_Enhanced",
  donor_col   = NULL,   # 自动识别
  dx_col      = NULL,   # 自动识别
  barcode_col = NULL,   # 自动识别
  celltype_col = NULL,  # 自动识别
  assay       = "RNA",
  min_genes_per_set = 5
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

strip_bc <- function(x) {
  x <- as.character(x)
  x <- sub("-1$", "", x)
  x
}

read_counts_auto <- function(path) {
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    x <- readRDS(path)
    if (inherits(x, "Seurat")) {
      return(GetAssayData(x, assay = DefaultAssay(x), slot = "counts"))
    }
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
  find_first_col(
    meta,
    c("^donor$", "donor_id", "individual", "subject", "patient", "sample", "orig.ident", "library", "case_id")
  )
}

detect_dx_col <- function(meta, user_col = NULL) {
  if (!is.null(user_col) && user_col %in% names(meta)) return(user_col)
  find_first_col(
    meta,
    c("^dx$", "^group$", "diagnosis", "condition", "status", "phenotype", "disease")
  )
}

flatten_gene_sets <- function(x, path = "root") {
  out <- list()

  if (is.null(x)) return(out)

  if (is.factor(x)) x <- as.character(x)

  if (is.character(x)) {
    genes <- unique(x)
    genes <- genes[!is.na(genes) & nzchar(genes)]
    if (length(genes) > 0) out[[path]] <- genes
    return(out)
  }

  if (is.data.frame(x)) {
    cn <- names(x)
    cn_low <- tolower(cn)
    gene_col <- cn[cn_low %in% c("gene", "genes", "symbol", "gene_symbol", "feature", "features")]
    genes <- character(0)

    if (length(gene_col) > 0) {
      genes <- unique(as.character(x[[gene_col[1]]]))
    } else if (!is.null(rownames(x)) && !all(grepl("^\\d+$", rownames(x)))) {
      genes <- unique(rownames(x))
    }

    genes <- genes[!is.na(genes) & nzchar(genes)]
    if (length(genes) > 0) out[[path]] <- genes
    return(out)
  }

  if (is.list(x)) {
    nms <- names(x)
    if (is.null(nms)) nms <- paste0("idx", seq_along(x))
    for (i in seq_along(x)) {
      child_path <- paste(path, nms[i], sep = " / ")
      out <- c(out, flatten_gene_sets(x[[i]], child_path))
    }
  }

  out
}

pick_first_match <- function(gs, patterns) {
  if (length(gs) == 0) return(NULL)
  nms <- names(gs)
  for (p in patterns) {
    idx <- grepl(p, nms, ignore.case = TRUE)
    if (any(idx)) {
      ii <- which(idx)[1]
      return(list(path = nms[ii], genes = gs[[ii]]))
    }
  }
  NULL
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

# -----------------------------
# 1. 读取原始计数和 meta
# -----------------------------
msg("Reading counts ...")
counts <- read_counts_auto(cfg$raw_input)
msg("Counts dim: ", paste(dim(counts), collapse = " x "))

msg("Reading metadata ...")
meta <- fread(cfg$meta_file, data.table = FALSE)
meta <- as.data.frame(meta)

barcode_col <- detect_barcode_col(meta, colnames(counts), cfg$barcode_col)
if (is.null(barcode_col)) stop("Could not identify barcode column in metadata.")

meta2 <- align_meta_to_cells(meta, barcode_col, colnames(counts))
keep <- !is.na(meta2[[barcode_col]])
meta2 <- meta2[keep, , drop = FALSE]
counts <- counts[, rownames(meta2), drop = FALSE]

msg("Matched cells after meta alignment: ", ncol(counts))

celltype_col <- detect_celltype_col(meta2, cfg$celltype_col)
if (is.null(celltype_col)) stop("Could not identify a cell-type annotation column containing microglia labels.")

celltype_values <- as.character(meta2[[celltype_col]])
mg_idx <- grepl("microgl", celltype_values, ignore.case = TRUE) | grepl("^MG$", celltype_values, ignore.case = TRUE)

if (sum(mg_idx) == 0) {
  stop("No microglia cells matched in column: ", celltype_col)
}

meta_mg <- meta2[mg_idx, , drop = FALSE]
counts_mg <- counts[, rownames(meta_mg), drop = FALSE]

msg("Microglia cells retained: ", ncol(counts_mg))

donor_col <- detect_donor_col(meta_mg, cfg$donor_col)
dx_col    <- detect_dx_col(meta_mg, cfg$dx_col)

if (is.null(donor_col)) stop("Could not identify donor column. Please set cfg$donor_col manually.")
if (is.null(dx_col))    stop("Could not identify diagnosis column. Please set cfg$dx_col manually.")

audit_meta <- data.frame(
  item = c("barcode_col", "celltype_col", "donor_col", "dx_col", "n_all_cells", "n_mg_cells"),
  value = c(barcode_col, celltype_col, donor_col, dx_col, ncol(counts), ncol(counts_mg))
)
fwrite(audit_meta, file.path(cfg$outdir, "Figure2_Metadata_Audit.tsv"), sep = "\t")

# -----------------------------
# 2. 仅对 microglia 建对象
# -----------------------------
msg("Building Seurat object for microglia only ...")
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

# -----------------------------
# 3. 读取并审计 signatures
# -----------------------------
msg("Reading Dual_Signatures_List.rds ...")
sig_obj <- readRDS(cfg$sig_rds)
gs <- flatten_gene_sets(sig_obj)

if (length(gs) == 0) stop("No gene-set-like objects could be extracted from Dual_Signatures_List.rds")

gs_summary <- data.frame(
  path = names(gs),
  n_genes = sapply(gs, length),
  preview = sapply(gs, function(x) paste(head(unique(x), 8), collapse = ", "))
)
gs_summary <- gs_summary[order(-gs_summary$n_genes, gs_summary$path), ]
fwrite(gs_summary, file.path(cfg$outdir, "Signature_Audit_Summary.tsv"), sep = "\t")

risk_pick <- pick_first_match(gs, c("(^|/)Risk($|/)", "risk", "spp1", "cluster.?2", "cluster_2", "cluster 2"))
homeo_direct <- pick_first_match(gs, c("homeost", "homeo"))
cluster0_pick <- pick_first_match(gs, c("cluster.?0", "cluster_0", "cluster 0"))
cluster1_pick <- pick_first_match(gs, c("cluster.?1", "cluster_1", "cluster 1"))
cluster4_pick <- pick_first_match(gs, c("cluster.?4", "cluster_4", "cluster 4"))

if (is.null(risk_pick)) stop("Could not identify risk signature from Dual_Signatures_List.rds")

homeo_source <- NULL
homeo_genes_raw <- NULL
homeo_path <- NULL

if (!is.null(homeo_direct)) {
  homeo_source <- "direct_homeostatic_match"
  homeo_genes_raw <- homeo_direct$genes
  homeo_path <- homeo_direct$path
} else if (!is.null(cluster0_pick) && !is.null(cluster1_pick)) {
  homeo_source <- "cluster0_plus_cluster1_union"
  homeo_genes_raw <- union(cluster0_pick$genes, cluster1_pick$genes)
  homeo_path <- paste(cluster0_pick$path, " + ", cluster1_pick$path)
} else if (!is.null(cluster4_pick)) {
  homeo_source <- "cluster4_fallback"
  homeo_genes_raw <- cluster4_pick$genes
  homeo_path <- cluster4_pick$path
} else if (!is.null(cluster0_pick)) {
  homeo_source <- "cluster0_only_fallback"
  homeo_genes_raw <- cluster0_pick$genes
  homeo_path <- cluster0_pick$path
} else if (!is.null(cluster1_pick)) {
  homeo_source <- "cluster1_only_fallback"
  homeo_genes_raw <- cluster1_pick$genes
  homeo_path <- cluster1_pick$path
} else {
  stop("Could not identify any homeostatic candidate signature from Dual_Signatures_List.rds")
}

features_use <- rownames(mg)
risk_genes_use  <- harmonize_genes(risk_pick$genes, features_use)
homeo_genes_use <- harmonize_genes(homeo_genes_raw, features_use)

choice_tbl <- data.frame(
  signature = c("Risk", "Homeostatic"),
  source = c(risk_pick$path, homeo_source),
  path_or_note = c(risk_pick$path, homeo_path),
  n_raw_genes = c(length(unique(risk_pick$genes)), length(unique(homeo_genes_raw))),
  n_matched_genes = c(length(risk_genes_use), length(homeo_genes_use))
)
fwrite(choice_tbl, file.path(cfg$outdir, "Chosen_Signatures.tsv"), sep = "\t")

writeLines(unique(risk_genes_use), file.path(cfg$outdir, "Chosen_Risk_Genes.txt"))
writeLines(unique(homeo_genes_use), file.path(cfg$outdir, "Chosen_Homeostatic_Genes.txt"))

if (length(risk_genes_use) < cfg$min_genes_per_set) {
  stop("Matched risk genes < ", cfg$min_genes_per_set, ". Please inspect signature extraction.")
}
if (length(homeo_genes_use) < cfg$min_genes_per_set) {
  stop("Matched homeostatic genes < ", cfg$min_genes_per_set, ". Please inspect signature extraction.")
}

# -----------------------------
# 4. AddModuleScore
# -----------------------------
msg("Scoring risk/homeostatic signatures ...")
mg <- AddModuleScore(mg, features = list(risk_genes_use), name = "RiskScore", assay = cfg$assay)
mg <- AddModuleScore(mg, features = list(homeo_genes_use), name = "HomeoScore", assay = cfg$assay)

md <- mg@meta.data
md$cell_id <- colnames(mg)
md$donor   <- as.character(md[[donor_col]])
md$dx      <- normalize_dx(md[[dx_col]])
md$Risk    <- md$RiskScore1
md$Homeo   <- md$HomeoScore1

# -----------------------------
# 5. donor-level 汇总
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

fwrite(donor_df, file.path(cfg$outdir, "Validation_DonorLevel_Scores.csv"))

# 同时导出细胞级分数，改为未压缩 tsv，避免 zlib/gzip 报错
fwrite(
  md[, c("cell_id", "donor", "dx", "Risk", "Homeo"), drop = FALSE],
  file.path(cfg$outdir, "Validation_CellLevel_Scores.tsv"),
  sep = "\t"
)

# -----------------------------
# 6. donor-level 主图
# -----------------------------
order_dx <- c("Control", "ASD")
donor_df$dx <- factor(donor_df$dx, levels = unique(c(order_dx, as.character(donor_df$dx))))
p_risk <- plot_p_from_wilcox(donor_df, "Risk_mean", "dx")
p_homeo <- plot_p_from_wilcox(donor_df, "Homeo_mean", "dx")

p_lab_risk  <- ifelse(is.na(p_risk),  "Wilcoxon p = NA", paste0("Wilcoxon p = ", signif(p_risk, 3)))
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

ggsave(
  filename = file.path(cfg$outdir, "Fig2B_Risk_DonorLevel.pdf"),
  plot = g_risk, width = 5.5, height = 4.8
)

ggsave(
  filename = file.path(cfg$outdir, "Fig2C_Homeo_DonorLevel.pdf"),
  plot = g_homeo, width = 5.5, height = 4.8
)

# -----------------------------
# 7. supplementary 候选：per-cell violin
# -----------------------------
md$dx <- factor(md$dx, levels = unique(c(order_dx, as.character(md$dx))))

v_risk <- ggplot(md, aes(x = dx, y = Risk, fill = dx)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.8) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 20, hjust = 1)
  ) +
  labs(title = "Per-cell risk score (supplementary candidate)", x = NULL, y = "Module score")

v_homeo <- ggplot(md, aes(x = dx, y = Homeo, fill = dx)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.8) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 20, hjust = 1)
  ) +
  labs(title = "Per-cell homeostatic score (supplementary candidate)", x = NULL, y = "Module score")

ggsave(
  filename = file.path(cfg$outdir, "Supp_Fig2_Risk_PerCell_Violin.pdf"),
  plot = v_risk, width = 5.5, height = 4.8
)

ggsave(
  filename = file.path(cfg$outdir, "Supp_Fig2_Homeo_PerCell_Violin.pdf"),
  plot = v_homeo, width = 5.5, height = 4.8
)

save_session_info(file.path(cfg$outdir, "sessionInfo_Figure2.txt"))
msg("Figure 2 enhanced analysis finished.")
