#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(data.table)
  library(Matrix)
  library(ggplot2)
  library(edgeR)
  library(limma)
  library(statmod)
})

# ============================================================
# Package 4 v2
# Primary validation is NOT restricted to mapped clusters.
# Structure:
#   A. descriptive mapping
#   B. all-microglia frozen-score replication
#   C. score-defined candidate/reference-like cell replication
#   D. replication decision
# ============================================================

option_list <- list(
  make_option("--disc_rds", type = "character",
              default = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Package2_strict_input_fixed_from_original.rds"),
  make_option("--val_rds", type = "character",
              default = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/Velmeshev_scRNA/Velmeshev_Object.rds"),
  make_option("--outdir", type = "character",
              default = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Package4_CrossCohortValidation_Velmeshev_v2"),

  make_option("--disc_assay", type = "character", default = "RNA"),
  make_option("--val_assay", type = "character", default = "RNA"),

  make_option("--disc_cluster_col", type = "character", default = "cluster_use"),
  make_option("--disc_donor_col", type = "character", default = "individual_ID"),
  make_option("--disc_dx_col", type = "character", default = "Diagnosis"),

  make_option("--val_global_cluster_col", type = "character", default = "cluster"),
  make_option("--val_donor_col", type = "character", default = "sample"),
  make_option("--val_dx_col", type = "character", default = "diagnosis"),
  make_option("--val_microglia_label", type = "character", default = "Microglia"),

  make_option("--case_label", type = "character", default = "ASD"),
  make_option("--control_label", type = "character", default = "Control"),

  make_option("--disc_cluster0", type = "character", default = "0"),
  make_option("--disc_cluster2", type = "character", default = "2"),

  make_option("--min_cells_sample", type = "integer", default = 20),
  make_option("--min_cells_group", type = "integer", default = 10),

  make_option("--top_n_pos", type = "integer", default = 100),
  make_option("--top_n_neg", type = "integer", default = 100),

  make_option("--val_npcs", type = "integer", default = 20),
  make_option("--val_resolution", type = "double", default = 0.4),

  make_option("--mapping_corr_delta_threshold", type = "double", default = 0.05),
  make_option("--mapping_state_delta_threshold", type = "double", default = 0.00),

  make_option("--candidate_quantile", type = "double", default = 0.75),
  make_option("--reference_quantile", type = "double", default = 0.25),

  make_option("--alpha", type = "double", default = 0.05),
  make_option("--threads", type = "integer", default = 8)
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
logfile <- file.path(opt$outdir, "00_package4_v2_run.log")

log_msg <- function(...) {
  msg <- sprintf(...)
  line <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg)
  cat(line, "\n")
  cat(line, "\n", file = logfile, append = TRUE)
}

safe_fwrite <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  ok <- tryCatch({
    if (grepl("\\.gz$", file)) {
      data.table::fwrite(x, file = file, sep = "\t", quote = FALSE, na = "NA", compress = "gzip")
    } else {
      data.table::fwrite(x, file = file, sep = "\t", quote = FALSE, na = "NA")
    }
    TRUE
  }, error = function(e) FALSE)

  if (!ok) {
    tmp <- sub("\\.gz$", "", file)
    data.table::fwrite(x, file = tmp, sep = "\t", quote = FALSE, na = "NA")
    if (grepl("\\.gz$", file)) {
      system(sprintf("gzip -f %s", shQuote(tmp)))
    }
  }
}

normalize_dx <- function(x, case_label = "ASD", control_label = "Control") {
  x0 <- trimws(as.character(x))
  xl <- tolower(x0)
  out <- rep(NA_character_, length(x0))
  out[grepl("^asd$|autism|case", xl)] <- case_label
  out[grepl("^control$|^ctrl$|^ctl$|normal|neurotyp|^hc$", xl)] <- control_label
  out[x0 %in% c(case_label, control_label)] <- x0[x0 %in% c(case_label, control_label)]
  out
}

get_assay_mat <- function(obj, assay = "RNA", layer_name = c("counts", "data")) {
  layer_name <- match.arg(layer_name)
  mat <- tryCatch({
    GetAssayData(obj, assay = assay, slot = layer_name)
  }, error = function(e) {
    GetAssayData(obj, assay = assay, layer = layer_name)
  })
  mat
}

aggregate_counts_sparse <- function(count_mat, groups) {
  groups <- as.character(groups)
  keep <- !is.na(groups) & groups != ""
  count_mat <- count_mat[, keep, drop = FALSE]
  groups <- groups[keep]
  gf <- factor(groups, levels = unique(groups))
  mm <- Matrix::sparse.model.matrix(~ 0 + gf)
  colnames(mm) <- sub("^gf", "", colnames(mm))
  pb <- count_mat %*% mm
  pb
}

safe_wilcox <- function(df, value_col, dx_col = "dx", case_label = "ASD", control_label = "Control") {
  d <- copy(df)
  d <- d[get(dx_col) %in% c(control_label, case_label) & is.finite(get(value_col))]
  if (nrow(d) == 0) return(list(p = NA_real_, statistic = NA_real_))
  n_case <- sum(d[[dx_col]] == case_label)
  n_ctrl <- sum(d[[dx_col]] == control_label)
  if (n_case < 1 || n_ctrl < 1) return(list(p = NA_real_, statistic = NA_real_))
  wt <- tryCatch(
    wilcox.test(d[[value_col]] ~ factor(d[[dx_col]], levels = c(control_label, case_label)), exact = FALSE),
    error = function(e) NULL
  )
  if (is.null(wt)) return(list(p = NA_real_, statistic = NA_real_))
  list(p = unname(wt$p.value), statistic = unname(wt$statistic))
}

safe_lm_dx <- function(df, value_col, weight_col = NULL, dx_col = "dx", case_label = "ASD", control_label = "Control") {
  d <- copy(df)
  d <- d[get(dx_col) %in% c(control_label, case_label) & is.finite(get(value_col))]
  if (nrow(d) < 3) return(list(effect = NA_real_, p = NA_real_))
  d[, dx_factor := factor(get(dx_col), levels = c(control_label, case_label))]
  fit <- tryCatch({
    if (!is.null(weight_col) && weight_col %in% names(d)) {
      lm(reformulate("dx_factor", response = value_col), data = d, weights = pmax(d[[weight_col]], 1))
    } else {
      lm(reformulate("dx_factor", response = value_col), data = d)
    }
  }, error = function(e) NULL)
  if (is.null(fit)) return(list(effect = NA_real_, p = NA_real_))
  sm <- summary(fit)$coefficients
  coef_name <- grep("^dx_factor", rownames(sm), value = TRUE)
  if (length(coef_name) == 0) return(list(effect = NA_real_, p = NA_real_))
  list(effect = unname(sm[coef_name[1], "Estimate"]), p = unname(sm[coef_name[1], "Pr(>|t|)"]))
}

plot_metric_box <- function(df, metric_col, ylab_text, out_pdf, title_text = NULL,
                            dx_col = "dx", case_label = "ASD", control_label = "Control") {
  d <- copy(df)
  d <- d[get(dx_col) %in% c(control_label, case_label) & is.finite(get(metric_col))]
  if (nrow(d) == 0) return(invisible(NULL))
  d[, dx_factor := factor(get(dx_col), levels = c(control_label, case_label))]
  p <- ggplot(d, aes(x = dx_factor, y = get(metric_col))) +
    geom_violin(trim = FALSE, scale = "width", alpha = 0.25) +
    geom_boxplot(width = 0.18, outlier.shape = NA) +
    geom_jitter(width = 0.12, height = 0, alpha = 0.75, size = 1.6) +
    theme_bw(base_size = 12) +
    labs(x = NULL, y = ylab_text, title = title_text)
  ggsave(out_pdf, p, width = 4.8, height = 4.2)
}

run_edger_case_control <- function(pb_counts, sample_meta, out_prefix, dx_col = "dx",
                                   case_label = "ASD", control_label = "Control") {
  res_empty <- data.table(
    gene = rownames(pb_counts),
    logFC = NA_real_,
    logCPM = NA_real_,
    F = NA_real_,
    PValue = NA_real_,
    FDR = NA_real_,
    analysis = basename(out_prefix)
  )

  sample_meta <- as.data.table(sample_meta)
  if (ncol(pb_counts) != nrow(sample_meta)) stop("pseudobulk counts 列数与 sample_meta 行数不一致")

  keep_samples <- sample_meta[[dx_col]] %in% c(control_label, case_label)
  pb_counts <- pb_counts[, keep_samples, drop = FALSE]
  sample_meta <- sample_meta[keep_samples]

  if (ncol(pb_counts) < 4) {
    safe_fwrite(res_empty, paste0(out_prefix, ".tsv.gz"))
    return(res_empty)
  }

  dx_factor <- factor(sample_meta[[dx_col]], levels = c(control_label, case_label))
  if (sum(dx_factor == control_label) < 2 || sum(dx_factor == case_label) < 2) {
    safe_fwrite(res_empty, paste0(out_prefix, ".tsv.gz"))
    return(res_empty)
  }

  y <- edgeR::DGEList(counts = pb_counts)
  keep_genes <- edgeR::filterByExpr(y, group = dx_factor)
  y <- y[keep_genes, , keep.lib.sizes = FALSE]
  if (nrow(y) == 0) {
    safe_fwrite(res_empty, paste0(out_prefix, ".tsv.gz"))
    return(res_empty)
  }

  y <- edgeR::calcNormFactors(y)
  design <- model.matrix(~ dx_factor)
  y <- edgeR::estimateDisp(y, design, robust = TRUE)
  fit <- edgeR::glmQLFit(y, design, robust = TRUE)
  qlf <- edgeR::glmQLFTest(fit, coef = ncol(design))
  tt <- as.data.table(edgeR::topTags(qlf, n = Inf, sort.by = "none")$table, keep.rownames = "gene")
  tt[, analysis := basename(out_prefix)]
  safe_fwrite(tt, paste0(out_prefix, ".tsv.gz"))
  tt
}

run_paired_state_de <- function(pb_counts, sample_meta, out_prefix,
                                cluster0 = "0", cluster2 = "2") {
  res_empty <- data.table(
    gene = rownames(pb_counts),
    logFC = NA_real_,
    logCPM = NA_real_,
    F = NA_real_,
    PValue = NA_real_,
    FDR = NA_real_,
    analysis = basename(out_prefix)
  )

  meta <- as.data.table(sample_meta)
  if (ncol(pb_counts) != nrow(meta)) stop("paired state DE: pseudobulk counts 列数与 sample_meta 行数不一致")

  meta <- copy(meta)
  meta[, cluster := factor(as.character(cluster), levels = c(cluster0, cluster2))]
  meta[, sample := as.character(sample)]

  sample_keep <- meta[, .N, by = .(sample, cluster)][, dcast(.SD, sample ~ cluster, value.var = "N", fill = 0)]
  sample_keep <- sample_keep[get(cluster0) > 0 & get(cluster2) > 0, sample]
  meta <- meta[sample %in% sample_keep]

  if (nrow(meta) < 4 || uniqueN(meta$sample) < 2) {
    safe_fwrite(res_empty, paste0(out_prefix, ".tsv.gz"))
    return(res_empty)
  }

  pb_counts <- pb_counts[, meta$sample_id, drop = FALSE]

  sample_factor <- factor(meta$sample)
  cluster_factor <- factor(meta$cluster, levels = c(cluster0, cluster2))
  design <- model.matrix(~ sample_factor + cluster_factor)
  coef_idx <- grep("^cluster_factor", colnames(design))

  if (length(coef_idx) != 1) {
    safe_fwrite(res_empty, paste0(out_prefix, ".tsv.gz"))
    return(res_empty)
  }

  y <- edgeR::DGEList(counts = pb_counts)
  keep_genes <- edgeR::filterByExpr(y, design = design)
  y <- y[keep_genes, , keep.lib.sizes = FALSE]
  if (nrow(y) == 0) {
    safe_fwrite(res_empty, paste0(out_prefix, ".tsv.gz"))
    return(res_empty)
  }

  y <- edgeR::calcNormFactors(y)
  y <- edgeR::estimateDisp(y, design, robust = TRUE)
  fit <- edgeR::glmQLFit(y, design, robust = TRUE)
  qlf <- edgeR::glmQLFTest(fit, coef = coef_idx)

  tt <- as.data.table(edgeR::topTags(qlf, n = Inf, sort.by = "none")$table, keep.rownames = "gene")
  tt[, analysis := basename(out_prefix)]
  safe_fwrite(tt, paste0(out_prefix, ".tsv.gz"))
  tt
}

run_pseudobulk_sample_delta <- function(pb_counts, sample_meta, out_prefix,
                                        group0_label, group2_label,
                                        case_label = "ASD", control_label = "Control") {
  res_empty <- data.table(
    gene = rownames(pb_counts),
    logFC = NA_real_,
    AveExpr = NA_real_,
    t = NA_real_,
    PValue = NA_real_,
    FDR = NA_real_,
    analysis = basename(out_prefix)
  )

  meta <- as.data.table(sample_meta)
  if (ncol(pb_counts) != nrow(meta)) stop("paired sample-delta: pseudobulk counts 列数与 sample_meta 行数不一致")

  meta <- copy(meta)
  meta[, group := as.character(group)]
  meta[, dx := as.character(dx)]

  meta <- meta[group %in% c(group0_label, group2_label) & dx %in% c(control_label, case_label)]
  if (nrow(meta) == 0) {
    safe_fwrite(res_empty, paste0(out_prefix, ".tsv.gz"))
    return(res_empty)
  }

  pair_map <- dcast(
    meta[, .(sample, dx, group, sample_id)],
    sample + dx ~ group,
    value.var = "sample_id"
  )

  if (!all(c(group0_label, group2_label) %in% colnames(pair_map))) {
    safe_fwrite(res_empty, paste0(out_prefix, ".tsv.gz"))
    return(res_empty)
  }

  pair_map <- pair_map[!is.na(get(group0_label)) & !is.na(get(group2_label))]
  if (nrow(pair_map) < 4) {
    safe_fwrite(res_empty, paste0(out_prefix, ".tsv.gz"))
    return(res_empty)
  }

  n_case <- sum(pair_map$dx == case_label)
  n_ctrl <- sum(pair_map$dx == control_label)
  if (n_case < 2 || n_ctrl < 2) {
    safe_fwrite(res_empty, paste0(out_prefix, ".tsv.gz"))
    return(res_empty)
  }

  s0 <- pair_map[[group0_label]]
  s2 <- pair_map[[group2_label]]

  y <- edgeR::DGEList(counts = pb_counts)
  y <- edgeR::calcNormFactors(y)
  logcpm <- edgeR::cpm(y, log = TRUE, prior.count = 2, normalized.lib.sizes = TRUE)

  idx0 <- match(s0, colnames(logcpm))
  idx2 <- match(s2, colnames(logcpm))
  delta_mat <- logcpm[, idx2, drop = FALSE] - logcpm[, idx0, drop = FALSE]
  colnames(delta_mat) <- pair_map$sample

  keep_genes <- rowSums(is.finite(delta_mat)) == ncol(delta_mat)
  delta_mat <- delta_mat[keep_genes, , drop = FALSE]
  if (nrow(delta_mat) == 0) {
    safe_fwrite(res_empty, paste0(out_prefix, ".tsv.gz"))
    return(res_empty)
  }

  design_df <- data.frame(
    sample = pair_map$sample,
    dx_factor = factor(pair_map$dx, levels = c(control_label, case_label)),
    stringsAsFactors = FALSE
  )
  design <- model.matrix(~ dx_factor, data = design_df)

  fit <- limma::lmFit(delta_mat, design)
  fit <- limma::eBayes(fit, robust = TRUE)

  coef_idx <- grep("^dx_factor", colnames(design))
  if (length(coef_idx) != 1) {
    safe_fwrite(res_empty, paste0(out_prefix, ".tsv.gz"))
    return(res_empty)
  }

  tt <- as.data.table(
    limma::topTable(fit, coef = coef_idx, n = Inf, sort.by = "none"),
    keep.rownames = "gene"
  )
  if ("P.Value" %in% colnames(tt)) setnames(tt, "P.Value", "PValue")
  if ("adj.P.Val" %in% colnames(tt)) setnames(tt, "adj.P.Val", "FDR")
  tt[, analysis := basename(out_prefix)]

  safe_fwrite(tt, paste0(out_prefix, ".tsv.gz"))
  safe_fwrite(pair_map, paste0(out_prefix, ".pairs.tsv"))
  tt
}

extract_top_state_genes <- function(state_de, n_up = 100, n_down = 100) {
  d <- as.data.table(state_de)
  d <- d[is.finite(logFC) & is.finite(PValue)]
  if (nrow(d) == 0) {
    return(list(pos = character(), neg = character(), table = data.table()))
  }

  up <- d[logFC > 0][order(PValue, -abs(logFC))]
  down <- d[logFC < 0][order(PValue, -abs(logFC))]
  up <- up[seq_len(min(nrow(up), n_up))]
  down <- down[seq_len(min(nrow(down), n_down))]

  if (nrow(up) > 0) {
    up[, direction := "disc_cluster2_up"]
    up[, rank := seq_len(.N)]
  }
  if (nrow(down) > 0) {
    down[, direction := "disc_cluster0_up"]
    down[, rank := seq_len(.N)]
  }

  out <- rbindlist(list(
    up[, .(direction, rank, gene, logFC, PValue, FDR)],
    down[, .(direction, rank, gene, logFC, PValue, FDR)]
  ), fill = TRUE)

  list(pos = up$gene, neg = down$gene, table = out)
}

score_signature_matrix <- function(mat_data, pos_genes, neg_genes) {
  pos_genes <- intersect(pos_genes, rownames(mat_data))
  neg_genes <- intersect(neg_genes, rownames(mat_data))

  pos_score <- if (length(pos_genes) == 0) rep(NA_real_, ncol(mat_data)) else Matrix::colMeans(mat_data[pos_genes, , drop = FALSE])
  neg_score <- if (length(neg_genes) == 0) rep(NA_real_, ncol(mat_data)) else Matrix::colMeans(mat_data[neg_genes, , drop = FALSE])

  data.table(
    candidate_state_score = as.numeric(pos_score),
    reference_like_score = as.numeric(neg_score),
    candidate_minus_reference = as.numeric(pos_score - neg_score),
    candidate_over_reference = as.numeric((pos_score + 1e-6) / (neg_score + 1e-6))
  )
}

make_centroids <- function(mat_data, groups) {
  groups <- as.character(groups)
  keep <- !is.na(groups) & groups != ""
  mat_data <- mat_data[, keep, drop = FALSE]
  groups <- groups[keep]
  grp_levels <- unique(groups)
  out <- lapply(grp_levels, function(g) Matrix::rowMeans(mat_data[, groups == g, drop = FALSE]))
  names(out) <- grp_levels
  do.call(cbind, out)
}

make_cluster_mapping <- function(disc_centroids, val_centroids, pos_genes, neg_genes,
                                 disc_cluster0 = "0", disc_cluster2 = "2") {
  shared <- intersect(rownames(disc_centroids), rownames(val_centroids))
  if (length(shared) < 50) stop("discovery 与 validation centroid 共享基因过少，无法做 cluster mapping")

  disc_centroids <- disc_centroids[shared, , drop = FALSE]
  val_centroids <- val_centroids[shared, , drop = FALSE]

  pos_use <- intersect(pos_genes, shared)
  neg_use <- intersect(neg_genes, shared)

  dist_fun <- function(v, c) sqrt(mean((v - c)^2, na.rm = TRUE))

  res <- rbindlist(lapply(colnames(val_centroids), function(vc) {
    vv <- val_centroids[, vc]
    c0 <- disc_centroids[, disc_cluster0]
    c2 <- disc_centroids[, disc_cluster2]

    cor0 <- suppressWarnings(cor(vv, c0, method = "spearman"))
    cor2 <- suppressWarnings(cor(vv, c2, method = "spearman"))
    d0 <- dist_fun(vv, c0)
    d2 <- dist_fun(vv, c2)

    pos_score <- if (length(pos_use) > 0) mean(vv[pos_use], na.rm = TRUE) else NA_real_
    neg_score <- if (length(neg_use) > 0) mean(vv[neg_use], na.rm = TRUE) else NA_real_

    data.table(
      validation_cluster = vc,
      corr_to_disc_cluster0 = cor0,
      corr_to_disc_cluster2 = cor2,
      corr_delta_c2_minus_c0 = cor2 - cor0,
      dist_to_disc_cluster0 = d0,
      dist_to_disc_cluster2 = d2,
      nearest_centroid_assignment = ifelse(d2 < d0, "disc_cluster2_like", "disc_cluster0_like"),
      correlation_assignment = ifelse(cor2 > cor0, "disc_cluster2_like", "disc_cluster0_like"),
      candidate_state_score = pos_score,
      reference_like_score = neg_score,
      candidate_minus_reference = pos_score - neg_score
    )
  }))
  res
}

summarize_metric <- function(df, value_col, label, weight_col = NULL,
                             case_label = "ASD", control_label = "Control") {
  if (nrow(df) == 0 || !value_col %in% names(df)) {
    return(data.table(
      metric = label,
      mean_case = NA_real_,
      mean_control = NA_real_,
      effect_case_minus_control = NA_real_,
      wilcox_p = NA_real_,
      lm_p = NA_real_,
      n_case = 0L,
      n_control = 0L
    ))
  }
  w <- safe_wilcox(df, value_col, case_label = case_label, control_label = control_label)
  l <- safe_lm_dx(df, value_col, weight_col = weight_col, case_label = case_label, control_label = control_label)
  data.table(
    metric = label,
    mean_case = mean(df[dx == case_label, get(value_col)], na.rm = TRUE),
    mean_control = mean(df[dx == control_label, get(value_col)], na.rm = TRUE),
    effect_case_minus_control = l$effect,
    wilcox_p = w$p,
    lm_p = l$p,
    n_case = uniqueN(df[dx == case_label, sample]),
    n_control = uniqueN(df[dx == control_label, sample])
  )
}

try(data.table::setDTthreads(opt$threads), silent = TRUE)

log_msg("Reading discovery object: %s", opt$disc_rds)
disc_obj <- readRDS(opt$disc_rds)

log_msg("Reading validation object: %s", opt$val_rds)
val_obj <- readRDS(opt$val_rds)

DefaultAssay(disc_obj) <- if (opt$disc_assay %in% Assays(disc_obj)) opt$disc_assay else DefaultAssay(disc_obj)
DefaultAssay(val_obj)  <- if (opt$val_assay  %in% Assays(val_obj))  opt$val_assay  else DefaultAssay(val_obj)

disc_meta <- as.data.table(disc_obj@meta.data, keep.rownames = "cell")
val_meta0 <- as.data.table(val_obj@meta.data, keep.rownames = "cell")

disc_meta[, cluster := as.character(get(opt$disc_cluster_col))]
disc_meta[, sample := as.character(get(opt$disc_donor_col))]
disc_meta[, dx := normalize_dx(get(opt$disc_dx_col), opt$case_label, opt$control_label)]
disc_meta <- disc_meta[cluster %in% c(opt$disc_cluster0, opt$disc_cluster2) & !is.na(sample) & sample != "" & !is.na(dx)]

val_meta0[, global_cluster := as.character(get(opt$val_global_cluster_col))]
val_meta0[, sample := as.character(get(opt$val_donor_col))]
val_meta0[, dx := normalize_dx(get(opt$val_dx_col), opt$case_label, opt$control_label)]
val_meta0 <- val_meta0[!is.na(sample) & sample != "" & !is.na(dx) & !is.na(global_cluster) & global_cluster != ""]

safe_fwrite(
  data.table(
    disc_n_cells = nrow(disc_meta),
    disc_n_samples = uniqueN(disc_meta$sample),
    val_global_n_cells = nrow(val_meta0),
    val_global_n_samples = uniqueN(val_meta0$sample),
    val_global_n_clusters = uniqueN(val_meta0$global_cluster)
  ),
  file.path(opt$outdir, "01_run_parameters.tsv")
)

safe_fwrite(val_meta0[, .N, by = .(dx)][order(dx)], file.path(opt$outdir, "02_validation_global_dx_counts.tsv"))
safe_fwrite(val_meta0[global_cluster == opt$val_microglia_label, .N, by = .(dx)][order(dx)],
            file.path(opt$outdir, "03_validation_microglia_dx_counts.tsv"))

val_meta_mg <- val_meta0[global_cluster == opt$val_microglia_label]
if (nrow(val_meta_mg) == 0) stop("没有找到 validation microglia")

val_obj_mg <- subset(val_obj, cells = val_meta_mg$cell)
val_meta_mg <- as.data.table(val_obj_mg@meta.data, keep.rownames = "cell")
val_meta_mg[, sample := as.character(val_meta0$sample[match(cell, val_meta0$cell)])]
val_meta_mg[, dx := as.character(val_meta0$dx[match(cell, val_meta0$cell)])]

safe_fwrite(val_meta_mg[, .(cell, sample, dx)], file.path(opt$outdir, "04_validation_microglia_manifest.tsv.gz"))

DefaultAssay(val_obj_mg) <- opt$val_assay
val_obj_mg <- NormalizeData(val_obj_mg, assay = opt$val_assay, verbose = FALSE)
val_obj_mg <- FindVariableFeatures(val_obj_mg, assay = opt$val_assay, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
val_obj_mg <- ScaleData(val_obj_mg, assay = opt$val_assay, verbose = FALSE)
val_obj_mg <- RunPCA(val_obj_mg, assay = opt$val_assay, npcs = max(opt$val_npcs, 20), verbose = FALSE)
val_obj_mg <- FindNeighbors(val_obj_mg, dims = 1:opt$val_npcs, verbose = FALSE)
val_obj_mg <- FindClusters(val_obj_mg, resolution = opt$val_resolution, verbose = FALSE)
val_obj_mg <- RunUMAP(val_obj_mg, dims = 1:opt$val_npcs, umap.method = "uwot", metric = "cosine", verbose = FALSE)

val_meta <- as.data.table(val_obj_mg@meta.data, keep.rownames = "cell")
val_meta[, sample := as.character(val_meta_mg$sample[match(cell, val_meta_mg$cell)])]
val_meta[, dx := as.character(val_meta_mg$dx[match(cell, val_meta_mg$cell)])]
val_meta[, val_subcluster := as.character(seurat_clusters)]
safe_fwrite(val_meta[, .(cell, sample, dx, val_subcluster)], file.path(opt$outdir, "05_validation_microglia_reclustered_manifest.tsv.gz"))

disc_obj <- subset(disc_obj, cells = disc_meta$cell)
disc_counts <- get_assay_mat(disc_obj, assay = opt$disc_assay, layer_name = "counts")
disc_counts <- disc_counts[, disc_meta$cell, drop = FALSE]

disc_pair_meta <- disc_meta[, .N, by = .(sample, dx, cluster)]
disc_pair_meta <- disc_pair_meta[N >= opt$min_cells_group]
disc_pair_wide <- dcast(disc_pair_meta, sample + dx ~ cluster, value.var = "N", fill = 0)
disc_paired_samples <- disc_pair_wide[get(opt$disc_cluster0) > 0 & get(opt$disc_cluster2) > 0, sample]
disc_meta_pair <- disc_meta[sample %in% disc_paired_samples]

disc_pb_meta <- unique(disc_meta_pair[, .(sample_id = paste(sample, cluster, sep = "__"), sample, dx, cluster)])
disc_pb_n <- disc_meta_pair[, .(n_cells = .N), by = .(sample_id = paste(sample, cluster, sep = "__"))]
disc_pb_meta <- merge(disc_pb_meta, disc_pb_n, by = "sample_id", all.x = TRUE)
disc_pb_meta <- disc_pb_meta[n_cells >= opt$min_cells_group][order(sample, cluster)]
disc_pb_counts <- aggregate_counts_sparse(
  disc_counts[, disc_meta_pair$cell, drop = FALSE],
  paste(disc_meta_pair$sample, disc_meta_pair$cluster, sep = "__")
)
disc_pb_counts <- disc_pb_counts[, disc_pb_meta$sample_id, drop = FALSE]
safe_fwrite(disc_pb_meta, file.path(opt$outdir, "06_discovery_paired_pseudobulk_manifest.tsv"))

disc_state_de <- run_paired_state_de(
  disc_pb_counts,
  disc_pb_meta[, .(sample_id, sample, dx, cluster)],
  file.path(opt$outdir, "07_discovery_state_definition_cluster2_vs_cluster0_paired"),
  cluster0 = opt$disc_cluster0,
  cluster2 = opt$disc_cluster2
)

sig_info <- extract_top_state_genes(disc_state_de, opt$top_n_pos, opt$top_n_neg)
safe_fwrite(sig_info$table, file.path(opt$outdir, "08_discovery_frozen_state_signature.tsv"))

disc_data <- tryCatch(get_assay_mat(disc_obj, assay = opt$disc_assay, layer_name = "data"), error = function(e) NULL)
if (is.null(disc_data) || nrow(disc_data) == 0 || ncol(disc_data) == 0) {
  disc_obj <- NormalizeData(disc_obj, assay = opt$disc_assay, verbose = FALSE)
  disc_data <- get_assay_mat(disc_obj, assay = opt$disc_assay, layer_name = "data")
}
disc_data <- disc_data[, disc_meta$cell, drop = FALSE]
disc_centroids <- make_centroids(disc_data, disc_meta$cluster)
safe_fwrite(as.data.table(disc_centroids, keep.rownames = "gene"), file.path(opt$outdir, "09_discovery_centroids.tsv.gz"))

val_data <- tryCatch(get_assay_mat(val_obj_mg, assay = opt$val_assay, layer_name = "data"), error = function(e) NULL)
if (is.null(val_data) || nrow(val_data) == 0 || ncol(val_data) == 0) {
  val_obj_mg <- NormalizeData(val_obj_mg, assay = opt$val_assay, verbose = FALSE)
  val_data <- get_assay_mat(val_obj_mg, assay = opt$val_assay, layer_name = "data")
}
val_data <- val_data[, val_meta$cell, drop = FALSE]

val_centroids <- make_centroids(val_data, val_meta$val_subcluster)
safe_fwrite(as.data.table(val_centroids, keep.rownames = "gene"), file.path(opt$outdir, "10_validation_subcluster_centroids.tsv.gz"))

mapping_dt <- make_cluster_mapping(
  disc_centroids = disc_centroids,
  val_centroids = val_centroids,
  pos_genes = sig_info$pos,
  neg_genes = sig_info$neg,
  disc_cluster0 = opt$disc_cluster0,
  disc_cluster2 = opt$disc_cluster2
)
mapping_dt[, candidate_like_descriptive := correlation_assignment == "disc_cluster2_like" &
                                    corr_delta_c2_minus_c0 >= opt$mapping_corr_delta_threshold &
                                    candidate_minus_reference > opt$mapping_state_delta_threshold]
safe_fwrite(mapping_dt, file.path(opt$outdir, "11_validation_subcluster_mapping_descriptive.tsv"))

val_scores <- score_signature_matrix(val_data, sig_info$pos, sig_info$neg)
val_scores[, cell := val_meta$cell]
val_scores <- cbind(val_meta[, .(cell, sample, dx, val_subcluster)], val_scores)

cand_thr <- as.numeric(quantile(val_scores$candidate_minus_reference, probs = opt$candidate_quantile, na.rm = TRUE))
ref_thr  <- as.numeric(quantile(val_scores$candidate_minus_reference, probs = opt$reference_quantile, na.rm = TRUE))

val_scores[, score_defined_candidate := candidate_minus_reference >= cand_thr]
val_scores[, score_defined_reference := candidate_minus_reference <= ref_thr]
safe_fwrite(val_scores, file.path(opt$outdir, "12_validation_cell_scores.tsv.gz"))

safe_fwrite(
  data.table(candidate_threshold = cand_thr, reference_threshold = ref_thr,
             candidate_quantile = opt$candidate_quantile, reference_quantile = opt$reference_quantile),
  file.path(opt$outdir, "13_score_defined_thresholds.tsv")
)

sample_all <- val_scores[, .(
  n_cells = .N,
  candidate_state_score = mean(candidate_state_score, na.rm = TRUE),
  reference_like_score = mean(reference_like_score, na.rm = TRUE),
  candidate_minus_reference = mean(candidate_minus_reference, na.rm = TRUE),
  log2_candidate_over_reference = mean(log2(candidate_over_reference), na.rm = TRUE)
), by = .(sample, dx)]
sample_all <- sample_all[n_cells >= opt$min_cells_sample]
safe_fwrite(sample_all, file.path(opt$outdir, "14_validation_sample_scores_all_microglia.tsv"))

sample_prop <- val_scores[, .(
  n_total = .N,
  n_candidate = sum(score_defined_candidate),
  n_reference = sum(score_defined_reference)
), by = .(sample, dx)]
sample_prop <- sample_prop[n_total >= opt$min_cells_sample]
sample_prop[, prop_candidate := n_candidate / pmax(n_total, 1)]
sample_prop[, prop_reference := n_reference / pmax(n_total, 1)]
sample_prop[, log2_ratio_candidate_over_reference := log2((n_candidate + 0.5) / (n_reference + 0.5))]
safe_fwrite(sample_prop, file.path(opt$outdir, "15_validation_sample_score_defined_proportions.tsv"))

sample_candidate_scores <- val_scores[score_defined_candidate == TRUE, .(
  n_cells = .N,
  mean_score = mean(candidate_minus_reference, na.rm = TRUE)
), by = .(sample, dx)]
sample_reference_scores <- val_scores[score_defined_reference == TRUE, .(
  n_cells = .N,
  mean_score = mean(candidate_minus_reference, na.rm = TRUE)
), by = .(sample, dx)]

sample_delta <- merge(
  sample_candidate_scores[, .(sample, dx, cand_n = n_cells, cand_score = mean_score)],
  sample_reference_scores[, .(sample, dx, ref_n = n_cells, ref_score = mean_score)],
  by = c("sample", "dx"),
  all = FALSE
)
if (nrow(sample_delta) > 0) {
  sample_delta[, score_defined_shift := cand_score - ref_score]
  sample_delta[, paired_weight := pmin(cand_n, ref_n)]
}
safe_fwrite(sample_delta, file.path(opt$outdir, "16_validation_sample_score_defined_shift.tsv"))

rep_stats <- rbindlist(list(
  summarize_metric(sample_all, "candidate_minus_reference",
                   "all_microglia_candidate_minus_reference", "n_cells",
                   opt$case_label, opt$control_label),
  summarize_metric(sample_all, "candidate_state_score",
                   "all_microglia_candidate_state_score", "n_cells",
                   opt$case_label, opt$control_label),
  summarize_metric(sample_all, "reference_like_score",
                   "all_microglia_reference_like_score", "n_cells",
                   opt$case_label, opt$control_label),
  summarize_metric(sample_prop, "prop_candidate",
                   "score_defined_candidate_proportion", "n_total",
                   opt$case_label, opt$control_label),
  summarize_metric(sample_prop, "prop_reference",
                   "score_defined_reference_proportion", "n_total",
                   opt$case_label, opt$control_label),
  summarize_metric(sample_prop, "log2_ratio_candidate_over_reference",
                   "score_defined_log2_candidate_over_reference", "n_total",
                   opt$case_label, opt$control_label),
  summarize_metric(sample_delta, "score_defined_shift",
                   "score_defined_candidate_vs_reference_shift", "paired_weight",
                   opt$case_label, opt$control_label)
), fill = TRUE)
rep_stats[, p_adj_bh := p.adjust(wilcox_p, method = "BH")]
safe_fwrite(rep_stats, file.path(opt$outdir, "17_validation_v2_replication_stats.tsv"))

plot_metric_box(sample_all, "candidate_minus_reference",
                "Candidate minus reference",
                file.path(opt$outdir, "18_plot_all_microglia_candidate_minus_reference.pdf"),
                "All validation microglia: frozen score")
plot_metric_box(sample_prop, "prop_candidate",
                "Candidate-like cell proportion",
                file.path(opt$outdir, "19_plot_score_defined_candidate_proportion.pdf"),
                "Score-defined candidate-like cells")
if (nrow(sample_delta) > 0) {
  plot_metric_box(sample_delta, "score_defined_shift",
                  "Candidate minus reference shift",
                  file.path(opt$outdir, "20_plot_score_defined_shift.pdf"),
                  "Score-defined candidate vs reference shift")
}

val_counts_mg <- get_assay_mat(val_obj_mg, assay = opt$val_assay, layer_name = "counts")
val_counts_mg <- val_counts_mg[, val_meta$cell, drop = FALSE]

pb_all <- aggregate_counts_sparse(val_counts_mg, val_scores$sample)
pb_all_meta <- unique(val_scores[, .(sample_id = sample, sample, dx)])
pb_all_meta <- pb_all_meta[match(colnames(pb_all), sample_id)]
pb_all_meta[, n_cells := val_scores[, .N, by = sample][match(sample_id, sample), N]]
pb_all_meta <- pb_all_meta[n_cells >= opt$min_cells_sample]
pb_all <- pb_all[, pb_all_meta$sample_id, drop = FALSE]
safe_fwrite(pb_all_meta, file.path(opt$outdir, "21_validation_pseudobulk_manifest_all_microglia.tsv"))
run_edger_case_control(pb_all, pb_all_meta,
                       file.path(opt$outdir, "22_validation_pseudobulk_all_microglia_ASD_vs_Control"),
                       case_label = opt$case_label, control_label = opt$control_label)

val_meta_cand <- val_scores[score_defined_candidate == TRUE, .(cell, sample, dx)]
if (nrow(val_meta_cand) > 0) {
  pb_cand <- aggregate_counts_sparse(val_counts_mg[, val_meta_cand$cell, drop = FALSE], val_meta_cand$sample)
  pb_cand_meta <- unique(val_meta_cand[, .(sample_id = sample, sample, dx)])
  pb_cand_meta <- pb_cand_meta[match(colnames(pb_cand), sample_id)]
  pb_cand_meta[, n_cells := val_meta_cand[, .N, by = sample][match(sample_id, sample), N]]
  pb_cand_meta <- pb_cand_meta[n_cells >= opt$min_cells_group]
  pb_cand <- pb_cand[, pb_cand_meta$sample_id, drop = FALSE]
  safe_fwrite(pb_cand_meta, file.path(opt$outdir, "23_validation_pseudobulk_manifest_score_defined_candidate.tsv"))
  run_edger_case_control(pb_cand, pb_cand_meta,
                         file.path(opt$outdir, "24_validation_pseudobulk_score_defined_candidate_ASD_vs_Control"),
                         case_label = opt$case_label, control_label = opt$control_label)
}

val_meta_ref <- val_scores[score_defined_reference == TRUE, .(cell, sample, dx)]
if (nrow(val_meta_ref) > 0) {
  pb_ref <- aggregate_counts_sparse(val_counts_mg[, val_meta_ref$cell, drop = FALSE], val_meta_ref$sample)
  pb_ref_meta <- unique(val_meta_ref[, .(sample_id = sample, sample, dx)])
  pb_ref_meta <- pb_ref_meta[match(colnames(pb_ref), sample_id)]
  pb_ref_meta[, n_cells := val_meta_ref[, .N, by = sample][match(sample_id, sample), N]]
  pb_ref_meta <- pb_ref_meta[n_cells >= opt$min_cells_group]
  pb_ref <- pb_ref[, pb_ref_meta$sample_id, drop = FALSE]
  safe_fwrite(pb_ref_meta, file.path(opt$outdir, "25_validation_pseudobulk_manifest_score_defined_reference.tsv"))
  run_edger_case_control(pb_ref, pb_ref_meta,
                         file.path(opt$outdir, "26_validation_pseudobulk_score_defined_reference_ASD_vs_Control"),
                         case_label = opt$case_label, control_label = opt$control_label)
}

if (nrow(val_meta_cand) > 0 && nrow(val_meta_ref) > 0) {
  pair_meta <- rbindlist(list(
    val_meta_cand[, .(cell, sample, dx, group = "candidate")],
    val_meta_ref[, .(cell, sample, dx, group = "reference")]
  ), fill = TRUE)

  pair_group <- unique(pair_meta[, .(sample_id = paste(sample, group, sep = "__"), sample, dx, group)])
  pair_n <- pair_meta[, .(n_cells = .N), by = .(sample_id = paste(sample, group, sep = "__"))]
  pair_group <- merge(pair_group, pair_n, by = "sample_id", all.x = TRUE)
  pair_group <- pair_group[n_cells >= opt$min_cells_group][order(sample, group)]

  pb_pair <- aggregate_counts_sparse(
    val_counts_mg[, pair_meta$cell, drop = FALSE],
    paste(pair_meta$sample, pair_meta$group, sep = "__")
  )
  pb_pair <- pb_pair[, pair_group$sample_id, drop = FALSE]

  safe_fwrite(pair_group, file.path(opt$outdir, "27_validation_pseudobulk_manifest_score_defined_pairs.tsv"))
  run_pseudobulk_sample_delta(
    pb_pair, pair_group,
    file.path(opt$outdir, "28_validation_pseudobulk_score_defined_delta_candidate_minus_reference_ASD_vs_Control"),
    group0_label = "reference",
    group2_label = "candidate",
    case_label = opt$case_label,
    control_label = opt$control_label
  )
}

mapping_support <- FALSE
if (nrow(mapping_dt) > 0) {
  mapping_support <- any(mapping_dt$corr_delta_c2_minus_c0 >= opt$mapping_corr_delta_threshold, na.rm = TRUE)
}

primary_stats <- rep_stats[
  metric %in% c(
    "all_microglia_candidate_minus_reference",
    "score_defined_candidate_proportion",
    "score_defined_log2_candidate_over_reference",
    "score_defined_candidate_vs_reference_shift"
  )
]
primary_stats[, supported_direction := effect_case_minus_control > 0]
primary_stats[, significant := is.finite(wilcox_p) & wilcox_p < opt$alpha]
n_primary_supported <- sum(primary_stats$significant & primary_stats$supported_direction, na.rm = TRUE)

decision <- if (n_primary_supported >= 2) {
  "supported cross-cohort replication"
} else if (n_primary_supported >= 1 || mapping_support) {
  "partial cross-cohort replication"
} else {
  "limited cross-cohort replication"
}

decision_dt <- data.table(
  mapping_support_descriptive = mapping_support,
  n_primary_supported_metrics = n_primary_supported,
  replication_decision = decision
)
safe_fwrite(decision_dt, file.path(opt$outdir, "29_replication_decision_v2.tsv"))

decision_lines <- c(
  sprintf("Replication decision: %s", decision),
  sprintf("Mapping support (descriptive only): %s", mapping_support),
  sprintf("Number of primary supported metrics: %s", n_primary_supported),
  "Primary decision in v2 is based on all-microglia / score-defined replication, not mapped clusters alone."
)
writeLines(decision_lines, con = file.path(opt$outdir, "30_replication_decision_notes_v2.txt"))

capture.output(sessionInfo(), file = file.path(opt$outdir, "31_sessionInfo.txt"))
log_msg("Package 4 v2 completed. Output directory: %s", opt$outdir)
