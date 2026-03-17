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
# Package 4: cross-cohort validation
#   4A. state mapping
#   4B. donor-level replication
#   4C. replication decision
# ============================================================

option_list <- list(
  make_option("--disc_rds", type = "character",
              default = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Package2_strict_input_fixed_from_original.rds"),
  make_option("--val_rds", type = "character",
              default = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/Velmeshev_scRNA/Velmeshev_Object.rds"),
  make_option("--outdir", type = "character",
              default = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Package4_CrossCohortValidation_Velmeshev"),

  make_option("--disc_assay", type = "character", default = "RNA"),
  make_option("--val_assay", type = "character", default = "RNA"),

  make_option("--disc_cluster_col", type = "character", default = "AUTO"),
  make_option("--disc_donor_col", type = "character", default = "AUTO"),
  make_option("--disc_dx_col", type = "character", default = "AUTO"),

  make_option("--val_cluster_col", type = "character", default = "AUTO"),
  make_option("--val_donor_col", type = "character", default = "AUTO"),
  make_option("--val_dx_col", type = "character", default = "AUTO"),

  make_option("--case_label", type = "character", default = "ASD"),
  make_option("--control_label", type = "character", default = "Control"),

  make_option("--disc_cluster0", type = "character", default = "0"),
  make_option("--disc_cluster2", type = "character", default = "2"),

  make_option("--min_cells_donor", type = "integer", default = 20),
  make_option("--min_cells_cluster", type = "integer", default = 10),

  make_option("--top_n_pos", type = "integer", default = 100),
  make_option("--top_n_neg", type = "integer", default = 100),

  make_option("--mapping_corr_delta_threshold", type = "double", default = 0.05),
  make_option("--mapping_state_delta_threshold", type = "double", default = 0.00),
  make_option("--alpha", type = "double", default = 0.05),

  make_option("--gmt", type = "character", default = "NONE"),
  make_option("--min_pathway_size", type = "integer", default = 10),
  make_option("--max_pathway_size", type = "integer", default = 500),
  make_option("--top_pathways_pos", type = "integer", default = 10),
  make_option("--top_pathways_neg", type = "integer", default = 10),

  make_option("--threads", type = "integer", default = 8)
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
logfile <- file.path(opt$outdir, "00_package4_run.log")

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

pick_first_existing <- function(nms, candidates) {
  hit <- candidates[candidates %in% nms]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

resolve_column <- function(meta, user_col = "AUTO", type = c("cluster", "donor", "dx")) {
  type <- match.arg(type)
  nms <- colnames(meta)
  if (!identical(user_col, "AUTO")) {
    if (!user_col %in% nms) stop(sprintf("指定的 %s 列不存在: %s", type, user_col))
    return(user_col)
  }
  candidates <- switch(
    type,
    cluster = c("cluster_use", "seurat_clusters", "cluster", "Cluster", "subcluster", "ident"),
    donor   = c("donor_id", "donor", "individualID", "individual_id", "subject_id", "subject", "projid", "sample_id", "orig.ident", "sample"),
    dx      = c("Diagnosis", "diagnosis", "Dx", "dx", "Group", "group", "condition", "Condition", "phenotype", "Phenotype")
  )
  hit <- pick_first_existing(nms, candidates)
  if (is.na(hit)) {
    stop(sprintf("无法自动识别 %s 列。可用 metadata 列包括：%s", type, paste(nms, collapse = ", ")))
  }
  hit
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
  meta[, donor := as.character(donor)]

  donor_keep <- meta[, .N, by = .(donor, cluster)][, dcast(.SD, donor ~ cluster, value.var = "N", fill = 0)]
  donor_keep <- donor_keep[get(cluster0) > 0 & get(cluster2) > 0, donor]
  meta <- meta[donor %in% donor_keep]

  if (nrow(meta) < 4 || uniqueN(meta$donor) < 2) {
    safe_fwrite(res_empty, paste0(out_prefix, ".tsv.gz"))
    return(res_empty)
  }

  pb_counts <- pb_counts[, meta$sample_id, drop = FALSE]

  donor_factor <- factor(meta$donor)
  cluster_factor <- factor(meta$cluster, levels = c(cluster0, cluster2))
  design <- model.matrix(~ donor_factor + cluster_factor)
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

run_pseudobulk_donor_delta <- function(pb_counts, sample_meta, out_prefix,
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
  if (ncol(pb_counts) != nrow(meta)) stop("paired donor-delta: pseudobulk counts 列数与 sample_meta 行数不一致")

  meta <- copy(meta)
  meta[, group := as.character(group)]
  meta[, dx := as.character(dx)]

  meta <- meta[group %in% c(group0_label, group2_label) & dx %in% c(control_label, case_label)]
  if (nrow(meta) == 0) {
    safe_fwrite(res_empty, paste0(out_prefix, ".tsv.gz"))
    return(res_empty)
  }

  pair_map <- dcast(
    meta[, .(donor, dx, group, sample_id)],
    donor + dx ~ group,
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
  colnames(delta_mat) <- pair_map$donor

  keep_genes <- rowSums(is.finite(delta_mat)) == ncol(delta_mat)
  delta_mat <- delta_mat[keep_genes, , drop = FALSE]
  if (nrow(delta_mat) == 0) {
    safe_fwrite(res_empty, paste0(out_prefix, ".tsv.gz"))
    return(res_empty)
  }

  design_df <- data.frame(
    donor = pair_map$donor,
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
  safe_fwrite(pair_map, paste0(out_prefix, ".donor_pairs.tsv"))
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

  list(
    pos = up$gene,
    neg = down$gene,
    table = out
  )
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

select_candidate_reference_clusters <- function(mapping_dt,
                                                corr_delta_threshold = 0.05,
                                                state_delta_threshold = 0.00) {
  d <- copy(mapping_dt)

  d[, candidate_like := correlation_assignment == "disc_cluster2_like" &
                        corr_delta_c2_minus_c0 >= corr_delta_threshold &
                        candidate_minus_reference > state_delta_threshold]

  d[, reference_like := correlation_assignment == "disc_cluster0_like" &
                        corr_delta_c2_minus_c0 <= -corr_delta_threshold &
                        candidate_minus_reference < -state_delta_threshold]

  if (sum(d$candidate_like) == 0) {
    fallback_cand <- d[correlation_assignment == "disc_cluster2_like"][order(-corr_delta_c2_minus_c0, -candidate_minus_reference)]
    if (nrow(fallback_cand) > 0) d[validation_cluster == fallback_cand$validation_cluster[1], candidate_like := TRUE]
  }

  if (sum(d$reference_like) == 0) {
    fallback_ref <- d[correlation_assignment == "disc_cluster0_like"][order(corr_delta_c2_minus_c0, candidate_minus_reference)]
    if (nrow(fallback_ref) > 0) d[validation_cluster == fallback_ref$validation_cluster[1], reference_like := TRUE]
  }

  d
}

parse_gmt <- function(gmt_file, min_size = 10, max_size = 500) {
  lines <- readLines(gmt_file)
  gmt <- lapply(lines, function(x) {
    parts <- strsplit(x, "\t", fixed = TRUE)[[1]]
    genes <- unique(parts[-c(1, 2)])
    genes
  })
  names(gmt) <- vapply(lines, function(x) strsplit(x, "\t", fixed = TRUE)[[1]][1], FUN.VALUE = character(1))
  gmt <- gmt[lengths(gmt) >= min_size & lengths(gmt) <= max_size]
  gmt
}

select_discovery_pathways <- function(state_de, gmt_list, top_pos = 10, top_neg = 10) {
  d <- copy(as.data.table(state_de))
  d <- d[is.finite(logFC) & is.finite(PValue)]
  d[, signed_stat := sign(logFC) * -log10(pmax(PValue, 1e-300))]
  ss <- d$signed_stat
  names(ss) <- d$gene

  path_dt <- rbindlist(lapply(names(gmt_list), function(pw) {
    genes <- intersect(gmt_list[[pw]], names(ss))
    if (length(genes) == 0) return(NULL)
    data.table(
      pathway = pw,
      n_overlap = length(genes),
      mean_signed_stat = mean(ss[genes], na.rm = TRUE),
      mean_logFC = mean(d[match(genes, gene), logFC], na.rm = TRUE)
    )
  }), fill = TRUE)

  if (nrow(path_dt) == 0) return(list(table = data.table(), pos = character(), neg = character()))

  pos <- path_dt[mean_signed_stat > 0][order(-mean_signed_stat)][seq_len(min(.N, top_pos))]
  neg <- path_dt[mean_signed_stat < 0][order(mean_signed_stat)][seq_len(min(.N, top_neg))]
  list(table = path_dt[order(-abs(mean_signed_stat))], pos = pos$pathway, neg = neg$pathway)
}

score_pathways_on_pseudobulk <- function(pb_logcpm, pathway_names, gmt_list) {
  if (length(pathway_names) == 0) return(data.table())
  genes_all <- rownames(pb_logcpm)
  zmat <- t(scale(t(pb_logcpm)))
  zmat[!is.finite(zmat)] <- NA_real_

  out <- rbindlist(lapply(pathway_names, function(pw) {
    genes <- intersect(gmt_list[[pw]], genes_all)
    if (length(genes) == 0) return(NULL)
    data.table(
      pathway = pw,
      sample_id = colnames(pb_logcpm),
      pathway_score = as.numeric(colMeans(t(zmat[genes, , drop = FALSE]), na.rm = TRUE))
    )
  }), fill = TRUE)

  out
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
      n_case_donor = 0L,
      n_control_donor = 0L
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
    n_case_donor = uniqueN(df[dx == case_label, donor]),
    n_control_donor = uniqueN(df[dx == control_label, donor])
  )
}

# ============================================================
# 0. Read objects
# ============================================================
try(data.table::setDTthreads(opt$threads), silent = TRUE)

log_msg("Reading discovery object: %s", opt$disc_rds)
if (!file.exists(opt$disc_rds)) stop("discovery RDS 不存在")
disc_obj <- readRDS(opt$disc_rds)

log_msg("Reading validation object: %s", opt$val_rds)
if (!file.exists(opt$val_rds)) stop("validation RDS 不存在")
val_obj <- readRDS(opt$val_rds)

disc_assay <- if (opt$disc_assay %in% Assays(disc_obj)) opt$disc_assay else DefaultAssay(disc_obj)
val_assay <- if (opt$val_assay %in% Assays(val_obj)) opt$val_assay else DefaultAssay(val_obj)
DefaultAssay(disc_obj) <- disc_assay
DefaultAssay(val_obj) <- val_assay

disc_meta <- as.data.table(disc_obj@meta.data, keep.rownames = "cell")
val_meta  <- as.data.table(val_obj@meta.data,  keep.rownames = "cell")

log_msg("Discovery metadata columns: %s", paste(colnames(disc_meta), collapse = ", "))
log_msg("Validation metadata columns: %s", paste(colnames(val_meta), collapse = ", "))

disc_cluster_col <- resolve_column(disc_meta, opt$disc_cluster_col, "cluster")
disc_donor_col   <- resolve_column(disc_meta, opt$disc_donor_col, "donor")
disc_dx_col      <- resolve_column(disc_meta, opt$disc_dx_col, "dx")

val_cluster_col <- resolve_column(val_meta, opt$val_cluster_col, "cluster")
val_donor_col   <- resolve_column(val_meta, opt$val_donor_col, "donor")
val_dx_col      <- resolve_column(val_meta, opt$val_dx_col, "dx")

disc_meta[, cluster := as.character(get(disc_cluster_col))]
disc_meta[, donor := as.character(get(disc_donor_col))]
disc_meta[, dx := normalize_dx(get(disc_dx_col), opt$case_label, opt$control_label)]

val_meta[, cluster := as.character(get(val_cluster_col))]
val_meta[, donor := as.character(get(val_donor_col))]
val_meta[, dx := normalize_dx(get(val_dx_col), opt$case_label, opt$control_label)]

disc_keep <- c(opt$disc_cluster0, opt$disc_cluster2)
disc_meta <- disc_meta[cluster %in% disc_keep & !is.na(donor) & donor != "" & !is.na(dx)]
val_meta  <- val_meta[!is.na(donor) & donor != "" & !is.na(dx) & !is.na(cluster) & cluster != ""]

if (nrow(disc_meta) == 0) stop("discovery 过滤后没有可用细胞")
if (nrow(val_meta) == 0) stop("validation 过滤后没有可用细胞")

disc_obj <- subset(disc_obj, cells = disc_meta$cell)
val_obj  <- subset(val_obj,  cells = val_meta$cell)

disc_meta <- as.data.table(disc_obj@meta.data, keep.rownames = "cell")
disc_meta[, cluster := as.character(get(disc_cluster_col))]
disc_meta[, donor := as.character(get(disc_donor_col))]
disc_meta[, dx := normalize_dx(get(disc_dx_col), opt$case_label, opt$control_label)]
disc_meta <- disc_meta[cluster %in% disc_keep & !is.na(donor) & donor != "" & !is.na(dx)]

val_meta <- as.data.table(val_obj@meta.data, keep.rownames = "cell")
val_meta[, cluster := as.character(get(val_cluster_col))]
val_meta[, donor := as.character(get(val_donor_col))]
val_meta[, dx := normalize_dx(get(val_dx_col), opt$case_label, opt$control_label)]
val_meta <- val_meta[!is.na(donor) & donor != "" & !is.na(dx) & !is.na(cluster) & cluster != ""]

safe_fwrite(
  data.table(
    disc_rds = opt$disc_rds,
    val_rds = opt$val_rds,
    disc_assay = disc_assay,
    val_assay = val_assay,
    disc_cluster_col = disc_cluster_col,
    disc_donor_col = disc_donor_col,
    disc_dx_col = disc_dx_col,
    val_cluster_col = val_cluster_col,
    val_donor_col = val_donor_col,
    val_dx_col = val_dx_col,
    disc_n_cells = nrow(disc_meta),
    disc_n_donors = uniqueN(disc_meta$donor),
    val_n_cells = nrow(val_meta),
    val_n_donors = uniqueN(val_meta$donor),
    val_n_clusters = uniqueN(val_meta$cluster)
  ),
  file.path(opt$outdir, "01_run_parameters.tsv")
)

safe_fwrite(disc_meta[, .(cell, donor, dx, cluster)], file.path(opt$outdir, "02_discovery_cell_manifest.tsv.gz"))
safe_fwrite(val_meta[, .(cell, donor, dx, cluster)],  file.path(opt$outdir, "03_validation_cell_manifest.tsv.gz"))

# ============================================================
# 1. Discovery frozen reference
# ============================================================
log_msg("Preparing discovery paired pseudobulk for cluster %s vs %s", opt$disc_cluster0, opt$disc_cluster2)

disc_counts <- get_assay_mat(disc_obj, assay = disc_assay, layer_name = "counts")
disc_counts <- disc_counts[, disc_meta$cell, drop = FALSE]

disc_pair_meta <- disc_meta[, .N, by = .(donor, dx, cluster)]
disc_pair_meta <- disc_pair_meta[N >= opt$min_cells_cluster]

disc_pair_wide <- dcast(disc_pair_meta, donor + dx ~ cluster, value.var = "N", fill = 0)
disc_paired_donors <- disc_pair_wide[get(opt$disc_cluster0) > 0 & get(opt$disc_cluster2) > 0, donor]
disc_meta_pair <- disc_meta[donor %in% disc_paired_donors]

disc_pb_meta <- unique(disc_meta_pair[, .(sample_id = paste(donor, cluster, sep = "__"), donor, dx, cluster)])
disc_pb_n <- disc_meta_pair[, .(n_cells = .N), by = .(sample_id = paste(donor, cluster, sep = "__"))]
disc_pb_meta <- merge(disc_pb_meta, disc_pb_n, by = "sample_id", all.x = TRUE)
disc_pb_meta <- disc_pb_meta[n_cells >= opt$min_cells_cluster][order(donor, cluster)]

disc_pb_counts <- aggregate_counts_sparse(
  disc_counts[, disc_meta_pair$cell, drop = FALSE],
  paste(disc_meta_pair$donor, disc_meta_pair$cluster, sep = "__")
)
disc_pb_counts <- disc_pb_counts[, disc_pb_meta$sample_id, drop = FALSE]
safe_fwrite(disc_pb_meta, file.path(opt$outdir, "04_discovery_paired_pseudobulk_manifest.tsv"))

log_msg("Running discovery paired state DE")
disc_state_de <- run_paired_state_de(
  disc_pb_counts,
  disc_pb_meta,
  file.path(opt$outdir, "05_discovery_state_definition_cluster2_vs_cluster0_paired"),
  cluster0 = opt$disc_cluster0,
  cluster2 = opt$disc_cluster2
)

sig_info <- extract_top_state_genes(disc_state_de, opt$top_n_pos, opt$top_n_neg)
safe_fwrite(sig_info$table, file.path(opt$outdir, "06_discovery_frozen_state_signature.tsv"))

log_msg("Building discovery centroids")
disc_data <- tryCatch(get_assay_mat(disc_obj, assay = disc_assay, layer_name = "data"), error = function(e) NULL)
if (is.null(disc_data) || nrow(disc_data) == 0 || ncol(disc_data) == 0) {
  disc_obj <- NormalizeData(disc_obj, assay = disc_assay, verbose = FALSE)
  disc_data <- get_assay_mat(disc_obj, assay = disc_assay, layer_name = "data")
}
disc_data <- disc_data[, disc_meta$cell, drop = FALSE]
disc_centroids <- make_centroids(disc_data, disc_meta$cluster)
safe_fwrite(as.data.table(disc_centroids, keep.rownames = "gene"), file.path(opt$outdir, "07_discovery_centroids.tsv.gz"))

# ============================================================
# 2. Validation cluster mapping
# ============================================================
log_msg("Preparing validation normalized data and centroids")
val_data <- tryCatch(get_assay_mat(val_obj, assay = val_assay, layer_name = "data"), error = function(e) NULL)
if (is.null(val_data) || nrow(val_data) == 0 || ncol(val_data) == 0) {
  val_obj <- NormalizeData(val_obj, assay = val_assay, verbose = FALSE)
  val_data <- get_assay_mat(val_obj, assay = val_assay, layer_name = "data")
}
val_data <- val_data[, val_meta$cell, drop = FALSE]

val_centroids <- make_centroids(val_data, val_meta$cluster)
safe_fwrite(as.data.table(val_centroids, keep.rownames = "gene"), file.path(opt$outdir, "08_validation_cluster_centroids.tsv.gz"))

log_msg("Running validation state mapping")
mapping_dt <- make_cluster_mapping(
  disc_centroids = disc_centroids,
  val_centroids = val_centroids,
  pos_genes = sig_info$pos,
  neg_genes = sig_info$neg,
  disc_cluster0 = opt$disc_cluster0,
  disc_cluster2 = opt$disc_cluster2
)

mapping_dt <- select_candidate_reference_clusters(
  mapping_dt,
  corr_delta_threshold = opt$mapping_corr_delta_threshold,
  state_delta_threshold = opt$mapping_state_delta_threshold
)
mapping_dt <- mapping_dt[order(-candidate_like, reference_like, -corr_delta_c2_minus_c0, -candidate_minus_reference)]
safe_fwrite(mapping_dt, file.path(opt$outdir, "09_validation_cluster_mapping.tsv"))

candidate_clusters <- mapping_dt[candidate_like == TRUE, validation_cluster]
reference_clusters <- mapping_dt[reference_like == TRUE, validation_cluster]

log_msg("Mapped candidate-like validation clusters: %s", paste(candidate_clusters, collapse = ", "))
log_msg("Mapped reference-like validation clusters: %s", paste(reference_clusters, collapse = ", "))

if (nrow(mapping_dt) > 0) {
  p_map <- ggplot(mapping_dt, aes(x = corr_to_disc_cluster0, y = corr_to_disc_cluster2, label = validation_cluster)) +
    geom_point(aes(shape = candidate_like, size = candidate_minus_reference)) +
    theme_bw(base_size = 11) +
    labs(x = "Correlation to discovery cluster 0",
         y = "Correlation to discovery cluster 2",
         title = "Validation cluster mapping to discovery retained reference")
  ggsave(file.path(opt$outdir, "10_validation_cluster_mapping_scatter.pdf"), p_map, width = 6.2, height = 5.2)
}

# ============================================================
# 3. Validation cell scores and donor-level replication
# ============================================================
log_msg("Scoring validation cells with frozen discovery signature")
val_scores <- score_signature_matrix(val_data, sig_info$pos, sig_info$neg)
val_scores[, cell := val_meta$cell]
val_scores <- cbind(val_meta[, .(cell, donor, dx, cluster)], val_scores)
val_scores[, is_candidate_cluster := cluster %in% candidate_clusters]
val_scores[, is_reference_cluster := cluster %in% reference_clusters]
val_scores[, in_retained_mapped_space := is_candidate_cluster | is_reference_cluster]
safe_fwrite(val_scores, file.path(opt$outdir, "11_validation_cell_scores.tsv.gz"))

# donor-level candidate
donor_candidate <- val_scores[is_candidate_cluster == TRUE, .(
  n_cells = .N,
  candidate_state_score = mean(candidate_state_score, na.rm = TRUE),
  reference_like_score = mean(reference_like_score, na.rm = TRUE),
  candidate_minus_reference = mean(candidate_minus_reference, na.rm = TRUE),
  log2_candidate_over_reference = mean(log2(candidate_over_reference), na.rm = TRUE)
), by = .(donor, dx)]

if (nrow(donor_candidate) > 0) donor_candidate <- donor_candidate[n_cells >= opt$min_cells_cluster]
safe_fwrite(donor_candidate, file.path(opt$outdir, "12_validation_donor_scores_candidateclusters.tsv"))

# donor-level reference
donor_reference <- val_scores[is_reference_cluster == TRUE, .(
  n_cells = .N,
  candidate_state_score = mean(candidate_state_score, na.rm = TRUE),
  reference_like_score = mean(reference_like_score, na.rm = TRUE),
  candidate_minus_reference = mean(candidate_minus_reference, na.rm = TRUE),
  log2_candidate_over_reference = mean(log2(candidate_over_reference), na.rm = TRUE)
), by = .(donor, dx)]

if (nrow(donor_reference) > 0) donor_reference <- donor_reference[n_cells >= opt$min_cells_cluster]
safe_fwrite(donor_reference, file.path(opt$outdir, "13_validation_donor_scores_referenceclusters.tsv"))

# donor-level retained mapped space
donor_retained <- val_scores[in_retained_mapped_space == TRUE, .(
  n_cells = .N,
  candidate_state_score = mean(candidate_state_score, na.rm = TRUE),
  reference_like_score = mean(reference_like_score, na.rm = TRUE),
  candidate_minus_reference = mean(candidate_minus_reference, na.rm = TRUE),
  log2_candidate_over_reference = mean(log2(candidate_over_reference), na.rm = TRUE)
), by = .(donor, dx)]

if (nrow(donor_retained) > 0) donor_retained <- donor_retained[n_cells >= opt$min_cells_donor]
safe_fwrite(donor_retained, file.path(opt$outdir, "14_validation_donor_scores_retainedspace.tsv"))

# donor-level candidate vs reference shift
donor_shift <- merge(
  donor_candidate[, .(donor, dx, cand_n = n_cells, cand_score = candidate_minus_reference)],
  donor_reference[, .(donor, dx, ref_n = n_cells, ref_score = candidate_minus_reference)],
  by = c("donor", "dx"),
  all = FALSE
)
if (nrow(donor_shift) > 0) {
  donor_shift[, shift_candidate_minus_reference := cand_score - ref_score]
  donor_shift[, paired_weight := pmin(cand_n, ref_n)]
}
safe_fwrite(donor_shift, file.path(opt$outdir, "15_validation_donor_candidate_vs_reference_shift.tsv"))

# donor-level abundance within retained mapped space
donor_abund <- val_scores[in_retained_mapped_space == TRUE, .(
  n_total_retained = .N,
  n_candidate = sum(is_candidate_cluster),
  n_reference = sum(is_reference_cluster)
), by = .(donor, dx)]

if (nrow(donor_abund) > 0) {
  donor_abund <- donor_abund[n_total_retained >= opt$min_cells_donor]
  donor_abund[, prop_candidate := n_candidate / pmax(n_total_retained, 1)]
  donor_abund[, log2_ratio_candidate_over_reference := log2((n_candidate + 0.5) / (n_reference + 0.5))]
}
safe_fwrite(donor_abund, file.path(opt$outdir, "16_validation_candidate_cluster_abundance.tsv"))

rep_stats <- rbindlist(list(
  summarize_metric(donor_candidate, "candidate_minus_reference",
                   "candidateclusters_candidate_minus_reference", "n_cells",
                   opt$case_label, opt$control_label),
  summarize_metric(donor_reference, "candidate_minus_reference",
                   "referenceclusters_candidate_minus_reference", "n_cells",
                   opt$case_label, opt$control_label),
  summarize_metric(donor_retained, "candidate_minus_reference",
                   "retainedspace_candidate_minus_reference", "n_cells",
                   opt$case_label, opt$control_label),
  summarize_metric(donor_shift, "shift_candidate_minus_reference",
                   "candidate_vs_reference_shift", "paired_weight",
                   opt$case_label, opt$control_label),
  summarize_metric(donor_abund, "prop_candidate",
                   "candidate_cluster_proportion_within_retainedspace", "n_total_retained",
                   opt$case_label, opt$control_label),
  summarize_metric(donor_abund, "log2_ratio_candidate_over_reference",
                   "log2_ratio_candidate_over_reference", "n_total_retained",
                   opt$case_label, opt$control_label)
), fill = TRUE)
rep_stats[, p_adj_bh := p.adjust(wilcox_p, method = "BH")]
safe_fwrite(rep_stats, file.path(opt$outdir, "17_validation_donor_replication_stats.tsv"))

if (nrow(donor_candidate) > 0) {
  plot_metric_box(donor_candidate, "candidate_minus_reference",
                  "Candidate minus reference",
                  file.path(opt$outdir, "18_plot_candidateclusters_candidate_minus_reference.pdf"),
                  "Validation donor-level score: candidate clusters")
}
if (nrow(donor_reference) > 0) {
  plot_metric_box(donor_reference, "candidate_minus_reference",
                  "Candidate minus reference",
                  file.path(opt$outdir, "19_plot_referenceclusters_candidate_minus_reference.pdf"),
                  "Validation donor-level score: reference clusters")
}
if (nrow(donor_retained) > 0) {
  plot_metric_box(donor_retained, "candidate_minus_reference",
                  "Candidate minus reference",
                  file.path(opt$outdir, "20_plot_retainedspace_candidate_minus_reference.pdf"),
                  "Validation donor-level score: retained mapped space")
}
if (nrow(donor_shift) > 0) {
  plot_metric_box(donor_shift, "shift_candidate_minus_reference",
                  "Candidate minus reference shift",
                  file.path(opt$outdir, "21_plot_candidate_vs_reference_shift.pdf"),
                  "Validation donor-level candidate vs reference shift")
}
if (nrow(donor_abund) > 0) {
  plot_metric_box(donor_abund, "prop_candidate",
                  "Candidate cluster proportion",
                  file.path(opt$outdir, "22_plot_candidate_cluster_proportion.pdf"),
                  "Validation donor-level candidate abundance")
}

# ============================================================
# 4. Validation pseudobulk replication
# ============================================================
log_msg("Running validation pseudobulk replication")
val_counts <- get_assay_mat(val_obj, assay = val_assay, layer_name = "counts")
val_counts <- val_counts[, val_meta$cell, drop = FALSE]

# retained mapped space
val_meta_retained <- val_meta[cluster %in% c(candidate_clusters, reference_clusters)]
if (nrow(val_meta_retained) > 0) {
  pb_ret <- aggregate_counts_sparse(val_counts[, val_meta_retained$cell, drop = FALSE], val_meta_retained$donor)
  pb_ret_meta <- unique(val_meta_retained[, .(sample_id = donor, donor, dx)])
  pb_ret_meta <- pb_ret_meta[match(colnames(pb_ret), sample_id)]
  pb_ret_meta[, n_cells := val_meta_retained[, .N, by = donor][match(sample_id, donor), N]]
  pb_ret_meta <- pb_ret_meta[n_cells >= opt$min_cells_donor]
  pb_ret <- pb_ret[, pb_ret_meta$sample_id, drop = FALSE]
  safe_fwrite(pb_ret_meta, file.path(opt$outdir, "23_validation_pseudobulk_manifest_retainedspace.tsv"))
  run_edger_case_control(pb_ret, pb_ret_meta,
                         file.path(opt$outdir, "24_validation_pseudobulk_retainedspace_ASD_vs_Control"),
                         case_label = opt$case_label, control_label = opt$control_label)
}

# candidate clusters
if (length(candidate_clusters) > 0) {
  val_meta_cand <- val_meta[cluster %in% candidate_clusters]
  if (nrow(val_meta_cand) > 0) {
    pb_cand <- aggregate_counts_sparse(val_counts[, val_meta_cand$cell, drop = FALSE], val_meta_cand$donor)
    pb_cand_meta <- unique(val_meta_cand[, .(sample_id = donor, donor, dx)])
    pb_cand_meta <- pb_cand_meta[match(colnames(pb_cand), sample_id)]
    pb_cand_meta[, n_cells := val_meta_cand[, .N, by = donor][match(sample_id, donor), N]]
    pb_cand_meta <- pb_cand_meta[n_cells >= opt$min_cells_cluster]
    pb_cand <- pb_cand[, pb_cand_meta$sample_id, drop = FALSE]
    safe_fwrite(pb_cand_meta, file.path(opt$outdir, "25_validation_pseudobulk_manifest_candidateclusters.tsv"))
    run_edger_case_control(pb_cand, pb_cand_meta,
                           file.path(opt$outdir, "26_validation_pseudobulk_candidateclusters_ASD_vs_Control"),
                           case_label = opt$case_label, control_label = opt$control_label)
  }
}

# reference clusters
if (length(reference_clusters) > 0) {
  val_meta_ref <- val_meta[cluster %in% reference_clusters]
  if (nrow(val_meta_ref) > 0) {
    pb_ref <- aggregate_counts_sparse(val_counts[, val_meta_ref$cell, drop = FALSE], val_meta_ref$donor)
    pb_ref_meta <- unique(val_meta_ref[, .(sample_id = donor, donor, dx)])
    pb_ref_meta <- pb_ref_meta[match(colnames(pb_ref), sample_id)]
    pb_ref_meta[, n_cells := val_meta_ref[, .N, by = donor][match(sample_id, donor), N]]
    pb_ref_meta <- pb_ref_meta[n_cells >= opt$min_cells_cluster]
    pb_ref <- pb_ref[, pb_ref_meta$sample_id, drop = FALSE]
    safe_fwrite(pb_ref_meta, file.path(opt$outdir, "27_validation_pseudobulk_manifest_referenceclusters.tsv"))
    run_edger_case_control(pb_ref, pb_ref_meta,
                           file.path(opt$outdir, "28_validation_pseudobulk_referenceclusters_ASD_vs_Control"),
                           case_label = opt$case_label, control_label = opt$control_label)
  }
}

# donor-delta candidate minus reference
if (length(candidate_clusters) > 0 && length(reference_clusters) > 0) {
  val_meta_pair <- val_meta[cluster %in% c(candidate_clusters, reference_clusters)]
  val_meta_pair[, group := ifelse(cluster %in% candidate_clusters, "candidate", "reference")]

  val_pair_group <- unique(val_meta_pair[, .(sample_id = paste(donor, group, sep = "__"), donor, dx, group)])
  val_pair_n <- val_meta_pair[, .(n_cells = .N), by = .(sample_id = paste(donor, group, sep = "__"))]
  val_pair_group <- merge(val_pair_group, val_pair_n, by = "sample_id", all.x = TRUE)
  val_pair_group <- val_pair_group[n_cells >= opt$min_cells_cluster][order(donor, group)]

  pb_pair <- aggregate_counts_sparse(
    val_counts[, val_meta_pair$cell, drop = FALSE],
    paste(val_meta_pair$donor, val_meta_pair$group, sep = "__")
  )
  pb_pair <- pb_pair[, val_pair_group$sample_id, drop = FALSE]

  safe_fwrite(val_pair_group, file.path(opt$outdir, "29_validation_pseudobulk_manifest_candidate_reference_pairs.tsv"))

  run_pseudobulk_donor_delta(
    pb_pair, val_pair_group,
    file.path(opt$outdir, "30_validation_pseudobulk_donorDelta_candidate_minus_reference_ASD_vs_Control"),
    group0_label = "reference",
    group2_label = "candidate",
    case_label = opt$case_label,
    control_label = opt$control_label
  )
}

# ============================================================
# 5. Optional pathway-level donor scores
# ============================================================
pathway_summary <- data.table()
if (!identical(opt$gmt, "NONE") && file.exists(opt$gmt)) {
  log_msg("Running optional pathway-level donor scoring using GMT: %s", opt$gmt)
  gmt_list <- parse_gmt(opt$gmt, opt$min_pathway_size, opt$max_pathway_size)
  disc_path <- select_discovery_pathways(
    disc_state_de,
    gmt_list,
    top_pos = opt$top_pathways_pos,
    top_neg = opt$top_pathways_neg
  )
  safe_fwrite(disc_path$table, file.path(opt$outdir, "31_discovery_selected_pathways.tsv"))

  if (nrow(donor_retained) > 0) {
    val_meta_retained <- val_meta[cluster %in% c(candidate_clusters, reference_clusters)]
    pb_ret <- aggregate_counts_sparse(val_counts[, val_meta_retained$cell, drop = FALSE], val_meta_retained$donor)
    pb_ret_meta <- unique(val_meta_retained[, .(sample_id = donor, donor, dx)])
    pb_ret_meta <- pb_ret_meta[match(colnames(pb_ret), sample_id)]

    y_ret <- edgeR::DGEList(counts = pb_ret)
    y_ret <- edgeR::calcNormFactors(y_ret)
    logcpm_ret <- edgeR::cpm(y_ret, log = TRUE, prior.count = 2, normalized.lib.sizes = TRUE)

    pw_pos <- score_pathways_on_pseudobulk(logcpm_ret, disc_path$pos, gmt_list)
    pw_neg <- score_pathways_on_pseudobulk(logcpm_ret, disc_path$neg, gmt_list)
    pw_all <- rbindlist(list(pw_pos, pw_neg), fill = TRUE)

    if (nrow(pw_all) > 0) {
      pw_all <- merge(pw_all, pb_ret_meta[, .(sample_id, donor, dx)], by = "sample_id", all.x = TRUE)
      safe_fwrite(pw_all, file.path(opt$outdir, "32_validation_retainedspace_pathway_scores.tsv.gz"))

      pathway_summary <- rbindlist(lapply(unique(pw_all$pathway), function(pw) {
        dsub <- pw_all[pathway == pw]
        w <- safe_wilcox(dsub, "pathway_score", case_label = opt$case_label, control_label = opt$control_label)
        l <- safe_lm_dx(dsub, "pathway_score", case_label = opt$case_label, control_label = opt$control_label)
        dir_label <- ifelse(pw %in% disc_path$pos, "disc_cluster2_pathway", "disc_cluster0_pathway")
        data.table(
          pathway = pw,
          direction = dir_label,
          mean_case = mean(dsub[dx == opt$case_label, pathway_score], na.rm = TRUE),
          mean_control = mean(dsub[dx == opt$control_label, pathway_score], na.rm = TRUE),
          effect_case_minus_control = l$effect,
          wilcox_p = w$p,
          lm_p = l$p
        )
      }), fill = TRUE)
      pathway_summary[, p_adj_bh := p.adjust(wilcox_p, method = "BH")]
      safe_fwrite(pathway_summary, file.path(opt$outdir, "33_validation_pathway_replication_stats.tsv"))
    }
  }
} else {
  log_msg("Skipping pathway-level donor scoring: no GMT provided")
}

# ============================================================
# 6. Replication decision
# ============================================================
log_msg("Making replication decision")

mapping_support <- FALSE
if (nrow(mapping_dt) > 0) {
  top_cand <- mapping_dt[candidate_like == TRUE][order(-corr_delta_c2_minus_c0, -candidate_minus_reference)]
  top_ref  <- mapping_dt[reference_like == TRUE][order(corr_delta_c2_minus_c0, candidate_minus_reference)]
  mapping_support <- nrow(top_cand) > 0 && nrow(top_ref) > 0 &&
    is.finite(top_cand$corr_delta_c2_minus_c0[1]) &&
    top_cand$corr_delta_c2_minus_c0[1] >= opt$mapping_corr_delta_threshold
}

donor_sig <- rep_stats[
  metric %in% c(
    "candidateclusters_candidate_minus_reference",
    "candidate_vs_reference_shift",
    "candidate_cluster_proportion_within_retainedspace",
    "log2_ratio_candidate_over_reference"
  )
]
donor_sig[, supported_direction := effect_case_minus_control > 0]
donor_sig[, significant := is.finite(wilcox_p) & wilcox_p < opt$alpha]
n_donor_sig <- sum(donor_sig$significant & donor_sig$supported_direction, na.rm = TRUE)

pathway_support <- FALSE
if (nrow(pathway_summary) > 0) {
  pw_pos <- pathway_summary[direction == "disc_cluster2_pathway"]
  pathway_support <- any(is.finite(pw_pos$wilcox_p) & pw_pos$wilcox_p < opt$alpha & pw_pos$effect_case_minus_control > 0, na.rm = TRUE)
}

decision <- if (mapping_support && n_donor_sig >= 2 && pathway_support) {
  "supported cross-cohort replication"
} else if (mapping_support && (n_donor_sig >= 1 || pathway_support)) {
  "partial cross-cohort replication"
} else {
  "limited cross-cohort replication"
}

decision_dt <- data.table(
  mapping_support = mapping_support,
  n_candidate_clusters = length(candidate_clusters),
  n_reference_clusters = length(reference_clusters),
  n_donor_supported_metrics = n_donor_sig,
  pathway_support = pathway_support,
  replication_decision = decision
)
safe_fwrite(decision_dt, file.path(opt$outdir, "34_replication_decision.tsv"))

decision_lines <- c(
  sprintf("Replication decision: %s", decision),
  sprintf("Mapping support: %s", mapping_support),
  sprintf("Candidate-like validation clusters: %s", paste(candidate_clusters, collapse = ", ")),
  sprintf("Reference-like validation clusters: %s", paste(reference_clusters, collapse = ", ")),
  sprintf("Number of donor-level supported metrics: %s", n_donor_sig),
  sprintf("Pathway support: %s", pathway_support),
  "Stopping rule:",
  "- If cluster mapping is weak, donor scores are not supportive, and pathway direction is inconsistent -> limited cross-cohort replication.",
  "- If mapping is present and at least one donor-aware or pathway-level feature is directionally supportive -> partial cross-cohort replication.",
  "- If mapping is present, donor-aware support is stronger, and pathway-level support is concordant -> supported cross-cohort replication."
)
writeLines(decision_lines, con = file.path(opt$outdir, "35_replication_decision_notes.txt"))

capture.output(sessionInfo(), file = file.path(opt$outdir, "36_sessionInfo.txt"))
log_msg("Package 4 completed. Output directory: %s", opt$outdir)
