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

option_list <- list(
  make_option("--rds", type = "character",
              default = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Package2_strict_input_fixed_from_original.rds"),
  make_option("--outdir", type = "character",
              default = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Package3_Discovery_DonorAwareDiseaseAssociation"),
  make_option("--assay", type = "character", default = "RNA"),
  make_option("--cluster_col", type = "character", default = "AUTO"),
  make_option("--donor_col", type = "character", default = "AUTO"),
  make_option("--dx_col", type = "character", default = "AUTO"),
  make_option("--case_label", type = "character", default = "ASD"),
  make_option("--control_label", type = "character", default = "Control"),
  make_option("--min_cells_donor", type = "integer", default = 20),
  make_option("--min_cells_cluster", type = "integer", default = 5),
  make_option("--oof_top_n_pos", type = "integer", default = 100),
  make_option("--oof_top_n_neg", type = "integer", default = 100),
  make_option("--threads", type = "integer", default = 8)
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

safe_message <- function(...) cat(sprintf(...), "\n")

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

get_assay_mat <- function(obj, assay = "RNA", layer_name = c("counts", "data")) {
  layer_name <- match.arg(layer_name)
  mat <- tryCatch({
    GetAssayData(obj, assay = assay, slot = layer_name)
  }, error = function(e) {
    GetAssayData(obj, assay = assay, layer = layer_name)
  })
  mat
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

safe_quasibinomial <- function(df, success_col, failure_col, dx_col = "dx", case_label = "ASD", control_label = "Control") {
  d <- copy(df)
  d <- d[get(dx_col) %in% c(control_label, case_label)]
  if (nrow(d) < 3) return(list(logit_effect = NA_real_, p = NA_real_))
  d[, dx_factor := factor(get(dx_col), levels = c(control_label, case_label))]
  fit <- tryCatch(
    glm(cbind(get(success_col), get(failure_col)) ~ dx_factor, family = quasibinomial(), data = d),
    error = function(e) NULL
  )
  if (is.null(fit)) return(list(logit_effect = NA_real_, p = NA_real_))
  sm <- summary(fit)$coefficients
  coef_name <- grep("^dx_factor", rownames(sm), value = TRUE)
  if (length(coef_name) == 0) return(list(logit_effect = NA_real_, p = NA_real_))
  list(logit_effect = unname(sm[coef_name[1], "Estimate"]), p = unname(sm[coef_name[1], "Pr(>|t|)"]))
}

plot_metric_box <- function(df, metric_col, ylab_text, out_pdf, title_text = NULL, dx_col = "dx", case_label = "ASD", control_label = "Control") {
  d <- copy(df)
  d <- d[get(dx_col) %in% c(control_label, case_label) & is.finite(get(metric_col))]
  if (nrow(d) == 0) return(invisible(NULL))
  d[, dx_factor := factor(get(dx_col), levels = c(control_label, case_label))]
  p <- ggplot(d, aes(x = dx_factor, y = get(metric_col))) +
    geom_violin(trim = FALSE, scale = "width", alpha = 0.25) +
    geom_boxplot(width = 0.18, outlier.shape = NA) +
    geom_jitter(width = 0.12, height = 0, alpha = 0.75, size = 1.8) +
    theme_bw(base_size = 12) +
    labs(x = NULL, y = ylab_text, title = title_text)
  ggsave(out_pdf, p, width = 4.8, height = 4.2)
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

run_edger_case_control <- function(pb_counts, sample_meta, out_prefix, dx_col = "dx", case_label = "ASD", control_label = "Control") {
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

run_paired_cluster_state_de <- function(pb_counts, sample_meta, out_prefix) {
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
  meta[, cluster := factor(as.character(cluster), levels = c("0", "2"))]
  meta[, donor := as.character(donor)]

  donor_keep <- meta[, .N, by = .(donor, cluster)][, dcast(.SD, donor ~ cluster, value.var = "N", fill = 0)]
  donor_keep <- donor_keep[`0` > 0 & `2` > 0, donor]
  meta <- meta[donor %in% donor_keep]

  if (nrow(meta) < 4 || uniqueN(meta$donor) < 2) {
    safe_fwrite(res_empty, paste0(out_prefix, ".tsv.gz"))
    return(res_empty)
  }

  pb_counts <- pb_counts[, meta$sample_id, drop = FALSE]

  donor_factor <- factor(meta$donor)
  cluster_factor <- factor(meta$cluster, levels = c("0", "2"))
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

extract_top_state_genes <- function(state_de, n_up = 100, n_down = 100) {
  d <- as.data.table(state_de)
  d <- d[is.finite(logFC) & is.finite(PValue)]
  if (nrow(d) == 0) {
    return(data.table(
      direction = character(),
      rank = integer(),
      gene = character(),
      logFC = numeric(),
      PValue = numeric(),
      FDR = numeric()
    ))
  }

  up <- d[logFC > 0][order(PValue, -abs(logFC))]
  down <- d[logFC < 0][order(PValue, -abs(logFC))]

  up <- up[seq_len(min(nrow(up), n_up))]
  down <- down[seq_len(min(nrow(down), n_down))]

  if (nrow(up) > 0) {
    up[, direction := "Cluster2_up"]
    up[, rank := seq_len(.N)]
  }
  if (nrow(down) > 0) {
    down[, direction := "Cluster0_up"]
    down[, rank := seq_len(.N)]
  }

  out <- rbindlist(list(
    up[, .(direction, rank, gene, logFC, PValue, FDR)],
    down[, .(direction, rank, gene, logFC, PValue, FDR)]
  ), fill = TRUE)

  out
}

compute_oof_state_concordance <- function(pb_counts, sample_meta, out_prefix,
                                          top_n_pos = 100, top_n_neg = 100,
                                          case_label = "ASD", control_label = "Control") {
  res_empty <- data.table(
    donor = character(),
    dx = character(),
    sample0 = character(),
    sample2 = character(),
    n_pos = integer(),
    n_neg = integer(),
    n_genes_used = integer(),
    state_concordance = numeric()
  )

  meta <- as.data.table(sample_meta)
  if (ncol(pb_counts) != nrow(meta)) stop("OOF: pseudobulk counts 列数与 sample_meta 行数不一致")

  meta <- copy(meta)
  meta[, cluster := as.character(cluster)]
  meta[, dx := as.character(dx)]

  pair_map <- dcast(
    meta[, .(donor, dx, cluster, sample_id)],
    donor + dx ~ cluster,
    value.var = "sample_id"
  )

  if (!all(c("0", "2") %in% colnames(pair_map))) {
    safe_fwrite(res_empty, paste0(out_prefix, ".tsv"))
    return(res_empty)
  }

  pair_map <- pair_map[!is.na(`0`) & !is.na(`2`)]
  pair_map <- pair_map[dx %in% c(control_label, case_label)]

  if (nrow(pair_map) < 4) {
    safe_fwrite(res_empty, paste0(out_prefix, ".tsv"))
    return(res_empty)
  }

  y_all <- edgeR::DGEList(counts = pb_counts)
  y_all <- edgeR::calcNormFactors(y_all)
  logcpm_all <- edgeR::cpm(y_all, log = TRUE, prior.count = 2, normalized.lib.sizes = TRUE)

  donor_results <- vector("list", nrow(pair_map))

  for (i in seq_len(nrow(pair_map))) {
    donor_holdout <- pair_map$donor[i]

    train_meta <- meta[donor != donor_holdout]
    train_pair_map <- dcast(
      train_meta[, .(donor, dx, cluster, sample_id)],
      donor + dx ~ cluster,
      value.var = "sample_id"
    )

    if (!all(c("0", "2") %in% colnames(train_pair_map))) {
      donor_results[[i]] <- NULL
      next
    }

    train_pair_map <- train_pair_map[!is.na(`0`) & !is.na(`2`)]
    if (nrow(train_pair_map) < 4) {
      donor_results[[i]] <- NULL
      next
    }

    train_donors <- unique(train_pair_map$donor)
    train_meta <- train_meta[donor %in% train_donors]
    train_meta <- unique(train_meta[, .(sample_id, donor, dx, cluster)])
    train_meta <- train_meta[order(donor, cluster)]

    if (!all(train_meta$sample_id %in% colnames(pb_counts))) {
      donor_results[[i]] <- NULL
      next
    }

    train_counts <- pb_counts[, train_meta$sample_id, drop = FALSE]
    train_de <- run_paired_cluster_state_de(
      train_counts,
      train_meta,
      out_prefix = tempfile(pattern = "tmp_oof_state_", tmpdir = tempdir())
    )

    train_de <- as.data.table(train_de)
    train_de <- train_de[is.finite(logFC) & is.finite(PValue)]

    pos_genes <- train_de[logFC > 0][order(PValue, -abs(logFC)), head(gene, top_n_pos)]
    neg_genes <- train_de[logFC < 0][order(PValue, -abs(logFC)), head(gene, top_n_neg)]
    use_genes <- unique(c(pos_genes, neg_genes))

    if (length(use_genes) < 10) {
      donor_results[[i]] <- data.table(
        donor = donor_holdout,
        dx = pair_map$dx[i],
        sample0 = pair_map$`0`[i],
        sample2 = pair_map$`2`[i],
        n_pos = length(pos_genes),
        n_neg = length(neg_genes),
        n_genes_used = length(use_genes),
        state_concordance = NA_real_
      )
      next
    }

    template_dt <- train_de[match(use_genes, gene)]
    template_vec <- template_dt$logFC
    names(template_vec) <- template_dt$gene

    s0 <- pair_map$`0`[i]
    s2 <- pair_map$`2`[i]

    if (!(s0 %in% colnames(logcpm_all) && s2 %in% colnames(logcpm_all))) {
      donor_results[[i]] <- data.table(
        donor = donor_holdout,
        dx = pair_map$dx[i],
        sample0 = s0,
        sample2 = s2,
        n_pos = length(pos_genes),
        n_neg = length(neg_genes),
        n_genes_used = length(use_genes),
        state_concordance = NA_real_
      )
      next
    }

    delta_vec <- logcpm_all[use_genes, s2, drop = TRUE] - logcpm_all[use_genes, s0, drop = TRUE]
    delta_vec <- delta_vec[names(template_vec)]
    keep <- is.finite(delta_vec) & is.finite(template_vec)

    score <- if (sum(keep) >= 10) {
      suppressWarnings(cor(delta_vec[keep], template_vec[keep], method = "spearman"))
    } else {
      NA_real_
    }

    donor_results[[i]] <- data.table(
      donor = donor_holdout,
      dx = pair_map$dx[i],
      sample0 = s0,
      sample2 = s2,
      n_pos = length(pos_genes),
      n_neg = length(neg_genes),
      n_genes_used = sum(keep),
      state_concordance = score
    )
  }

  out <- rbindlist(donor_results, fill = TRUE)
  safe_fwrite(out, paste0(out_prefix, ".tsv"))
  out
}

run_pseudobulk_donor_delta <- function(pb_counts, sample_meta, out_prefix,
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
  if (ncol(pb_counts) != nrow(meta)) {
    stop("paired donor-delta: pseudobulk counts 列数与 sample_meta 行数不一致")
  }

  meta <- copy(meta)
  meta[, cluster := as.character(cluster)]
  meta[, dx := as.character(dx)]

  meta <- meta[cluster %in% c("0", "2") & dx %in% c(control_label, case_label)]
  if (nrow(meta) == 0) {
    safe_fwrite(res_empty, paste0(out_prefix, ".tsv.gz"))
    return(res_empty)
  }

  pair_map <- dcast(
    meta[, .(donor, dx, cluster, sample_id)],
    donor + dx ~ cluster,
    value.var = "sample_id"
  )

  if (!all(c("0", "2") %in% colnames(pair_map))) {
    safe_fwrite(res_empty, paste0(out_prefix, ".tsv.gz"))
    return(res_empty)
  }

  pair_map <- pair_map[!is.na(`0`) & !is.na(`2`)]

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

  if (!all(c(pair_map$`0`, pair_map$`2`) %in% colnames(pb_counts))) {
    stop("pair_map 中的 sample_id 有部分不在 pseudobulk counts 列名中")
  }

  y <- edgeR::DGEList(counts = pb_counts)
  y <- edgeR::calcNormFactors(y)
  logcpm <- edgeR::cpm(y, log = TRUE, prior.count = 2, normalized.lib.sizes = TRUE)

  idx0 <- match(pair_map$`0`, colnames(logcpm))
  idx2 <- match(pair_map$`2`, colnames(logcpm))

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

# -------------------------------------------------------------------
# 0. 读取对象与 metadata 解析
# -------------------------------------------------------------------
set.seed(1)
try(data.table::setDTthreads(opt$threads), silent = TRUE)

obj <- readRDS(opt$rds)
assay_use <- if (opt$assay %in% Assays(obj)) opt$assay else DefaultAssay(obj)
DefaultAssay(obj) <- assay_use
meta <- as.data.table(obj@meta.data, keep.rownames = "cell")

cluster_col <- resolve_column(meta, opt$cluster_col, "cluster")
donor_col   <- resolve_column(meta, opt$donor_col,   "donor")
dx_col      <- resolve_column(meta, opt$dx_col,      "dx")

meta[, cluster_raw := as.character(get(cluster_col))]
meta[, donor := as.character(get(donor_col))]
meta[, dx_raw := as.character(get(dx_col))]
meta[, dx := normalize_dx(dx_raw, case_label = opt$case_label, control_label = opt$control_label)]
meta[, cluster := cluster_raw]

cluster_keep <- c("0", "2")
meta <- meta[cluster %in% cluster_keep]
meta <- meta[!is.na(donor) & donor != "" & !is.na(dx)]

if (nrow(meta) == 0) stop("过滤后没有可用细胞。请检查 cluster / donor / dx 列解析是否正确。")

obj <- subset(obj, cells = meta$cell)
meta <- as.data.table(obj@meta.data, keep.rownames = "cell")
meta[, cluster_raw := as.character(get(cluster_col))]
meta[, donor := as.character(get(donor_col))]
meta[, dx_raw := as.character(get(dx_col))]
meta[, dx := normalize_dx(dx_raw, case_label = opt$case_label, control_label = opt$control_label)]
meta[, cluster := cluster_raw]
meta <- meta[cluster %in% cluster_keep & !is.na(donor) & donor != "" & !is.na(dx)]
setkey(meta, cell)

safe_fwrite(
  data.table(
    input_rds = opt$rds,
    assay = assay_use,
    cluster_col = cluster_col,
    donor_col = donor_col,
    dx_col = dx_col,
    n_cells = nrow(meta),
    n_donors = uniqueN(meta$donor),
    n_asd = uniqueN(meta[dx == opt$case_label, donor]),
    n_control = uniqueN(meta[dx == opt$control_label, donor]),
    min_cells_donor = opt$min_cells_donor,
    min_cells_cluster = opt$min_cells_cluster,
    oof_top_n_pos = opt$oof_top_n_pos,
    oof_top_n_neg = opt$oof_top_n_neg
  ),
  file.path(opt$outdir, "00_run_parameters.tsv")
)

safe_fwrite(
  meta[, .(cell, donor, dx, cluster)],
  file.path(opt$outdir, "01_cell_manifest_used.tsv.gz")
)

# -------------------------------------------------------------------
# 1. donor-level abundance / state composition
# -------------------------------------------------------------------
abund <- meta[, .(
  n_total_retained = .N,
  n_c0 = sum(cluster == "0"),
  n_c2 = sum(cluster == "2")
), by = .(donor, dx)]

abund <- abund[n_total_retained >= opt$min_cells_donor]
abund[, prop_c0 := n_c0 / n_total_retained]
abund[, prop_c2 := n_c2 / n_total_retained]
abund[, log2_ratio_c2_over_c0 := log2((n_c2 + 0.5) / (n_c0 + 0.5))]
abund[, odds_c2_over_c0 := (n_c2 + 0.5) / (n_c0 + 0.5)]

safe_fwrite(abund, file.path(opt$outdir, "02_abundance_per_donor.tsv"))

abund_stats <- rbindlist(list(
  {
    w <- safe_wilcox(abund, "prop_c2", case_label = opt$case_label, control_label = opt$control_label)
    l <- safe_lm_dx(abund, "prop_c2", dx_col = "dx", case_label = opt$case_label, control_label = opt$control_label)
    q <- safe_quasibinomial(abund, "n_c2", "n_c0", dx_col = "dx", case_label = opt$case_label, control_label = opt$control_label)
    data.table(
      metric = "prop_c2_within_retained_core",
      mean_case = mean(abund[dx == opt$case_label, prop_c2], na.rm = TRUE),
      mean_control = mean(abund[dx == opt$control_label, prop_c2], na.rm = TRUE),
      effect_case_minus_control = l$effect,
      wilcox_p = w$p,
      lm_p = l$p,
      quasibinomial_logit_effect = q$logit_effect,
      quasibinomial_p = q$p
    )
  },
  {
    w <- safe_wilcox(abund, "log2_ratio_c2_over_c0", case_label = opt$case_label, control_label = opt$control_label)
    l <- safe_lm_dx(abund, "log2_ratio_c2_over_c0", dx_col = "dx", case_label = opt$case_label, control_label = opt$control_label)
    data.table(
      metric = "log2_ratio_c2_over_c0",
      mean_case = mean(abund[dx == opt$case_label, log2_ratio_c2_over_c0], na.rm = TRUE),
      mean_control = mean(abund[dx == opt$control_label, log2_ratio_c2_over_c0], na.rm = TRUE),
      effect_case_minus_control = l$effect,
      wilcox_p = w$p,
      lm_p = l$p,
      quasibinomial_logit_effect = NA_real_,
      quasibinomial_p = NA_real_
    )
  }
), fill = TRUE)

safe_fwrite(abund_stats, file.path(opt$outdir, "03_abundance_stats.tsv"))
plot_metric_box(abund, "prop_c2", "Cluster 2 proportion within retained core",
                file.path(opt$outdir, "04_plot_prop_c2.pdf"),
                "Discovery donor-level abundance")
plot_metric_box(abund, "log2_ratio_c2_over_c0", "log2((C2+0.5)/(C0+0.5))",
                file.path(opt$outdir, "05_plot_log2ratio_c2_over_c0.pdf"),
                "Discovery donor-level state ratio")

# -------------------------------------------------------------------
# 2. donor-cluster pseudobulk manifest（用于 paired state / donor-delta）
# -------------------------------------------------------------------
mat_counts <- get_assay_mat(obj, assay = assay_use, layer_name = "counts")
mat_counts <- mat_counts[, meta$cell, drop = FALSE]

paired_meta <- meta[, .N, by = .(donor, dx, cluster)]
paired_meta <- paired_meta[N >= opt$min_cells_cluster]

paired_wide <- dcast(
  paired_meta,
  donor + dx ~ cluster,
  value.var = "N",
  fill = 0
)

paired_donors <- paired_wide[`0` > 0 & `2` > 0, donor]
meta_pair <- meta[donor %in% paired_donors]

meta_pair_group <- unique(
  meta_pair[, .(sample_id = paste(donor, cluster, sep = "__"), donor, dx, cluster)]
)

meta_pair_n <- meta_pair[, .(n_cells = .N), by = .(sample_id = paste(donor, cluster, sep = "__"))]
meta_pair_group <- merge(meta_pair_group, meta_pair_n, by = "sample_id", all.x = TRUE)
meta_pair_group <- meta_pair_group[n_cells >= opt$min_cells_cluster]
meta_pair_group <- meta_pair_group[order(donor, cluster)]

counts_pair <- aggregate_counts_sparse(
  mat_counts[, meta_pair$cell, drop = FALSE],
  paste(meta_pair$donor, meta_pair$cluster, sep = "__")
)

counts_pair <- counts_pair[, meta_pair_group$sample_id, drop = FALSE]

safe_fwrite(
  meta_pair_group,
  file.path(opt$outdir, "06_paired_pseudobulk_manifest_donor_cluster.tsv")
)

# -------------------------------------------------------------------
# 3. data-driven state definition: paired cluster2 vs cluster0 pseudobulk DE
# -------------------------------------------------------------------
state_de <- run_paired_cluster_state_de(
  counts_pair,
  meta_pair_group,
  file.path(opt$outdir, "07_state_definition_cluster2_vs_cluster0_paired")
)

top_state_genes <- extract_top_state_genes(state_de, n_up = 100, n_down = 100)
safe_fwrite(top_state_genes, file.path(opt$outdir, "08_state_definition_top_genes.tsv"))

# -------------------------------------------------------------------
# 4. out-of-fold donor state concordance（data-driven, no hardcoded modules）
# -------------------------------------------------------------------
oof_scores <- compute_oof_state_concordance(
  counts_pair,
  meta_pair_group,
  file.path(opt$outdir, "09_oof_state_concordance_per_donor"),
  top_n_pos = opt$oof_top_n_pos,
  top_n_neg = opt$oof_top_n_neg,
  case_label = opt$case_label,
  control_label = opt$control_label
)

oof_stats <- {
  w <- safe_wilcox(oof_scores, "state_concordance",
                   case_label = opt$case_label,
                   control_label = opt$control_label)
  l <- safe_lm_dx(oof_scores, "state_concordance",
                  case_label = opt$case_label,
                  control_label = opt$control_label)
  data.table(
    metric = "oof_state_concordance_cluster2_minus_cluster0",
    mean_case = mean(oof_scores[dx == opt$case_label, state_concordance], na.rm = TRUE),
    mean_control = mean(oof_scores[dx == opt$control_label, state_concordance], na.rm = TRUE),
    effect_case_minus_control = l$effect,
    wilcox_p = w$p,
    lm_p = l$p,
    n_case_donor = uniqueN(oof_scores[dx == opt$case_label, donor]),
    n_control_donor = uniqueN(oof_scores[dx == opt$control_label, donor])
  )
}
safe_fwrite(oof_stats, file.path(opt$outdir, "10_oof_state_concordance_stats.tsv"))

plot_metric_box(
  oof_scores,
  "state_concordance",
  "Out-of-fold state concordance",
  file.path(opt$outdir, "11_plot_oof_state_concordance.pdf"),
  "Discovery donor-level state concordance"
)

# -------------------------------------------------------------------
# 5. pooled retained pseudobulk ASD vs Control
# -------------------------------------------------------------------
pb_pooled_groups <- meta$donor
pb_pooled_counts <- aggregate_counts_sparse(mat_counts, pb_pooled_groups)
pb_pooled_meta <- unique(meta[, .(sample_id = donor, donor, dx)])
pb_pooled_meta <- pb_pooled_meta[match(colnames(pb_pooled_counts), sample_id)]
pb_pooled_meta[, n_cells := meta[, .N, by = donor][match(sample_id, donor), N]]
pb_pooled_meta <- pb_pooled_meta[n_cells >= opt$min_cells_donor]
pb_pooled_counts <- pb_pooled_counts[, pb_pooled_meta$sample_id, drop = FALSE]
safe_fwrite(pb_pooled_meta, file.path(opt$outdir, "12_pseudobulk_manifest_pooled.tsv"))
run_edger_case_control(
  pb_pooled_counts,
  pb_pooled_meta,
  file.path(opt$outdir, "13_pseudobulk_all_retained_ASD_vs_Control"),
  dx_col = "dx",
  case_label = opt$case_label,
  control_label = opt$control_label
)

# -------------------------------------------------------------------
# 6. cluster-specific pseudobulk ASD vs Control
# -------------------------------------------------------------------
for (cl in c("0", "2")) {
  meta_cl <- meta[cluster == cl]
  if (nrow(meta_cl) == 0) next

  grp <- meta_cl$donor
  pb_counts_cl <- aggregate_counts_sparse(mat_counts[, meta_cl$cell, drop = FALSE], grp)
  pb_meta_cl <- unique(meta_cl[, .(sample_id = donor, donor, dx)])
  pb_meta_cl <- pb_meta_cl[match(colnames(pb_counts_cl), sample_id)]
  pb_meta_cl[, n_cells := meta_cl[, .N, by = donor][match(sample_id, donor), N]]
  pb_meta_cl <- pb_meta_cl[n_cells >= opt$min_cells_cluster]
  pb_counts_cl <- pb_counts_cl[, pb_meta_cl$sample_id, drop = FALSE]
  pb_meta_cl[, cluster := cl]

  if (cl == "0") {
    safe_fwrite(pb_meta_cl, file.path(opt$outdir, "14_pseudobulk_manifest_cluster0.tsv"))
    run_edger_case_control(
      pb_counts_cl,
      pb_meta_cl,
      file.path(opt$outdir, "15_pseudobulk_cluster0_ASD_vs_Control"),
      dx_col = "dx",
      case_label = opt$case_label,
      control_label = opt$control_label
    )
  } else if (cl == "2") {
    safe_fwrite(pb_meta_cl, file.path(opt$outdir, "16_pseudobulk_manifest_cluster2.tsv"))
    run_edger_case_control(
      pb_counts_cl,
      pb_meta_cl,
      file.path(opt$outdir, "17_pseudobulk_cluster2_ASD_vs_Control"),
      dx_col = "dx",
      case_label = opt$case_label,
      control_label = opt$control_label
    )
  }
}

# -------------------------------------------------------------------
# 7. donor-delta pseudobulk ASD vs Control
# -------------------------------------------------------------------
run_pseudobulk_donor_delta(
  counts_pair,
  meta_pair_group,
  file.path(opt$outdir, "18_pseudobulk_donorDelta_cluster2_minus_cluster0_ASD_vs_Control"),
  case_label = opt$case_label,
  control_label = opt$control_label
)

# -------------------------------------------------------------------
# 8. analysis notes
# -------------------------------------------------------------------
summary_lines <- c(
  "Package 3 = discovery cohort donor-aware disease association (re-written to minimize heuristic hardcoded modules).",
  "Abundance main endpoint: donor-level Cluster2 proportion within retained core microglia, plus Cluster2:Cluster0 ratio.",
  "State-definition main endpoint: paired donor-aware pseudobulk DE for Cluster2 vs Cluster0 across all paired donors (dx-agnostic).",
  "Program/state-shift main endpoint: out-of-fold donor state concordance score computed from training-donor state template and held-out donor [Cluster2 - Cluster0] delta.",
  "Disease-expression endpoints: pooled retained pseudobulk ASD vs Control, cluster0-specific ASD vs Control, cluster2-specific ASD vs Control, and donor-level pseudobulk delta model [Cluster2 - Cluster0] ~ dx.",
  "No hand-curated SPP1-axis or activation/homeostatic heuristic modules are used as primary disease-defining signatures in this version."
)
writeLines(summary_lines, con = file.path(opt$outdir, "19_analysis_notes.txt"))

capture.output(sessionInfo(), file = file.path(opt$outdir, "20_sessionInfo.txt"))
safe_message("Package 3 完成。输出目录：%s", opt$outdir)
