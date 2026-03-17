#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
  stop("Usage: Rscript Step08_Gandal2022_bulk_program_validation_v3.R <RData> <candidate_genes.txt> <reference_genes.txt> <outdir>")
}

rdata_path     <- args[1]
candidate_path <- args[2]
reference_path <- args[3]
outdir         <- args[4]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(nlme)
  library(ggplot2)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

logfile <- file.path(outdir, "00_step08_run.log")
log_msg <- function(...) {
  msg <- sprintf(...)
  line <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg)
  cat(line, "\n")
  cat(line, "\n", file = logfile, append = TRUE)
}

safe_fwrite <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, file = file, sep = "\t", quote = FALSE, na = "NA")
}

read_gene_list <- function(path) {
  x <- readLines(path, warn = FALSE)
  x <- toupper(trimws(x))
  x <- x[x != "" & !is.na(x)]
  unique(x)
}

standardize_dx <- function(x) {
  x0 <- trimws(as.character(x))
  xl <- tolower(x0)
  out <- rep(NA_character_, length(x0))
  out[xl %in% c("asd", "autism", "case")] <- "ASD"
  out[xl %in% c("ctl", "control", "ctrl")] <- "Control"
  out[xl %in% c("dup15q")] <- "Dup15q"
  factor(out, levels = c("Control", "ASD", "Dup15q"))
}

usable_factor <- function(x) {
  x <- x[!is.na(x)]
  length(x) > 0 && length(unique(x)) > 1
}

usable_numeric <- function(x) {
  x <- x[is.finite(x)]
  length(x) > 0 && length(unique(x)) > 1
}

score_gene_set <- function(gset, mat_z) {
  hit <- intersect(gset, rownames(mat_z))
  if (length(hit) == 0) {
    stop("No overlap between supplied gene set and mapped expression matrix.")
  }
  list(
    score = colMeans(mat_z[hit, , drop = FALSE], na.rm = TRUE),
    hit = hit
  )
}

fit_one_program <- function(df, score_col) {
  # 候选协变量
  covar_candidates <- c("region", "SeqBatch", "Sex", "Ancestry", "Age", "Age_sqd", "PMI", "RIN")

  keep_covars <- c()
  for (cc in covar_candidates) {
    if (!cc %in% colnames(df)) next
    if (is.numeric(df[[cc]])) {
      if (usable_numeric(df[[cc]])) keep_covars <- c(keep_covars, cc)
    } else {
      if (usable_factor(df[[cc]])) keep_covars <- c(keep_covars, cc)
    }
  }

  # 分层 fallback，避免模型太脆
  covar_sets <- list(
    keep_covars,
    setdiff(keep_covars, "Age_sqd"),
    intersect(keep_covars, c("region", "SeqBatch", "Sex", "Age", "PMI", "RIN")),
    intersect(keep_covars, c("region", "Sex", "Age", "PMI", "RIN")),
    intersect(keep_covars, c("region", "Sex", "Age"))
  )

  covar_sets <- Filter(function(x) TRUE, covar_sets)

  for (covars_use in covar_sets) {
    need_cols <- unique(c(score_col, "dx_std", "Subject", covars_use))
    dat <- copy(df)[, ..need_cols]
    dat <- dat[complete.cases(dat)]

    if (nrow(dat) < 8) next
    if (length(unique(dat$dx_std)) < 2) next

    dat[, dx_std := factor(dx_std, levels = c("Control", "ASD"))]
    if ("Subject" %in% colnames(dat)) dat[, Subject := as.factor(Subject)]
    for (cc in covars_use) {
      if (!is.numeric(dat[[cc]])) dat[[cc]] <- as.factor(dat[[cc]])
    }

    rhs <- c("dx_std", covars_use)
    form <- as.formula(paste(score_col, "~", paste(rhs, collapse = " + ")))

    use_random_subject <- ("Subject" %in% colnames(dat)) && (length(unique(dat$Subject)) < nrow(dat))

    if (use_random_subject) {
      fit <- tryCatch(
        nlme::lme(
          fixed = form,
          random = ~1 | Subject,
          data = dat,
          method = "REML",
          na.action = na.omit,
          control = nlme::lmeControl(returnObject = TRUE, opt = "optim")
        ),
        error = function(e) NULL
      )

      if (!is.null(fit)) {
        tt <- summary(fit)$tTable
        term <- "dx_stdASD"
        if (term %in% rownames(tt)) {
          return(data.table(
            program = score_col,
            model_type = "lme_random_intercept_subject",
            covariates_used = paste(covars_use, collapse = ";"),
            n_samples = nrow(dat),
            n_subjects = uniqueN(dat$Subject),
            beta_ASD_vs_Control = unname(tt[term, "Value"]),
            se = unname(tt[term, "Std.Error"]),
            t_value = unname(tt[term, "t-value"]),
            p_value = unname(tt[term, "p-value"])
          ))
        }
      }
    }

    fit_lm <- tryCatch(
      lm(form, data = dat),
      error = function(e) NULL
    )
    if (!is.null(fit_lm)) {
      sm <- summary(fit_lm)$coefficients
      term <- "dx_stdASD"
      if (term %in% rownames(sm)) {
        return(data.table(
          program = score_col,
          model_type = "lm_fallback",
          covariates_used = paste(covars_use, collapse = ";"),
          n_samples = nrow(dat),
          n_subjects = if ("Subject" %in% colnames(dat)) uniqueN(dat$Subject) else NA_integer_,
          beta_ASD_vs_Control = unname(sm[term, "Estimate"]),
          se = unname(sm[term, "Std. Error"]),
          t_value = unname(sm[term, "t value"]),
          p_value = unname(sm[term, grep("^Pr", colnames(sm))])
        ))
      }
    }
  }

  data.table(
    program = score_col,
    model_type = NA_character_,
    covariates_used = NA_character_,
    n_samples = NA_integer_,
    n_subjects = NA_integer_,
    beta_ASD_vs_Control = NA_real_,
    se = NA_real_,
    t_value = NA_real_,
    p_value = NA_real_
  )
}

make_boxplot <- function(df, score_col, ylab_text, title_text, out_pdf, out_png) {
  d <- copy(df)[!is.na(dx_std) & is.finite(get(score_col))]
  if (nrow(d) == 0) return(invisible(NULL))
  d[, dx_std := factor(dx_std, levels = c("Control", "ASD"))]

  p <- ggplot(d, aes(x = dx_std, y = get(score_col))) +
    geom_violin(trim = FALSE, scale = "width", alpha = 0.25) +
    geom_boxplot(width = 0.18, outlier.shape = NA) +
    geom_jitter(width = 0.12, height = 0, alpha = 0.75, size = 1.4) +
    theme_bw(base_size = 12) +
    labs(x = NULL, y = ylab_text, title = title_text)

  ggsave(out_pdf, p, width = 4.8, height = 4.2)
  ggsave(out_png, p, width = 4.8, height = 4.2, dpi = 300)
}

log_msg("Loading RData: %s", rdata_path)
e <- new.env()
load(rdata_path, envir = e)

req_objs <- c("datExpr", "datMeta", "datMeta_model")
if (!all(req_objs %in% ls(e))) {
  stop(sprintf("RData missing required objects: %s", paste(setdiff(req_objs, ls(e)), collapse = ", ")))
}

expr <- e$datExpr
meta <- as.data.table(e$datMeta)
meta_model <- as.data.table(e$datMeta_model)

# -----------------------------
# 1. 对齐
# -----------------------------
samples <- colnames(expr)
if (is.null(samples)) stop("datExpr has no colnames.")
if (!("sample_id" %in% colnames(meta))) stop("datMeta lacks sample_id.")

ord_meta <- match(samples, meta$sample_id)
if (any(is.na(ord_meta))) stop("Failed to align datMeta to datExpr colnames by sample_id.")
meta <- meta[ord_meta]
stopifnot(all(meta$sample_id == samples))

# datMeta_model 优先按 rownames 对齐；否则按当前顺序强行对齐
if (!is.null(rownames(as.data.frame(e$datMeta_model))) &&
    all(samples %in% rownames(as.data.frame(e$datMeta_model)))) {
  meta_model <- as.data.table(e$datMeta_model[ samples, , drop = FALSE ], keep.rownames = "sample_id_model")
} else {
  if (nrow(meta_model) != length(samples)) stop("datMeta_model row count != sample count.")
  meta_model[, sample_id_model := samples]
}

safe_fwrite(
  data.table(
    sample = samples,
    meta_sample_id = meta$sample_id,
    model_sample_id = meta_model$sample_id_model
  ),
  file.path(outdir, "01_alignment_audit.tsv")
)

# -----------------------------
# 2. 诊断标准化，仅保留 ASD / Control
# -----------------------------
if (!("Diagnosis" %in% colnames(meta))) stop("datMeta lacks Diagnosis.")

meta[, dx_std := standardize_dx(Diagnosis)]

dx_audit <- data.table(
  sample = meta$sample_id,
  Diagnosis_raw = as.character(meta$Diagnosis),
  dx_std = as.character(meta$dx_std),
  region = if ("region" %in% colnames(meta)) as.character(meta$region) else NA_character_,
  DxReg = if ("DxReg" %in% colnames(meta_model)) as.character(meta_model$DxReg) else NA_character_
)
safe_fwrite(dx_audit, file.path(outdir, "02_diagnosis_audit.tsv"))

keep_samples <- !is.na(meta$dx_std) & meta$dx_std %in% c("Control", "ASD")
expr <- expr[, keep_samples, drop = FALSE]
meta <- meta[keep_samples]
meta_model <- meta_model[keep_samples]

meta[, dx_std := factor(as.character(dx_std), levels = c("Control", "ASD"))]
safe_fwrite(meta[, .N, by = .(dx_std)][order(dx_std)], file.path(outdir, "03_filtered_dx_counts.tsv"))

# -----------------------------
# 3. 读取 frozen gene sets
# -----------------------------
candidate_genes <- read_gene_list(candidate_path)
reference_genes <- read_gene_list(reference_path)

writeLines(candidate_genes, con = file.path(outdir, "04_candidate_genes_used.txt"))
writeLines(reference_genes, con = file.path(outdir, "04_reference_genes_used.txt"))

# -----------------------------
# 4. ENSG -> SYMBOL 映射
# -----------------------------
raw_ids <- rownames(expr)
if (is.null(raw_ids)) stop("datExpr has no rownames.")

# ENSG00000000003.14_1 -> ENSG00000000003
ensg_core <- sub("\\..*$", "", raw_ids)

log_msg("Mapping ENSG to SYMBOL via org.Hs.eg.db ...")
symbol_map <- AnnotationDbi::mapIds(
  org.Hs.eg.db,
  keys = unique(ensg_core),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

mapped_symbol <- as.character(symbol_map[ensg_core])

row_annot <- data.table(
  raw_id = raw_ids,
  ensg_core = ensg_core,
  symbol = mapped_symbol
)
safe_fwrite(row_annot, file.path(outdir, "05_row_annotation.tsv"))

keep_gene <- !is.na(mapped_symbol) & mapped_symbol != ""
expr_use <- expr[keep_gene, , drop = FALSE]
mapped_symbol <- toupper(mapped_symbol[keep_gene])

# 同一 symbol 多行时取平均
expr_sum <- rowsum(expr_use, group = mapped_symbol, reorder = FALSE)
dup_n <- as.numeric(table(mapped_symbol)[rownames(expr_sum)])
expr_symbol <- expr_sum / dup_n

# gene-wise z-score across samples
expr_z <- t(scale(t(expr_symbol)))
expr_z[!is.finite(expr_z)] <- 0

# -----------------------------
# 6. program score
# -----------------------------
cand_res <- score_gene_set(candidate_genes, expr_z)
ref_res  <- score_gene_set(reference_genes, expr_z)

score_dt <- data.table(
  sample = colnames(expr_z),
  candidate_program_score = as.numeric(cand_res$score),
  reference_program_score = as.numeric(ref_res$score),
  candidate_minus_reference = as.numeric(cand_res$score - ref_res$score)
)

safe_fwrite(
  data.table(
    program = c("candidate", "reference"),
    n_input = c(length(candidate_genes), length(reference_genes)),
    n_hit = c(length(cand_res$hit), length(ref_res$hit))
  ),
  file.path(outdir, "06_gene_set_overlap.tsv")
)

writeLines(cand_res$hit, con = file.path(outdir, "06_candidate_hit_genes.txt"))
writeLines(ref_res$hit,  con = file.path(outdir, "06_reference_hit_genes.txt"))

# -----------------------------
# 7. 组装分析表
# -----------------------------
analysis_dt <- copy(score_dt)

analysis_dt[, dx_std := meta$dx_std[match(sample, meta$sample_id)]]
analysis_dt[, Subject := if ("Subject" %in% colnames(meta_model)) as.character(meta_model$Subject[match(sample, meta_model$sample_id_model)]) else
                        if ("subject" %in% colnames(meta)) as.character(meta$subject[match(sample, meta$sample_id)]) else NA_character_]
analysis_dt[, region := if ("region" %in% colnames(meta)) as.character(meta$region[match(sample, meta$sample_id)]) else NA_character_]
analysis_dt[, SeqBatch := if ("SeqBatch" %in% colnames(meta_model)) as.character(meta_model$SeqBatch[match(sample, meta_model$sample_id_model)]) else
                         if ("seq_batch" %in% colnames(meta)) as.character(meta$seq_batch[match(sample, meta$sample_id)]) else NA_character_]
analysis_dt[, Sex := if ("Sex" %in% colnames(meta_model)) as.character(meta_model$Sex[match(sample, meta_model$sample_id_model)]) else
                    if ("Sex" %in% colnames(meta)) as.character(meta$Sex[match(sample, meta$sample_id)]) else NA_character_]
analysis_dt[, Ancestry := if ("Ancestry" %in% colnames(meta_model)) as.character(meta_model$Ancestry[match(sample, meta_model$sample_id_model)]) else NA_character_]
analysis_dt[, Age := if ("Age" %in% colnames(meta_model)) as.numeric(meta_model$Age[match(sample, meta_model$sample_id_model)]) else
                    if ("Age" %in% colnames(meta)) as.numeric(meta$Age[match(sample, meta$sample_id)]) else NA_real_]
analysis_dt[, Age_sqd := if ("Age_sqd" %in% colnames(meta_model)) as.numeric(meta_model$Age_sqd[match(sample, meta_model$sample_id_model)]) else NA_real_]
analysis_dt[, PMI := if ("PMI" %in% colnames(meta_model)) as.numeric(meta_model$PMI[match(sample, meta_model$sample_id_model)]) else
                    if ("PMI" %in% colnames(meta)) as.numeric(meta$PMI[match(sample, meta$sample_id)]) else NA_real_]
analysis_dt[, RIN := if ("RIN" %in% colnames(meta_model)) as.numeric(meta_model$RIN[match(sample, meta_model$sample_id_model)]) else NA_real_]

safe_fwrite(analysis_dt, file.path(outdir, "07_sample_level_program_scores.tsv"))

# -----------------------------
# 8. mixed model
# -----------------------------
log_msg("Fitting mixed models ...")

res1 <- fit_one_program(copy(analysis_dt), "candidate_program_score")
res2 <- fit_one_program(copy(analysis_dt), "reference_program_score")
res3 <- fit_one_program(copy(analysis_dt), "candidate_minus_reference")

res <- rbindlist(list(res1, res2, res3), fill = TRUE)
res[, direction_expected := c("ASD_higher", "ASD_lower", "ASD_higher")]
res[, direction_observed := fifelse(beta_ASD_vs_Control > 0, "ASD_higher", "ASD_lower")]
res[, direction_match := direction_expected == direction_observed]
res[, BH_FDR := p.adjust(p_value, method = "BH")]

safe_fwrite(res, file.path(outdir, "08_program_mixed_model_results.tsv"))

# -----------------------------
# 9. plot
# -----------------------------
make_boxplot(
  analysis_dt, "candidate_program_score",
  "Candidate program score", "Gandal 2022 bulk: candidate program",
  file.path(outdir, "09_candidate_program_score.pdf"),
  file.path(outdir, "09_candidate_program_score.png")
)

make_boxplot(
  analysis_dt, "reference_program_score",
  "Reference program score", "Gandal 2022 bulk: reference program",
  file.path(outdir, "10_reference_program_score.pdf"),
  file.path(outdir, "10_reference_program_score.png")
)

make_boxplot(
  analysis_dt, "candidate_minus_reference",
  "Candidate - reference", "Gandal 2022 bulk: candidate minus reference",
  file.path(outdir, "11_candidate_minus_reference.pdf"),
  file.path(outdir, "11_candidate_minus_reference.png")
)

# -----------------------------
# 10. summary
# -----------------------------
summary_lines <- c(
  "Gandal 2022 bulk cortex program-level validation",
  "================================================",
  "",
  sprintf("Input RData: %s", rdata_path),
  sprintf("Candidate genes file: %s", candidate_path),
  sprintf("Reference genes file: %s", reference_path),
  sprintf("n samples retained (ASD + Control): %d", nrow(analysis_dt)),
  sprintf("n candidate hit genes: %d", length(cand_res$hit)),
  sprintf("n reference hit genes: %d", length(ref_res$hit)),
  "",
  "Program-level results:"
)

for (i in seq_len(nrow(res))) {
  summary_lines <- c(
    summary_lines,
    sprintf(
      "%s | model=%s | beta=%.4f | p=%.4g | BH=%.4g | expected=%s | observed=%s | match=%s | covars=%s",
      res$program[i],
      res$model_type[i],
      res$beta_ASD_vs_Control[i],
      res$p_value[i],
      res$BH_FDR[i],
      res$direction_expected[i],
      res$direction_observed[i],
      res$direction_match[i],
      res$covariates_used[i]
    )
  )
}

writeLines(summary_lines, con = file.path(outdir, "12_summary.txt"))
writeLines(capture.output(sessionInfo()), con = file.path(outdir, "13_sessionInfo.txt"))

log_msg("Done.")
