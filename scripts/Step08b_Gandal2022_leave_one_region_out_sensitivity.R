#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript Step08b_Gandal2022_leave_one_region_out_sensitivity.R <sample_level_scores.tsv> <outdir>")
}

infile <- args[1]
outdir <- args[2]

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(nlme)
  library(ggplot2)
})

logfile <- file.path(outdir, "00_step08b_run.log")
log_msg <- function(...) {
  msg <- sprintf(...)
  line <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg)
  cat(line, "\n")
  cat(line, "\n", file = logfile, append = TRUE)
}

safe_fwrite <- function(x, file) {
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  fwrite(x, file = file, sep = "\t", quote = FALSE, na = "NA")
}

usable_factor <- function(x) {
  x <- x[!is.na(x)]
  length(x) > 0 && length(unique(x)) > 1
}

usable_numeric <- function(x) {
  x <- x[is.finite(x)]
  length(x) > 0 && length(unique(x)) > 1
}

fit_one_program <- function(df, score_col) {
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

  covar_sets <- list(
    keep_covars,
    setdiff(keep_covars, "Age_sqd"),
    intersect(keep_covars, c("region", "SeqBatch", "Sex", "Ancestry", "Age", "PMI", "RIN")),
    intersect(keep_covars, c("region", "SeqBatch", "Sex", "Age", "PMI", "RIN")),
    intersect(keep_covars, c("region", "Sex", "Age", "PMI", "RIN")),
    intersect(keep_covars, c("region", "Sex", "Age"))
  )

  for (covars_use in covar_sets) {
    need_cols <- unique(c(score_col, "dx_std", "Subject", covars_use))
    dat <- copy(df)[, ..need_cols]
    dat <- dat[complete.cases(dat)]

    if (nrow(dat) < 8) next
    if (length(unique(dat$dx_std)) < 2) next

    dat[, dx_std := factor(dx_std, levels = c("Control", "ASD"))]
    dat[, Subject := as.factor(Subject)]

    for (cc in covars_use) {
      if (!is.numeric(dat[[cc]])) dat[[cc]] <- as.factor(dat[[cc]])
    }

    rhs <- c("dx_std", covars_use)
    form <- as.formula(paste(score_col, "~", paste(rhs, collapse = " + ")))

    use_random_subject <- length(unique(dat$Subject)) < nrow(dat)

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

    fit_lm <- tryCatch(lm(form, data = dat), error = function(e) NULL)
    if (!is.null(fit_lm)) {
      sm <- summary(fit_lm)$coefficients
      term <- "dx_stdASD"
      if (term %in% rownames(sm)) {
        return(data.table(
          model_type = "lm_fallback",
          covariates_used = paste(covars_use, collapse = ";"),
          n_samples = nrow(dat),
          n_subjects = uniqueN(dat$Subject),
          beta_ASD_vs_Control = unname(sm[term, "Estimate"]),
          se = unname(sm[term, "Std. Error"]),
          t_value = unname(sm[term, "t value"]),
          p_value = unname(sm[term, grep("^Pr", colnames(sm))])
        ))
      }
    }
  }

  data.table(
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

make_forest_plot <- function(dt, program_name, out_pdf, out_png) {
  d <- copy(dt)[program == program_name]
  if (nrow(d) == 0) return(invisible(NULL))

  d[, comparison := factor(comparison, levels = rev(unique(comparison)))]
  d[, ci_low := beta_ASD_vs_Control - 1.96 * se]
  d[, ci_high := beta_ASD_vs_Control + 1.96 * se]

  p <- ggplot(d, aes(x = beta_ASD_vs_Control, y = comparison)) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.18) +
    geom_point(size = 2) +
    theme_bw(base_size = 12) +
    labs(
      x = "ASD vs Control beta (95% CI)",
      y = NULL,
      title = paste0(program_name, ": leave-one-region-out sensitivity")
    )

  ggsave(out_pdf, p, width = 6.8, height = 5.2)
  ggsave(out_png, p, width = 6.8, height = 5.2, dpi = 300)
}

log_msg("Reading input: %s", infile)
dt <- fread(infile)
dt <- as.data.table(dt)

required_cols <- c(
  "sample", "dx_std", "Subject", "region",
  "candidate_program_score", "reference_program_score", "candidate_minus_reference"
)
missing_cols <- setdiff(required_cols, colnames(dt))
if (length(missing_cols) > 0) {
  stop(sprintf("Missing required columns: %s", paste(missing_cols, collapse = ", ")))
}

# 标准化 dx
dt[, dx_std := as.character(dx_std)]
dt[dx_std %in% c("CTL", "Control", "ctrl", "ctl"), dx_std := "Control"]
dt[dx_std %in% c("ASD", "autism", "case"), dx_std := "ASD"]
dt <- dt[dx_std %in% c("Control", "ASD")]

# 标准化部分字段类型
dt[, Subject := as.character(Subject)]
dt[, region := as.character(region)]
if ("SeqBatch" %in% colnames(dt)) dt[, SeqBatch := as.character(SeqBatch)]
if ("Sex" %in% colnames(dt)) dt[, Sex := as.character(Sex)]
if ("Ancestry" %in% colnames(dt)) dt[, Ancestry := as.character(Ancestry)]

safe_fwrite(dt[, .N, by = .(dx_std)][order(dx_std)], file.path(outdir, "01_input_dx_counts.tsv"))
safe_fwrite(dt[, .N, by = .(region)][order(region)], file.path(outdir, "02_input_region_counts.tsv"))
safe_fwrite(dt[, .N, by = .(dx_std, region)][order(region, dx_std)], file.path(outdir, "03_input_dx_by_region_counts.tsv"))

programs <- c("candidate_program_score", "reference_program_score", "candidate_minus_reference")
expected_direction <- c(
  candidate_program_score = "ASD_higher",
  reference_program_score = "ASD_lower",
  candidate_minus_reference = "ASD_higher"
)

# 先跑 full-data baseline
log_msg("Fitting full-data baseline models")
full_res <- rbindlist(lapply(programs, function(pgm) {
  rr <- fit_one_program(copy(dt), pgm)
  rr[, program := pgm]
  rr[, omitted_region := "NONE"]
  rr[, comparison := "Full data"]
  rr
}), fill = TRUE)

# 再逐个 leave-one-region-out
regions <- sort(unique(dt$region))
log_msg("Detected %d regions for leave-one-region-out", length(regions))

loo_res <- rbindlist(lapply(regions, function(rg) {
  dsub <- dt[region != rg]
  rbindlist(lapply(programs, function(pgm) {
    rr <- fit_one_program(copy(dsub), pgm)
    rr[, program := pgm]
    rr[, omitted_region := rg]
    rr[, comparison := paste0("Drop ", rg)]
    rr
  }), fill = TRUE)
}), fill = TRUE)

res <- rbindlist(list(full_res, loo_res), fill = TRUE)

res[, direction_expected := expected_direction[program]]
res[, direction_observed := fifelse(beta_ASD_vs_Control > 0, "ASD_higher", "ASD_lower")]
res[, direction_match := direction_expected == direction_observed]

# leave-one-region-out 内部单独做 FDR
res[, p_value_BH_within_all := p.adjust(p_value, method = "BH")]
res[omitted_region != "NONE", p_value_BH_within_LOO := p.adjust(p_value, method = "BH")]
res[omitted_region == "NONE", p_value_BH_within_LOO := NA_real_]

safe_fwrite(res, file.path(outdir, "04_leave_one_region_out_results.tsv"))

# 稳定性摘要
baseline <- res[omitted_region == "NONE", .(
  program, baseline_beta = beta_ASD_vs_Control, baseline_p = p_value,
  baseline_model = model_type, baseline_covars = covariates_used
)]

stab <- merge(
  res[omitted_region != "NONE"],
  baseline,
  by = "program",
  all.x = TRUE
)

stab[, beta_diff_from_baseline := beta_ASD_vs_Control - baseline_beta]
stab[, abs_beta_diff_from_baseline := abs(beta_diff_from_baseline)]
stab[, sign_same_as_baseline := sign(beta_ASD_vs_Control) == sign(baseline_beta)]

stab_summary <- stab[, .(
  n_regions_tested = .N,
  n_direction_match = sum(direction_match, na.rm = TRUE),
  direction_match_fraction = mean(direction_match, na.rm = TRUE),
  n_same_sign_as_baseline = sum(sign_same_as_baseline, na.rm = TRUE),
  same_sign_fraction = mean(sign_same_as_baseline, na.rm = TRUE),
  min_beta = min(beta_ASD_vs_Control, na.rm = TRUE),
  max_beta = max(beta_ASD_vs_Control, na.rm = TRUE),
  min_p = min(p_value, na.rm = TRUE),
  max_p = max(p_value, na.rm = TRUE),
  max_abs_beta_shift = max(abs_beta_diff_from_baseline, na.rm = TRUE)
), by = program]

safe_fwrite(stab, file.path(outdir, "05_leave_one_region_out_vs_baseline.tsv"))
safe_fwrite(stab_summary, file.path(outdir, "06_leave_one_region_out_stability_summary.tsv"))

# 画图
make_forest_plot(
  res, "candidate_program_score",
  file.path(outdir, "07_candidate_program_leave_one_region_out_forest.pdf"),
  file.path(outdir, "07_candidate_program_leave_one_region_out_forest.png")
)

make_forest_plot(
  res, "reference_program_score",
  file.path(outdir, "08_reference_program_leave_one_region_out_forest.pdf"),
  file.path(outdir, "08_reference_program_leave_one_region_out_forest.png")
)

make_forest_plot(
  res, "candidate_minus_reference",
  file.path(outdir, "09_candidate_minus_reference_leave_one_region_out_forest.pdf"),
  file.path(outdir, "09_candidate_minus_reference_leave_one_region_out_forest.png")
)

# summary 文本
summary_lines <- c(
  "Gandal 2022 leave-one-region-out sensitivity",
  "============================================",
  "",
  sprintf("Input score table: %s", infile),
  sprintf("Total ASD + Control samples: %d", nrow(dt)),
  sprintf("Unique subjects: %d", uniqueN(dt$Subject)),
  sprintf("Unique regions: %d", uniqueN(dt$region)),
  "",
  "Baseline full-data results:"
)

for (i in seq_len(nrow(full_res))) {
  summary_lines <- c(
    summary_lines,
    sprintf(
      "%s | beta=%.4f | p=%.4g | expected=%s | observed=%s | match=%s | covars=%s",
      full_res$program[i],
      full_res$beta_ASD_vs_Control[i],
      full_res$p_value[i],
      expected_direction[full_res$program[i]],
      ifelse(full_res$beta_ASD_vs_Control[i] > 0, "ASD_higher", "ASD_lower"),
      expected_direction[full_res$program[i]] == ifelse(full_res$beta_ASD_vs_Control[i] > 0, "ASD_higher", "ASD_lower"),
      full_res$covariates_used[i]
    )
  )
}

summary_lines <- c(summary_lines, "", "Leave-one-region-out stability summary:")

for (i in seq_len(nrow(stab_summary))) {
  summary_lines <- c(
    summary_lines,
    sprintf(
      "%s | regions_tested=%d | direction_match_fraction=%.3f | same_sign_fraction=%.3f | min_beta=%.4f | max_beta=%.4f | min_p=%.4g | max_p=%.4g | max_abs_beta_shift=%.4f",
      stab_summary$program[i],
      stab_summary$n_regions_tested[i],
      stab_summary$direction_match_fraction[i],
      stab_summary$same_sign_fraction[i],
      stab_summary$min_beta[i],
      stab_summary$max_beta[i],
      stab_summary$min_p[i],
      stab_summary$max_p[i],
      stab_summary$max_abs_beta_shift[i]
    )
  )
}

writeLines(summary_lines, file.path(outdir, "10_summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(outdir, "11_sessionInfo.txt"))

log_msg("Done.")
