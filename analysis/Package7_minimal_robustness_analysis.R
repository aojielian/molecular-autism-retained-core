#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

DEFAULT_PROJECT_DIR <- normalizePath(Sys.getenv("MA_PROJECT_DIR", "."), mustWork = FALSE)
option_list <- list(
  make_option("--pkg4_dir", type = "character", default = file.path(DEFAULT_PROJECT_DIR, "results", "Package4_CrossCohortValidation_Velmeshev_v2")),
  make_option("--strict_dir", type = "character", default = file.path(DEFAULT_PROJECT_DIR, "results", "Package4b_DirectionConsistencySummary_v2_manual")),
  make_option("--pkg5_dir", type = "character", default = file.path(DEFAULT_PROJECT_DIR, "results", "Package5_DirectionalConcordance_GlobalSummary")),
  make_option("--outdir", type = "character", default = file.path(DEFAULT_PROJECT_DIR, "results", "Package7_MinimalRobustnessAnalysis"))
)
opt <- parse_args(OptionParser(option_list = option_list))
PKG4_DIR    <- opt$pkg4_dir
STRICT_DIR  <- opt$strict_dir
PKG5_DIR    <- opt$pkg5_dir
OUTDIR      <- opt$outdir
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

f_prop    <- file.path(PKG4_DIR, "15_validation_sample_score_defined_proportions.tsv")
f_shift   <- file.path(PKG4_DIR, "16_validation_sample_score_defined_shift.tsv")
f_strict  <- file.path(STRICT_DIR, "03_direction_consistency_summary_pretty.tsv")
f_shared  <- file.path(STRICT_DIR, "00_shared_comparable_samples.tsv")

need_files <- c(f_prop, f_shift, f_strict, f_shared)
for (f in need_files) {
  if (!file.exists(f)) stop("Missing input: ", f)
}

norm_name <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

find_col <- function(nms, patterns) {
  for (p in patterns) {
    hit <- which(grepl(p, nms, perl = TRUE))
    if (length(hit) > 0) return(nms[hit[1]])
  }
  return(NA_character_)
}

recode_dx <- function(x) {
  y <- tolower(trimws(as.character(x)))
  y2 <- gsub("[^a-z0-9]+", "_", y)
  out <- rep(NA_character_, length(y2))
  out[grepl("(^|_)(asd|case|patient|proband)(_|$)", y2, perl = TRUE)] <- "ASD"
  out[grepl("(^|_)(ctrl|control|td|typical|non_asd)(_|$)", y2, perl = TRUE)] <- "Control"
  out
}

read_dt <- function(f) {
  if (grepl("\\.csv(\\.gz)?$", f, perl = TRUE)) fread(f) else fread(f, sep = "\t")
}

standardize_sample_dx <- function(dt, file_tag = "") {
  setDT(dt)
  old0 <- colnames(dt)
  old_unique <- make.unique(old0, sep = "_origdup")
  setnames(dt, old = seq_along(old0), new = old_unique)
  new0 <- norm_name(old_unique)
  new_unique <- make.unique(new0, sep = "_dup")
  setnames(dt, old = seq_along(new_unique), new = new_unique)

  dx_col <- find_col(colnames(dt), c("^diagnosis$", "^dx$", "^group$", "^status$", "^condition$"))
  sample_col <- find_col(colnames(dt), c(
    "^sample$", "^sample_id$", "^sampleid$",
    "^donor$", "^donor_id$", "^donorid$",
    "^subject$", "^subject_id$", "^subjectid$",
    "^individual$", "^individual_id$",
    "^iid$", "^orig_ident$"
  ))

  if (is.na(dx_col)) stop("No dx column found in ", file_tag)
  if (is.na(sample_col)) stop("No sample column found in ", file_tag)

  dt[, dx2 := recode_dx(get(dx_col))]
  dt <- dt[!is.na(dx2)]
  dt[, sample2 := as.character(get(sample_col))]
  list(dt = dt, dx_col = dx_col, sample_col = sample_col)
}

summarize_metric <- function(dt, metric_col, metric_label, expected_direction, source_set) {
  if (is.na(metric_col) || !(metric_col %in% colnames(dt))) {
    return(data.table(
      source_set = source_set,
      metric_label = metric_label,
      expected_direction = expected_direction,
      observed_direction = NA_character_,
      direction_match = NA,
      n_asd = NA_integer_,
      n_control = NA_integer_,
      lm_beta_asd = NA_real_,
      lm_p = NA_real_,
      wilcox_p = NA_real_,
      aligned_effect = NA_real_
    ))
  }

  sub <- copy(dt)[, .(sample2, dx2, value = suppressWarnings(as.numeric(get(metric_col))))]
  sub <- sub[!is.na(value)]
  sub <- sub[, .(value = mean(value, na.rm = TRUE)), by = .(sample2, dx2)]

  asd_vals  <- sub[dx2 == "ASD", value]
  ctrl_vals <- sub[dx2 == "Control", value]

  sub[, dx2 := factor(dx2, levels = c("Control", "ASD"))]
  fit <- tryCatch(lm(value ~ dx2, data = sub), error = function(e) NULL)

  beta <- NA_real_
  p_lm <- NA_real_
  if (!is.null(fit)) {
    cf <- summary(fit)$coefficients
    if ("dx2ASD" %in% rownames(cf)) {
      beta <- unname(cf["dx2ASD", "Estimate"])
      p_lm <- unname(cf["dx2ASD", "Pr(>|t|)"])
    }
  }

  p_w <- tryCatch(wilcox.test(value ~ dx2, data = sub, exact = FALSE)$p.value, error = function(e) NA_real_)
  primary_effect <- ifelse(is.na(beta), mean(asd_vals) - mean(ctrl_vals), beta)

  observed_direction <- if (is.na(primary_effect)) {
    NA_character_
  } else if (primary_effect > 0) {
    "ASD > Control"
  } else if (primary_effect < 0) {
    "ASD < Control"
  } else {
    "No difference"
  }

  direction_match <- observed_direction == expected_direction
  expected_sign <- ifelse(expected_direction == "ASD > Control", 1, -1)
  aligned_effect <- primary_effect * expected_sign

  data.table(
    source_set = source_set,
    metric_label = metric_label,
    expected_direction = expected_direction,
    observed_direction = observed_direction,
    direction_match = direction_match,
    n_asd = length(asd_vals),
    n_control = length(ctrl_vals),
    lm_beta_asd = beta,
    lm_p = p_lm,
    wilcox_p = p_w,
    aligned_effect = aligned_effect
  )
}

calc_global <- function(dt) {
  x <- copy(dt)
  x[, nominal_p := ifelse(!is.na(lm_p), lm_p, wilcox_p)]
  x[, direction_match := as.logical(direction_match)]
  x <- x[!is.na(direction_match) & !is.na(nominal_p)]

  n_total <- nrow(x)
  n_matched <- sum(x$direction_match)
  n_mismatched <- n_total - n_matched

  sign_dir <- ifelse(x$direction_match, 1, -1)
  p_clamped <- pmax(pmin(x$nominal_p, 1 - 1e-15), 1e-15)
  z_abs <- qnorm(1 - p_clamped / 2)
  z_signed <- sign_dir * z_abs

  stouffer_z <- sum(z_signed) / sqrt(n_total)
  stouffer_one_sided_p <- 1 - pnorm(stouffer_z)
  stouffer_two_sided_p <- 2 * pnorm(-abs(stouffer_z))

  data.table(
    n_total_metrics = n_total,
    n_direction_matched = n_matched,
    n_direction_mismatched = n_mismatched,
    direction_match_fraction = n_matched / n_total,
    binomial_one_sided_p = pbinom(n_matched - 1, n_total, 0.5, lower.tail = FALSE),
    binomial_two_sided_p = binom.test(n_matched, n_total, 0.5)$p.value,
    stouffer_z = stouffer_z,
    stouffer_one_sided_p = stouffer_one_sided_p,
    stouffer_two_sided_p = stouffer_two_sided_p
  )
}

prop0  <- read_dt(f_prop)
shift0 <- read_dt(f_shift)

propx  <- standardize_sample_dx(prop0,  "15_validation_sample_score_defined_proportions")
shiftx <- standardize_sample_dx(shift0, "16_validation_sample_score_defined_shift")

prop  <- propx$dt
shift <- shiftx$dt

prop_candidate_col <- find_col(colnames(prop), c(
  "^prop_candidate$",
  "^candidate_proportion$",
  "^candidate_prop$",
  "^score_defined_candidate_proportion$"
))

prop_reference_col <- find_col(colnames(prop), c(
  "^prop_reference$",
  "^reference_proportion$",
  "^reference_prop$",
  "^score_defined_reference_proportion$"
))

shift_col <- find_col(colnames(shift), c(
  "^score_defined_shift$",
  "^candidate_vs_reference_shift$",
  "^score_defined_candidate_vs_reference_shift$"
))

if (!is.na(prop_candidate_col) && !is.na(prop_reference_col)) {
  prop[, log2_candidate_reference_ratio := log2((get(prop_candidate_col) + 1e-6) / (get(prop_reference_col) + 1e-6))]
}
log2_col <- "log2_candidate_reference_ratio"

strict_samples <- unique(fread(f_shared, sep = "\t")$sample2)
broader_samples <- intersect(unique(prop$sample2), unique(shift$sample2))

prop_broader  <- prop[sample2 %in% broader_samples]
shift_broader <- shift[sample2 %in% broader_samples]

prop_strict   <- prop[sample2 %in% strict_samples]
shift_strict  <- shift[sample2 %in% strict_samples]

sample_count_tbl <- rbindlist(list(
  unique(prop_broader[,  .(source_set = "broader_15_16_intersection", sample2, dx2)]),
  unique(prop_strict[,   .(source_set = "strict_12_15_16_intersection", sample2, dx2)])
), fill = TRUE)

sample_count_sum <- sample_count_tbl[, .N, by = .(source_set, dx2)]
fwrite(sample_count_sum, file.path(OUTDIR, "00_harmonization_sample_counts.tsv"), sep = "\t")

harm_tbl <- rbindlist(list(
  summarize_metric(prop_broader,  prop_candidate_col, "Score-defined cells: candidate proportion",      "ASD > Control", "broader_15_16_intersection"),
  summarize_metric(prop_broader,  prop_reference_col, "Score-defined cells: reference proportion",      "ASD < Control", "broader_15_16_intersection"),
  summarize_metric(prop_broader,  log2_col,           "Score-defined cells: log2(candidate/reference)", "ASD > Control", "broader_15_16_intersection"),
  summarize_metric(shift_broader, shift_col,          "Score-defined cells: candidate vs reference shift","ASD > Control", "broader_15_16_intersection"),
  summarize_metric(prop_strict,   prop_candidate_col, "Score-defined cells: candidate proportion",      "ASD > Control", "strict_12_15_16_intersection"),
  summarize_metric(prop_strict,   prop_reference_col, "Score-defined cells: reference proportion",      "ASD < Control", "strict_12_15_16_intersection"),
  summarize_metric(prop_strict,   log2_col,           "Score-defined cells: log2(candidate/reference)", "ASD > Control", "strict_12_15_16_intersection"),
  summarize_metric(shift_strict,  shift_col,          "Score-defined cells: candidate vs reference shift","ASD > Control", "strict_12_15_16_intersection")
), fill = TRUE)

fwrite(harm_tbl, file.path(OUTDIR, "01_score_defined_harmonization_sensitivity.tsv"), sep = "\t")

harm_wide <- dcast(
  harm_tbl,
  metric_label + expected_direction ~ source_set,
  value.var = c("observed_direction", "direction_match", "lm_beta_asd", "lm_p", "wilcox_p", "aligned_effect")
)

harm_wide[, direction_stable := observed_direction_broader_15_16_intersection == observed_direction_strict_12_15_16_intersection]
harm_wide[, match_stable := direction_match_broader_15_16_intersection == direction_match_strict_12_15_16_intersection]

fwrite(harm_wide, file.path(OUTDIR, "02_score_defined_direction_stability.tsv"), sep = "\t")

strict_pretty <- fread(f_strict, sep = "\t")
base_global <- calc_global(strict_pretty)
base_global[, omitted_metric := "none"]
base_global[, source_set := "strict_12_15_16_intersection"]

loo_tbl <- rbindlist(lapply(seq_len(nrow(strict_pretty)), function(i) {
  sub <- strict_pretty[-i]
  x <- calc_global(sub)
  x[, omitted_metric := strict_pretty$metric_label[i]]
  x[, source_set := "strict_12_15_16_intersection"]
  x
}), fill = TRUE)

loo_tbl <- rbindlist(list(base_global, loo_tbl), fill = TRUE)
fwrite(loo_tbl, file.path(OUTDIR, "03_leave_one_metric_out_global_summary.tsv"), sep = "\t")

sink(file.path(OUTDIR, "04_package7_summary.txt"))
cat("Package 7 summary\n")
cat("=================\n\n")
cat("Harmonization sensitivity (score-defined metrics only):\n")
print(sample_count_sum)
cat("\nDirection stability across broader vs strict sets:\n")
cat("- stable observed direction count = ", sum(harm_wide$direction_stable, na.rm = TRUE), "/", nrow(harm_wide), "\n", sep = "")
cat("- stable direction_match count = ", sum(harm_wide$match_stable, na.rm = TRUE), "/", nrow(harm_wide), "\n\n", sep = "")
cat("Leave-one-metric-out global summary:\n")
print(loo_tbl)
sink()

message("Package 7 done.")
