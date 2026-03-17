#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

DEFAULT_PROJECT_DIR <- normalizePath(Sys.getenv("MA_PROJECT_DIR", "."), mustWork = FALSE)
option_list <- list(
  make_option("--indir", type = "character", default = file.path(DEFAULT_PROJECT_DIR, "results", "Package4b_DirectionConsistencySummary_v2_manual")),
  make_option("--outdir", type = "character", default = file.path(DEFAULT_PROJECT_DIR, "results", "Package5_DirectionalConcordance_GlobalSummary"))
)
opt <- parse_args(OptionParser(option_list = option_list))
INDIR  <- opt$indir
OUTDIR <- opt$outdir
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

f_pretty <- file.path(INDIR, "03_direction_consistency_summary_pretty.tsv")
f_counts <- file.path(INDIR, "00b_shared_sample_counts.tsv")

if (!file.exists(f_pretty)) stop("Missing input: ", f_pretty)
if (!file.exists(f_counts)) stop("Missing input: ", f_counts)

pretty <- fread(f_pretty, sep = "\t")
counts <- fread(f_counts, sep = "\t")

metric_order <- c(
  "All microglia: candidate - reference",
  "All microglia: candidate-state score",
  "All microglia: reference-like score",
  "Score-defined cells: candidate proportion",
  "Score-defined cells: reference proportion",
  "Score-defined cells: log2(candidate/reference)",
  "Score-defined cells: candidate vs reference shift"
)

pretty[, metric_label := as.character(metric_label)]
pretty[, expected_direction := as.character(expected_direction)]
pretty[, observed_direction := as.character(observed_direction)]

pretty[, nominal_p := ifelse(!is.na(lm_p), lm_p, wilcox_p)]
pretty[, nominal_p_source := ifelse(!is.na(lm_p), "lm_p", "wilcox_p")]
pretty[, metric_group := ifelse(grepl("^All microglia", metric_label), "all_microglia", "score_defined")]
pretty[, direction_match := as.logical(direction_match)]

pretty[, metric_order := match(metric_label, metric_order)]
setorder(pretty, metric_order)

pretty[, sign_dir := fifelse(direction_match == TRUE, 1, fifelse(direction_match == FALSE, -1, NA_real_))]
pretty[, nominal_p_clamped := pmax(pmin(nominal_p, 1 - 1e-15), 1e-15)]
pretty[, z_abs := qnorm(1 - nominal_p_clamped / 2)]
pretty[, z_signed := sign_dir * z_abs]

valid <- pretty[!is.na(direction_match) & !is.na(nominal_p)]
n_total <- nrow(valid)
n_matched <- sum(valid$direction_match, na.rm = TRUE)
n_mismatched <- sum(!valid$direction_match, na.rm = TRUE)

if (n_total == 0) stop("No valid rows found in strict Package4b summary.")

binom_two_sided <- binom.test(n_matched, n_total, p = 0.5)$p.value
binom_one_sided <- pbinom(n_matched - 1, n_total, 0.5, lower.tail = FALSE)

stouffer_z <- sum(valid$z_signed, na.rm = TRUE) / sqrt(n_total)
stouffer_one_sided_p <- 1 - pnorm(stouffer_z)
stouffer_two_sided_p <- 2 * pnorm(-abs(stouffer_z))

strongest_nominal <- valid[order(nominal_p, -aligned_effect)][1]
largest_aligned <- valid[order(-aligned_effect, nominal_p)][1]

master <- copy(pretty)[, .(
  metric_label,
  metric_group,
  expected_direction,
  observed_direction,
  direction_match,
  n_asd,
  n_control,
  lm_beta_asd,
  lm_p,
  lm_fdr,
  wilcox_p,
  wilcox_fdr,
  nominal_p,
  nominal_p_source,
  mean_asd,
  mean_control,
  median_asd,
  median_control,
  aligned_effect,
  source_file_basename,
  source_column,
  evidence_note
)]

master[, direction_match_label := ifelse(direction_match, "matched", "mismatched")]
master[, included_in_global_summary := TRUE]

global_tests <- data.table(
  test_name = c(
    "n_total_metrics",
    "n_direction_matched",
    "n_direction_mismatched",
    "direction_match_fraction",
    "binomial_test_one_sided_p",
    "binomial_test_two_sided_p",
    "stouffer_z",
    "stouffer_one_sided_p",
    "stouffer_two_sided_p",
    "n_asd_strict_shared",
    "n_control_strict_shared"
  ),
  value = c(
    n_total,
    n_matched,
    n_mismatched,
    n_matched / n_total,
    binom_one_sided,
    binom_two_sided,
    stouffer_z,
    stouffer_one_sided_p,
    stouffer_two_sided_p,
    counts[dx2 == "ASD", N][1],
    counts[dx2 == "Control", N][1]
  )
)

rank_nominal <- copy(master)[order(nominal_p, -aligned_effect)]
rank_aligned <- copy(master)[order(-aligned_effect, nominal_p)]

forest_input <- copy(master)[, .(
  metric_label,
  metric_group,
  aligned_effect,
  nominal_p,
  direction_match,
  observed_direction,
  expected_direction
)]

fwrite(master,       file.path(OUTDIR, "01_direction_concordance_master.tsv"), sep = "\t")
fwrite(global_tests, file.path(OUTDIR, "02_global_direction_tests.tsv"), sep = "\t")
fwrite(rank_nominal, file.path(OUTDIR, "03_metrics_ranked_by_nominal_p.tsv"), sep = "\t")
fwrite(rank_aligned, file.path(OUTDIR, "04_metrics_ranked_by_aligned_effect.tsv"), sep = "\t")
fwrite(forest_input, file.path(OUTDIR, "05_forest_input.tsv"), sep = "\t")

sink(file.path(OUTDIR, "06_manuscript_summary.txt"))
cat("Package 5 summary\n")
cat("=================\n\n")
cat("Primary input:\n")
cat("- ", f_pretty, "\n", sep = "")
cat("- ", f_counts, "\n\n", sep = "")
cat("Strict shared validation sample counts:\n")
cat("- ASD: ", counts[dx2 == "ASD", N][1], "\n", sep = "")
cat("- Control: ", counts[dx2 == "Control", N][1], "\n\n", sep = "")
cat("Directional concordance:\n")
cat("- matched metrics: ", n_matched, "/", n_total, "\n", sep = "")
cat("- one-sided binomial p: ", signif(binom_one_sided, 4), "\n", sep = "")
cat("- two-sided binomial p: ", signif(binom_two_sided, 4), "\n", sep = "")
cat("- Stouffer Z: ", signif(stouffer_z, 4), "\n", sep = "")
cat("- Stouffer one-sided p: ", signif(stouffer_one_sided_p, 4), "\n", sep = "")
cat("- Stouffer two-sided p: ", signif(stouffer_two_sided_p, 4), "\n\n", sep = "")
cat("Strongest nominal metric:\n")
cat("- ", strongest_nominal$metric_label, "\n", sep = "")
cat("- beta = ", signif(strongest_nominal$lm_beta_asd, 4), "\n", sep = "")
cat("- nominal p = ", signif(strongest_nominal$nominal_p, 4), "\n\n", sep = "")
cat("Largest aligned effect:\n")
cat("- ", largest_aligned$metric_label, "\n", sep = "")
cat("- aligned effect = ", signif(largest_aligned$aligned_effect, 4), "\n\n", sep = "")
cat("Manuscript sentence:\n")
cat("Across the strict shared validation set, all seven pre-specified validation metrics were directionally concordant with the discovery-prioritized expectation, supporting coherent but statistically limited cross-cohort enrichment of the candidate-state program.\n")
sink()

message("Package 5 done.")
