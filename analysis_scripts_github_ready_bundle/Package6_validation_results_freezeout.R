#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

DEFAULT_PROJECT_DIR <- normalizePath(Sys.getenv("MA_PROJECT_DIR", "."), mustWork = FALSE)
option_list <- list(
  make_option("--strict_dir", type = "character", default = file.path(DEFAULT_PROJECT_DIR, "results", "Package4b_DirectionConsistencySummary_v2_manual")),
  make_option("--pkg5_dir", type = "character", default = file.path(DEFAULT_PROJECT_DIR, "results", "Package5_DirectionalConcordance_GlobalSummary")),
  make_option("--outdir", type = "character", default = file.path(DEFAULT_PROJECT_DIR, "results", "Package6_ValidationResults_Freezeout"))
)
opt <- parse_args(OptionParser(option_list = option_list))
STRICT_DIR  <- opt$strict_dir
PKG5_DIR    <- opt$pkg5_dir
OUTDIR      <- opt$outdir
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

f_master   <- file.path(PKG5_DIR, "01_direction_concordance_master.tsv")
f_global   <- file.path(PKG5_DIR, "02_global_direction_tests.tsv")
f_counts   <- file.path(STRICT_DIR, "00b_shared_sample_counts.tsv")
f_samples  <- file.path(STRICT_DIR, "00_shared_comparable_samples.tsv")
f_detect   <- file.path(STRICT_DIR, "01_detected_metric_columns.tsv")
f_pretty   <- file.path(STRICT_DIR, "03_direction_consistency_summary_pretty.tsv")

need_files <- c(f_master, f_global, f_counts, f_samples, f_detect, f_pretty)
for (f in need_files) {
  if (!file.exists(f)) stop("Missing input: ", f)
}

master  <- fread(f_master,  sep = "\t")
global  <- fread(f_global,  sep = "\t")
counts  <- fread(f_counts,  sep = "\t")
samples <- fread(f_samples, sep = "\t")
detect  <- fread(f_detect,  sep = "\t")
pretty  <- fread(f_pretty,  sep = "\t")

master[, primary_summary_set := "strict_shared_validation_set"]
master[, included_in_primary_harmonized_summary := TRUE]
master[, summary_version := "Package4b_v2_manual_strict_shared"]

supp_table <- copy(master)[, .(
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
  aligned_effect,
  source_file_basename,
  source_column,
  evidence_note,
  primary_summary_set,
  included_in_primary_harmonized_summary,
  summary_version
)]

source_trace <- copy(master)[, .(
  metric_label,
  metric_group,
  source_file_basename,
  source_column,
  strict_summary_file = basename(f_pretty),
  metric_detection_file = basename(f_detect),
  included_in_primary_harmonized_summary,
  primary_summary_set,
  summary_version
)]

sample_counts <- copy(counts)
sample_counts[, sample_set := "strict_shared_validation_set"]

sample_list <- copy(samples)
sample_list[, sample_set := "strict_shared_validation_set"]

global_wide <- dcast(global, . ~ test_name, value.var = "value")

get_val <- function(metric_name, col_name) {
  x <- master[metric_label == metric_name, get(col_name)]
  if (length(x) == 0) return(NA)
  x[1]
}

manu_numbers <- data.table(
  item = c(
    "strict_shared_n_asd",
    "strict_shared_n_control",
    "n_total_metrics",
    "n_direction_matched",
    "n_direction_mismatched",
    "direction_match_fraction",
    "binomial_test_one_sided_p",
    "binomial_test_two_sided_p",
    "stouffer_z",
    "stouffer_one_sided_p",
    "stouffer_two_sided_p",
    "strongest_metric_label",
    "strongest_metric_beta",
    "strongest_metric_nominal_p",
    "candidate_proportion_beta",
    "candidate_proportion_nominal_p",
    "candidate_state_score_beta",
    "candidate_state_score_nominal_p",
    "candidate_minus_reference_beta",
    "candidate_minus_reference_nominal_p",
    "log2_ratio_beta",
    "log2_ratio_nominal_p"
  ),
  value = c(
    counts[dx2 == "ASD", N][1],
    counts[dx2 == "Control", N][1],
    global[test_name == "n_total_metrics", value][1],
    global[test_name == "n_direction_matched", value][1],
    global[test_name == "n_direction_mismatched", value][1],
    global[test_name == "direction_match_fraction", value][1],
    global[test_name == "binomial_test_one_sided_p", value][1],
    global[test_name == "binomial_test_two_sided_p", value][1],
    global[test_name == "stouffer_z", value][1],
    global[test_name == "stouffer_one_sided_p", value][1],
    global[test_name == "stouffer_two_sided_p", value][1],
    master[order(nominal_p, -aligned_effect), metric_label][1],
    master[order(nominal_p, -aligned_effect), lm_beta_asd][1],
    master[order(nominal_p, -aligned_effect), nominal_p][1],
    get_val("Score-defined cells: candidate proportion", "lm_beta_asd"),
    get_val("Score-defined cells: candidate proportion", "nominal_p"),
    get_val("All microglia: candidate-state score", "lm_beta_asd"),
    get_val("All microglia: candidate-state score", "nominal_p"),
    get_val("All microglia: candidate - reference", "lm_beta_asd"),
    get_val("All microglia: candidate - reference", "nominal_p"),
    get_val("Score-defined cells: log2(candidate/reference)", "lm_beta_asd"),
    get_val("Score-defined cells: log2(candidate/reference)", "nominal_p")
  )
)

fwrite(master,       file.path(OUTDIR, "01_validation_master_metrics.tsv"), sep = "\t")
fwrite(supp_table,   file.path(OUTDIR, "02_validation_supplementary_table.tsv"), sep = "\t")
fwrite(source_trace, file.path(OUTDIR, "03_validation_source_trace.tsv"), sep = "\t")
fwrite(sample_counts,file.path(OUTDIR, "04_validation_sample_set_counts.tsv"), sep = "\t")
fwrite(sample_list,  file.path(OUTDIR, "05_validation_sample_set_list.tsv"), sep = "\t")
fwrite(manu_numbers, file.path(OUTDIR, "06_validation_manuscript_numbers.tsv"), sep = "\t")
fwrite(global_wide,  file.path(OUTDIR, "07_validation_global_tests_wide.tsv"), sep = "\t")

sink(file.path(OUTDIR, "08_validation_results_blurb.txt"))
cat("Package 6 freeze-out summary\n")
cat("============================\n\n")
cat("Primary summary set: strict_shared_validation_set\n")
cat("ASD samples: ", counts[dx2 == "ASD", N][1], "\n", sep = "")
cat("Control samples: ", counts[dx2 == "Control", N][1], "\n\n", sep = "")
cat("Global direction summary:\n")
cat("- matched metrics: ", global[test_name == "n_direction_matched", value][1], "/", global[test_name == "n_total_metrics", value][1], "\n", sep = "")
cat("- one-sided binomial p: ", signif(global[test_name == "binomial_test_one_sided_p", value][1], 4), "\n", sep = "")
cat("- Stouffer one-sided p: ", signif(global[test_name == "stouffer_one_sided_p", value][1], 4), "\n\n", sep = "")
cat("Highlighted metrics:\n")
cat("- candidate proportion: beta=", signif(get_val("Score-defined cells: candidate proportion", "lm_beta_asd"), 4),
    ", nominal p=", signif(get_val("Score-defined cells: candidate proportion", "nominal_p"), 4), "\n", sep = "")
cat("- candidate-state score: beta=", signif(get_val("All microglia: candidate-state score", "lm_beta_asd"), 4),
    ", nominal p=", signif(get_val("All microglia: candidate-state score", "nominal_p"), 4), "\n", sep = "")
cat("- candidate-minus-reference: beta=", signif(get_val("All microglia: candidate - reference", "lm_beta_asd"), 4),
    ", nominal p=", signif(get_val("All microglia: candidate - reference", "nominal_p"), 4), "\n", sep = "")
cat("- log2(candidate/reference): beta=", signif(get_val("Score-defined cells: log2(candidate/reference)", "lm_beta_asd"), 4),
    ", nominal p=", signif(get_val("Score-defined cells: log2(candidate/reference)", "nominal_p"), 4), "\n\n", sep = "")
cat("Recommended manuscript sentence:\n")
cat("Across the strict shared validation set, all seven pre-specified validation metrics were directionally concordant with the discovery-prioritized expectation, with the strongest nominal support observed for the score-defined candidate proportion, while the overall effects remained modest and were not robust after multiple-testing correction.\n")
sink()

message("Package 6 done.")
