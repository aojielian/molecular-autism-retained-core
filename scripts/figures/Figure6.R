#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(stringr)
  library(forcats)
  library(scales)
  library(optparse)
})


project_dir_default <- Sys.getenv("MA_PROJECT_DIR", ".")
option_list <- list(
  make_option("--pkg8", type = "character", default = file.path(project_dir_default, "results", "Package8_ExternalValidation_Gandal2022"), help = "Package8 directory [default: %default]"),
  make_option("--pkg8b", type = "character", default = file.path(project_dir_default, "results", "Package8b_Gandal2022_LOO_Sensitivity"), help = "Package8b directory [default: %default]"),
  make_option("--outdir", type = "character", default = file.path(project_dir_default, "results", "Figure6_Gandal2022_ExternalValidation"), help = "Output directory [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
cfg <- list(
  pkg8 = opt$pkg8,
  pkg8b = opt$pkg8b,
  outdir = opt$outdir
)

dir.create(cfg$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cfg$outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cfg$outdir, "plots"), recursive = TRUE, showWarnings = FALSE)

read_tsv_auto <- function(path) {
  if (grepl("\\.gz$", path)) {
    fread(cmd = paste("gzip -dc", shQuote(path)), sep = "\t", header = TRUE, data.table = TRUE)
  } else {
    fread(path, sep = "\t", header = TRUE, data.table = TRUE)
  }
}

find_col <- function(df, candidates, required = TRUE) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit) > 0) return(hit[1])
  if (required) stop("Cannot find column. Candidates tried: ", paste(candidates, collapse = ", "),
                     "\nAvailable: ", paste(colnames(df), collapse = ", "))
  NULL
}

std_dx <- function(x) {
  x <- as.character(x)
  dplyr::case_when(
    x %in% c("ASD", "Autism", "Case", "case", "1") ~ "ASD",
    x %in% c("Control", "CTL", "ctrl", "control", "0") ~ "Control",
    TRUE ~ x
  )
}

fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
}

fmt_p <- function(x) {
  out <- rep("NA", length(x))
  ok <- is.finite(x)
  out[ok & x < 1e-4] <- formatC(x[ok & x < 1e-4], format = "e", digits = 2)
  out[ok & x >= 1e-4] <- formatC(x[ok & x >= 1e-4], format = "f", digits = 3)
  out
}

clean_program_name <- function(x) {
  x2 <- tolower(gsub("[^a-z0-9]+", "_", x))
  dplyr::case_when(
    str_detect(x2, "candidate") & str_detect(x2, "minus") ~ "candidate_minus_reference",
    str_detect(x2, "candidate") ~ "candidate_program",
    str_detect(x2, "reference") ~ "reference_program",
    TRUE ~ x2
  )
}

program_pretty <- c(
  candidate_program = "Candidate-program score",
  reference_program = "Reference-program score",
  candidate_minus_reference = "Candidate-minus-reference score"
)

program_order <- c("candidate_program", "reference_program", "candidate_minus_reference")
program_colors <- c("Control" = "#CFCFCF", "ASD" = "#D8AEB6")

# ------------------------------
# Read main sample-level scores
# ------------------------------
scores <- read_tsv_auto(file.path(cfg$pkg8, "07_sample_level_program_scores.tsv"))
results <- read_tsv_auto(file.path(cfg$pkg8, "08_program_mixed_model_results.tsv"))
loo <- read_tsv_auto(file.path(cfg$pkg8b, "04_leave_one_region_out_results.tsv"))
loo_vs_base <- read_tsv_auto(file.path(cfg$pkg8b, "05_leave_one_region_out_vs_baseline.tsv"))
loo_sum <- read_tsv_auto(file.path(cfg$pkg8b, "06_leave_one_region_out_stability_summary.tsv"))

fwrite(scores, file.path(cfg$outdir, "tables", "07_sample_level_program_scores.copy.tsv"), sep = "\t")
fwrite(results, file.path(cfg$outdir, "tables", "08_program_mixed_model_results.copy.tsv"), sep = "\t")
fwrite(loo, file.path(cfg$outdir, "tables", "04_leave_one_region_out_results.copy.tsv"), sep = "\t")
fwrite(loo_vs_base, file.path(cfg$outdir, "tables", "05_leave_one_region_out_vs_baseline.copy.tsv"), sep = "\t")
fwrite(loo_sum, file.path(cfg$outdir, "tables", "06_leave_one_region_out_stability_summary.copy.tsv"), sep = "\t")

message("[read] scores columns: ", paste(colnames(scores), collapse = ", "))
message("[read] results columns: ", paste(colnames(results), collapse = ", "))
message("[read] loo columns: ", paste(colnames(loo), collapse = ", "))

# ------------------------------
# Standardize sample-level table
# ------------------------------
dx_col <- find_col(scores, c("dx_std", "dx", "Diagnosis", "diagnosis", "group", "Group", "Dx"))
region_col <- find_col(scores, c("region", "Region", "cortical_region", "brain_region", "BrainRegion"), required = FALSE)
id_col <- find_col(scores, c("sample_id", "SampleID", "sample", "donor", "subject_id", "IID", "id"), required = FALSE)

scores[[dx_col]] <- std_dx(scores[[dx_col]])
setnames(scores, dx_col, "dx_std")
if (!is.null(region_col)) setnames(scores, region_col, "region_std")
if (!is.null(id_col)) setnames(scores, id_col, "sample_std")

# identify score columns
score_cols <- colnames(scores)[vapply(scores, is.numeric, logical(1))]
# exclude obvious non-score columns
score_cols <- setdiff(score_cols, c("age", "RIN", "PMI", "sex", "n", "N"))

score_map <- data.frame(
  col = score_cols,
  clean = clean_program_name(score_cols),
  stringsAsFactors = FALSE
)
score_map <- score_map[score_map$clean %in% program_order, , drop = FALSE]
score_map <- score_map[!duplicated(score_map$clean), , drop = FALSE]

if (nrow(score_map) < 3) {
  stop("Could not identify the three program score columns from sample-level table. Available numeric columns: ",
       paste(score_cols, collapse = ", "))
}

scores_long <- rbindlist(lapply(seq_len(nrow(score_map)), function(i) {
  data.table(
    dx = scores$dx_std,
    region = if ("region_std" %in% colnames(scores)) scores$region_std else NA_character_,
    sample_id = if ("sample_std" %in% colnames(scores)) scores$sample_std else paste0("sample_", seq_len(nrow(scores))),
    program = score_map$clean[i],
    score = scores[[score_map$col[i]]]
  )
}), use.names = TRUE)

scores_long[, dx := factor(dx, levels = c("Control", "ASD"))]
scores_long[, program := factor(program, levels = program_order)]

# ------------------------------
# Standardize mixed model results
# ------------------------------
res_program_col <- find_col(results, c("program", "metric", "score", "outcome", "analysis", "program_name"))
res_beta_col <- find_col(results, c("beta_ASD_vs_Control", "beta", "estimate", "coef", "effect", "dx_beta", "Estimate"))
res_p_col <- find_col(results, c("p_value", "p", "pval", "PValue", "Pr(>|t|)", "pvalue"))
res_fdr_col <- find_col(results, c("BH_FDR", "fdr", "FDR", "padj", "adj_p", "q_value", "qvalue"), required = FALSE)

results_std <- copy(results)
results_std[, program_clean := clean_program_name(get(res_program_col))]
results_std <- results_std[program_clean %in% program_order]
results_std <- results_std[!duplicated(program_clean)]
results_std[, beta_std := as.numeric(get(res_beta_col))]
results_std[, p_std := as.numeric(get(res_p_col))]
results_std[, fdr_std := if (!is.null(res_fdr_col)) as.numeric(get(res_fdr_col)) else NA_real_]

# ------------------------------
# Standardize leave-one-region-out results
# ------------------------------
loo_program_col <- find_col(loo, c("program", "metric", "score", "outcome", "analysis", "program_name"))
loo_beta_col <- find_col(loo, c("beta_ASD_vs_Control", "beta", "estimate", "coef", "effect", "dx_beta", "Estimate"))
loo_region_col <- find_col(loo, c("omitted_region", "left_out_region", "leave_out_region", "excluded_region", "region", "held_out_region"))
loo_p_col <- find_col(loo, c("p_value", "p", "pval", "PValue", "pvalue"), required = FALSE)
loo_fdr_col <- find_col(loo, c("p_value_BH_within_all", "p_value_BH_within_LOO", "fdr", "FDR", "padj", "adj_p", "q_value"), required = FALSE)

loo_std <- copy(loo)
loo_std[, program_clean := clean_program_name(get(loo_program_col))]
loo_std <- loo_std[program_clean %in% program_order]
loo_std[, beta_std := as.numeric(get(loo_beta_col))]
loo_std[, region_left_out := as.character(get(loo_region_col))]
loo_std[, p_std := if (!is.null(loo_p_col)) as.numeric(get(loo_p_col)) else NA_real_]
loo_std[, fdr_std := if (!is.null(loo_fdr_col)) as.numeric(get(loo_fdr_col)) else NA_real_]

# get baseline beta from main results
base_beta <- results_std[, .(program_clean, base_beta = beta_std)]
loo_std <- merge(loo_std, base_beta, by = "program_clean", all.x = TRUE)

# ------------------------------
# Main score panels A-C
# ------------------------------
make_main_panel <- function(program_key, title_text) {
  dt <- scores_long[program == program_key]
  rr <- results_std[program_clean == program_key]
  ann <- paste0(
    "beta = ", fmt_num(rr$beta_std, 3), "\n",
    "P = ", fmt_p(rr$p_std),
    if (!all(is.na(rr$fdr_std))) paste0("\nFDR = ", fmt_p(rr$fdr_std)) else ""
  )

  y_top <- quantile(dt$score, 0.97, na.rm = TRUE)

  ggplot(dt, aes(x = dx, y = score, fill = dx)) +
    geom_violin(trim = FALSE, color = NA, alpha = 0.9) +
    geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white", color = "black") +
    geom_jitter(width = 0.08, size = 2.1, alpha = 0.95) +
    scale_fill_manual(values = program_colors) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      plot.title = element_text(face = "bold")
    ) +
    labs(
      title = title_text,
      x = NULL,
      y = unname(program_pretty[program_key])
    ) +
    annotate("label", x = 1.5, y = y_top, label = ann, size = 3.3)
}

pA <- make_main_panel("candidate_program", "A. Candidate-program score")
pB <- make_main_panel("reference_program", "B. Reference-program score")
pC <- make_main_panel("candidate_minus_reference", "C. Candidate-minus-reference score")

# ------------------------------
# LOO forest/lollipop panels D-F
# ------------------------------
make_loo_panel <- function(program_key, title_text) {
  dt <- loo_std[program_clean == program_key]
  dt <- dt[!(toupper(trimws(region_left_out)) %in% c("NONE", "BASELINE", "ALL"))]
  dt <- dt[order(beta_std)]
  dt[, region_label := factor(region_left_out, levels = region_left_out)]
  base <- unique(dt$base_beta)[1]

  ggplot(dt, aes(x = beta_std, y = region_label)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
    geom_vline(xintercept = base, linetype = 3, color = "grey55") +
    geom_segment(aes(x = 0, xend = beta_std, y = region_label, yend = region_label),
                 color = "grey70", linewidth = 0.7) +
    geom_point(size = 2.8, fill = "white", shape = 21, color = "black") +
    theme_bw(base_size = 11) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold")
    ) +
    labs(
      title = title_text,
      x = "Beta (ASD vs Control)",
      y = NULL
    )
}

pD <- make_loo_panel("candidate_program", "D. Candidate-program leave-one-region-out")
pE <- make_loo_panel("reference_program", "E. Reference-program leave-one-region-out")
pF <- make_loo_panel("candidate_minus_reference", "F. Candidate-minus-reference leave-one-region-out")

# ------------------------------
# Assemble
# ------------------------------
fig6 <- (pA | pB | pC) / (pD | pE | pF) +
  plot_layout(heights = c(1, 1))

ggsave(file.path(cfg$outdir, "plots", "Figure6_combined.pdf"), fig6, width = 18, height = 12)
ggsave(file.path(cfg$outdir, "plots", "Figure6_combined.png"), fig6, width = 18, height = 12, dpi = 300)

saveRDS(list(score_map = score_map,
             results_std = results_std,
             loo_std = loo_std),
        file.path(cfg$outdir, "tables", "Figure6_objects.rds"))

message("[done] Figure 6 written to: ", cfg$outdir)
