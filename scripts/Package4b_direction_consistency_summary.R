#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

parse_args <- function(args) {
  out <- list(
    pkg4_dir = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Package4_CrossCohortValidation_Velmeshev_v2",
    outdir   = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Package4b_DirectionConsistencySummary"
  )
  if (length(args) == 0) return(out)
  i <- 1
  while (i <= length(args)) {
    key <- args[i]
    if (startsWith(key, "--")) {
      key2 <- sub("^--", "", key)
      if (grepl("=", key2, fixed = TRUE)) {
        sp <- strsplit(key2, "=", fixed = TRUE)[[1]]
        out[[sp[1]]] <- sp[2]
      } else if (i < length(args)) {
        out[[key2]] <- args[i + 1]
        i <- i + 1
      }
    }
    i <- i + 1
  }
  out
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
pkg4_dir <- normalizePath(args$pkg4_dir, mustWork = FALSE)
outdir   <- args$outdir
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

message("[Package4b] pkg4_dir = ", pkg4_dir)
message("[Package4b] outdir   = ", outdir)

metric_config <- data.table(
  metric_id = c(
    "candidate_minus_reference",
    "candidate_state_score",
    "reference_like_score",
    "score_defined_candidate_proportion",
    "score_defined_reference_proportion",
    "log2_candidate_reference_ratio",
    "score_defined_candidate_vs_reference_shift"
  ),
  metric_label = c(
    "All microglia: candidate - reference",
    "All microglia: candidate-state score",
    "All microglia: reference-like score",
    "Score-defined cells: candidate proportion",
    "Score-defined cells: reference proportion",
    "Score-defined cells: log2(candidate/reference)",
    "Score-defined cells: candidate vs reference shift"
  ),
  expected_direction = c(
    "ASD > Control",
    "ASD > Control",
    "ASD < Control",
    "ASD > Control",
    "ASD < Control",
    "ASD > Control",
    "ASD > Control"
  )
)
metric_config[, expected_sign := ifelse(expected_direction == "ASD > Control", 1, -1)]

norm_name <- function(x) {
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

match_first_pattern <- function(x, patterns) {
  idx <- which(vapply(patterns, function(p) any(grepl(p, x, perl = TRUE)), logical(1)))
  if (length(idx) == 0) return(NA_character_)
  p <- patterns[idx[1]]
  hit <- which(grepl(p, x, perl = TRUE))
  if (length(hit) == 0) return(NA_character_)
  x[hit[1]]
}

find_metric_cols <- function(nm_cols) {
  pattern_list <- list(
    candidate_minus_reference = c(
      "^candidate_minus_reference$",
      "^all_microglia_candidate_minus_reference$",
      "^candidate_reference_difference$",
      "^candidate_reference_diff$",
      "^cand_minus_ref$",
      "^candidate_minus_ref$"
    ),
    candidate_state_score = c(
      "^candidate_state_score$",
      "^all_microglia_candidate_state_score$",
      "^candidate_score$",
      "^cand_state_score$",
      "^candidate_like_score$"
    ),
    reference_like_score = c(
      "^reference_like_score$",
      "^all_microglia_reference_like_score$",
      "^reference_score$",
      "^ref_like_score$"
    ),
    score_defined_candidate_proportion = c(
      "^score_defined_candidate_proportion$",
      "^score_defined_candidate_prop$",
      "^candidate_proportion$",
      "^candidate_prop$",
      "^prop_candidate$",
      "^fraction_candidate$"
    ),
    score_defined_reference_proportion = c(
      "^score_defined_reference_proportion$",
      "^score_defined_reference_prop$",
      "^reference_proportion$",
      "^reference_prop$",
      "^prop_reference$",
      "^fraction_reference$"
    ),
    log2_candidate_reference_ratio = c(
      "^log2_candidate_reference_ratio$",
      "^log2_candidate_reference$",
      "^log2_cand_ref_ratio$",
      "^log2_cand_ref$",
      "^log2_candidate_reference_proportion_ratio$"
    ),
    score_defined_candidate_vs_reference_shift = c(
      "^score_defined_candidate_vs_reference_shift$",
      "^candidate_vs_reference_shift$",
      "^score_defined_shift$",
      "^score_defined_candidate_reference_shift$",
      "^shift_candidate_vs_reference$"
    )
  )
  out <- lapply(pattern_list, function(pats) match_first_pattern(nm_cols, pats))
  out[!vapply(out, function(x) is.na(x) || length(x) == 0, logical(1))]
}

recode_dx <- function(x) {
  y <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(y))
  out[grepl("(^|_)(asd|case|patient|proband)(_|$)", gsub("[^a-z0-9]+", "_", y), perl = TRUE)] <- "ASD"
  out[grepl("(^|_)(ctrl|control|td|typical|non_asd)(_|$)", gsub("[^a-z0-9]+", "_", y), perl = TRUE)] <- "Control"
  out
}

safe_read <- function(file, nrows = Inf, select = NULL) {
  ext <- tolower(file)
  dt <- tryCatch({
    if (grepl("\\.csv(\\.gz)?$", ext, perl = TRUE)) {
      fread(file, nrows = nrows, select = select, showProgress = FALSE)
    } else {
      fread(file, nrows = nrows, select = select, sep = "\t", showProgress = FALSE)
    }
  }, error = function(e) NULL)
  dt
}

candidate_files <- list.files(
  pkg4_dir,
  recursive = TRUE,
  full.names = TRUE,
  pattern = "\\.(tsv|txt|csv)(\\.gz)?$"
)

if (length(candidate_files) == 0) {
  stop("No .tsv/.txt/.csv files found under: ", pkg4_dir)
}

message("[Package4b] candidate files found: ", length(candidate_files))

file_scan <- list()
metric_tables <- list()
metric_source_rows <- list()

for (f in candidate_files) {
  head_dt <- safe_read(f, nrows = 50)
  if (is.null(head_dt) || ncol(head_dt) == 0) next

  orig_cols <- colnames(head_dt)
  nm_cols <- norm_name(orig_cols)
  col_map <- setNames(orig_cols, nm_cols)

  dx_nm <- match_first_pattern(nm_cols, c("^diagnosis$", "^dx$", "^group$", "^status$", "^condition$"))
  sample_nm <- match_first_pattern(nm_cols, c(
    "^sample$", "^sample_id$", "^donor$", "^donor_id$", "^subject$", "^subject_id$",
    "^individual$", "^individual_id$", "^iid$", "^person$", "^orig_ident$"
  ))
  metric_hits_nm <- find_metric_cols(nm_cols)

  file_scan[[length(file_scan) + 1]] <- data.table(
    file = f,
    has_dx = !is.na(dx_nm),
    has_sample = !is.na(sample_nm),
    metric_hits = paste(names(metric_hits_nm), collapse = ";")
  )

  if (is.na(dx_nm) || length(metric_hits_nm) == 0) next

  need_nm <- unique(c(dx_nm, sample_nm, unlist(metric_hits_nm, use.names = FALSE)))
  need_orig <- unname(col_map[need_nm])
  dt <- safe_read(f, select = need_orig)
  if (is.null(dt) || nrow(dt) == 0) next

  setDT(dt)
  rename_vec <- setNames(names(col_map)[match(colnames(dt), col_map)], colnames(dt))
  setnames(dt, old = names(rename_vec), new = unname(rename_vec))

  dt[, .pkg4b_dx := recode_dx(get(dx_nm))]
  dt <- dt[!is.na(.pkg4b_dx)]
  if (nrow(dt) == 0) next

  if (!is.na(sample_nm)) {
    dt[, .pkg4b_sample := as.character(get(sample_nm))]
  } else {
    dt[, .pkg4b_sample := paste0("row_", seq_len(.N))]
  }

  for (metric_id in names(metric_hits_nm)) {
    metric_nm <- metric_hits_nm[[metric_id]]
    sub <- dt[, .(.pkg4b_sample, .pkg4b_dx, value = suppressWarnings(as.numeric(get(metric_nm))))]
    sub <- sub[!is.na(value)]
    if (nrow(sub) == 0) next

    sub <- sub[, .(value = mean(value, na.rm = TRUE)), by = .(.pkg4b_sample, .pkg4b_dx)]
    n_complete <- nrow(sub)
    n_samples <- uniqueN(sub$.pkg4b_sample)
    score <- n_complete + 10000 * as.integer(!is.na(sample_nm)) + 100 * as.integer(n_samples >= 4)

    metric_source_rows[[length(metric_source_rows) + 1]] <- data.table(
      metric_id = metric_id,
      source_file = f,
      source_column = metric_nm,
      n_complete = n_complete,
      n_samples = n_samples,
      score = score
    )

    current <- metric_tables[[metric_id]]
    if (is.null(current) || score > attr(current, "score")) {
      metric_tables[[metric_id]] <- copy(sub)
      attr(metric_tables[[metric_id]], "score") <- score
      attr(metric_tables[[metric_id]], "source_file") <- f
      attr(metric_tables[[metric_id]], "source_column") <- metric_nm
    }
  }
}

scan_dt <- rbindlist(file_scan, fill = TRUE)
source_candidates_dt <- rbindlist(metric_source_rows, fill = TRUE)
if (nrow(scan_dt) > 0) fwrite(scan_dt, file.path(outdir, "00_detected_files_scan.tsv"), sep = "\t")
if (nrow(source_candidates_dt) > 0) fwrite(source_candidates_dt, file.path(outdir, "01_metric_source_candidates.tsv"), sep = "\t")

summarize_metric <- function(dt, metric_id) {
  dt <- copy(dt)
  setnames(dt, c(".pkg4b_sample", ".pkg4b_dx"), c("sample", "dx"))
  dt <- dt[dx %in% c("ASD", "Control") & !is.na(value)]
  if (nrow(dt) == 0) {
    return(data.table(
      metric_id = metric_id,
      n_asd = NA_integer_,
      n_control = NA_integer_,
      mean_asd = NA_real_,
      mean_control = NA_real_,
      median_asd = NA_real_,
      median_control = NA_real_,
      mean_diff_asd_minus_control = NA_real_,
      median_diff_asd_minus_control = NA_real_,
      lm_beta_asd = NA_real_,
      lm_p = NA_real_,
      wilcox_p = NA_real_,
      observed_direction = NA_character_,
      primary_effect = NA_real_
    ))
  }

  dx_levels <- unique(dt$dx)
  if (!all(c("ASD", "Control") %in% dx_levels)) {
    return(data.table(
      metric_id = metric_id,
      n_asd = sum(dt$dx == "ASD"),
      n_control = sum(dt$dx == "Control"),
      mean_asd = NA_real_,
      mean_control = NA_real_,
      median_asd = NA_real_,
      median_control = NA_real_,
      mean_diff_asd_minus_control = NA_real_,
      median_diff_asd_minus_control = NA_real_,
      lm_beta_asd = NA_real_,
      lm_p = NA_real_,
      wilcox_p = NA_real_,
      observed_direction = NA_character_,
      primary_effect = NA_real_
    ))
  }

  asd_vals <- dt[dx == "ASD", value]
  ctrl_vals <- dt[dx == "Control", value]

  out <- data.table(
    metric_id = metric_id,
    n_asd = length(asd_vals),
    n_control = length(ctrl_vals),
    mean_asd = mean(asd_vals, na.rm = TRUE),
    mean_control = mean(ctrl_vals, na.rm = TRUE),
    median_asd = median(asd_vals, na.rm = TRUE),
    median_control = median(ctrl_vals, na.rm = TRUE)
  )
  out[, mean_diff_asd_minus_control := mean_asd - mean_control]
  out[, median_diff_asd_minus_control := median_asd - median_control]

  dt[, dx := factor(dx, levels = c("Control", "ASD"))]

  fit <- tryCatch(lm(value ~ dx, data = dt), error = function(e) NULL)
  if (!is.null(fit)) {
    cf <- summary(fit)$coefficients
    if ("dxASD" %in% rownames(cf)) {
      out[, lm_beta_asd := unname(cf["dxASD", "Estimate"])]
      out[, lm_p := unname(cf["dxASD", "Pr(>|t|)"])]
    } else {
      out[, `:=`(lm_beta_asd = NA_real_, lm_p = NA_real_)]
    }
  } else {
    out[, `:=`(lm_beta_asd = NA_real_, lm_p = NA_real_)]
  }

  out[, wilcox_p := tryCatch(wilcox.test(value ~ dx, data = dt, exact = FALSE)$p.value, error = function(e) NA_real_)]

  primary_effect <- out$lm_beta_asd
  if (is.na(primary_effect)) primary_effect <- out$mean_diff_asd_minus_control
  out[, primary_effect := primary_effect]

  out[, observed_direction := fifelse(
    is.na(primary_effect),
    NA_character_,
    fifelse(primary_effect > 0, "ASD > Control", fifelse(primary_effect < 0, "ASD < Control", "No difference"))
  )]
  out
}

summary_list <- list()
source_best_list <- list()

for (mid in metric_config$metric_id) {
  dt <- metric_tables[[mid]]
  if (is.null(dt)) {
    summary_list[[length(summary_list) + 1]] <- data.table(
      metric_id = mid,
      n_asd = NA_integer_,
      n_control = NA_integer_,
      mean_asd = NA_real_,
      mean_control = NA_real_,
      median_asd = NA_real_,
      median_control = NA_real_,
      mean_diff_asd_minus_control = NA_real_,
      median_diff_asd_minus_control = NA_real_,
      lm_beta_asd = NA_real_,
      lm_p = NA_real_,
      wilcox_p = NA_real_,
      observed_direction = NA_character_,
      primary_effect = NA_real_,
      source_file = NA_character_,
      source_column = NA_character_
    )
  } else {
    sum_dt <- summarize_metric(dt, mid)
    sum_dt[, source_file := attr(dt, "source_file")]
    sum_dt[, source_column := attr(dt, "source_column")]
    summary_list[[length(summary_list) + 1]] <- sum_dt
  }
}

summary_dt <- rbindlist(summary_list, fill = TRUE)
summary_dt <- merge(metric_config, summary_dt, by = "metric_id", all.x = TRUE, sort = FALSE)
summary_dt[, direction_match := fifelse(
  is.na(primary_effect),
  NA,
  sign(primary_effect) == expected_sign
)]
summary_dt[, aligned_effect := primary_effect * expected_sign]
summary_dt[, wilcox_fdr := p.adjust(wilcox_p, method = "BH")]
summary_dt[, lm_fdr := p.adjust(lm_p, method = "BH")]
summary_dt[, evidence_note := fifelse(
  is.na(primary_effect),
  "metric_not_found",
  fifelse(direction_match & !is.na(wilcox_p) & wilcox_p < 0.05,
          "direction_match_and_wilcox_p_lt_0.05",
          fifelse(direction_match,
                  "direction_match_but_not_significant",
                  "direction_mismatch"))
)]
summary_dt[, source_file_basename := fifelse(is.na(source_file), NA_character_, basename(source_file))]

pretty_dt <- copy(summary_dt)[, .(
  metric_label,
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
  mean_asd,
  mean_control,
  median_asd,
  median_control,
  aligned_effect,
  source_file_basename,
  source_column,
  evidence_note
)]

fwrite(summary_dt, file.path(outdir, "02_direction_consistency_summary_full.tsv"), sep = "\t")
fwrite(pretty_dt, file.path(outdir, "03_direction_consistency_summary_pretty.tsv"), sep = "\t")

# manuscript-style compact table
compact_dt <- copy(summary_dt)
compact_dt[, observed_short := fifelse(observed_direction == "ASD > Control", ">", fifelse(observed_direction == "ASD < Control", "<", "="))]
compact_dt[, expected_short := fifelse(expected_direction == "ASD > Control", ">", "<")]
compact_dt[, direction_call := fifelse(is.na(direction_match), "missing", fifelse(direction_match, "match", "mismatch"))]
compact_dt <- compact_dt[, .(
  metric = metric_label,
  expected = expected_short,
  observed = observed_short,
  direction_call,
  beta_asd = signif(lm_beta_asd, 4),
  wilcox_p = signif(wilcox_p, 4),
  n_asd,
  n_control
)]
fwrite(compact_dt, file.path(outdir, "04_direction_consistency_summary_manuscript.tsv"), sep = "\t")

plot_dt <- summary_dt[!is.na(aligned_effect)]
if (nrow(plot_dt) > 0) {
  plot_dt[, metric_label := factor(metric_label, levels = rev(metric_config$metric_label))]
  plot_dt[, sig_class := fifelse(!is.na(wilcox_p) & wilcox_p < 0.05, "Wilcox p < 0.05", "Wilcox p >= 0.05")]
  plot_dt[, match_class := fifelse(direction_match, "Direction matched", "Direction mismatched")]
  plot_dt[, label_text := fifelse(
    is.na(wilcox_p),
    observed_direction,
    sprintf("%s | p=%.3g", observed_direction, wilcox_p)
  )]

  xmax <- max(abs(plot_dt$aligned_effect), na.rm = TRUE)
  if (!is.finite(xmax) || xmax == 0) xmax <- 1
  xpad <- xmax * 0.45

  p <- ggplot(plot_dt, aes(x = aligned_effect, y = metric_label)) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.5, colour = "grey40") +
    geom_segment(aes(x = 0, xend = aligned_effect, yend = metric_label), linewidth = 0.7, colour = "grey60") +
    geom_point(aes(shape = sig_class, fill = match_class), size = 3.6, colour = "black", stroke = 0.4) +
    geom_text(
      aes(label = label_text, hjust = ifelse(aligned_effect >= 0, -0.05, 1.05)),
      size = 3.2
    ) +
    scale_shape_manual(values = c("Wilcox p < 0.05" = 21, "Wilcox p >= 0.05" = 24)) +
    scale_fill_manual(values = c("Direction matched" = "#4DAF4A", "Direction mismatched" = "#E41A1C")) +
    coord_cartesian(xlim = c(-xmax - xpad * 0.2, xmax + xpad), clip = "off") +
    labs(
      x = "Expected-direction aligned effect (positive = concordant with discovery-prioritized expectation)",
      y = NULL,
      title = "Package 4b: score-based validation direction-consistency summary",
      subtitle = "Velmeshev validation microglia, sample-aware donor-level metrics"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.title = element_blank(),
      plot.margin = margin(10, 90, 10, 10)
    )

  ggsave(file.path(outdir, "05_direction_consistency_summary_plot.pdf"), p, width = 11, height = 4.8)
  ggsave(file.path(outdir, "05_direction_consistency_summary_plot.png"), p, width = 11, height = 4.8, dpi = 320)
}

sink(file.path(outdir, "06_package4b_run_notes.txt"))
cat("Package 4b run notes\n")
cat("=====================\n")
cat("pkg4_dir: ", pkg4_dir, "\n", sep = "")
cat("outdir:   ", outdir, "\n\n", sep = "")
cat("Metrics requested:\n")
for (i in seq_len(nrow(metric_config))) {
  cat("- ", metric_config$metric_id[i], " | expected: ", metric_config$expected_direction[i], "\n", sep = "")
}
cat("\nDetected best sources:\n")
for (i in seq_len(nrow(summary_dt))) {
  cat("- ", summary_dt$metric_id[i], ": ",
      ifelse(is.na(summary_dt$source_file_basename[i]), "NOT FOUND", summary_dt$source_file_basename[i]),
      " | col=", ifelse(is.na(summary_dt$source_column[i]), "NA", summary_dt$source_column[i]), "\n", sep = "")
}
cat("\nInterpretation rule:\n")
cat("aligned_effect > 0 means the observed effect direction is concordant with the pre-specified expected direction.\n")
cat("direction_match is based on the sign of lm_beta_asd when available; otherwise the mean difference ASD-Control is used.\n")
cat("wilcox_p is provided as a non-parametric companion test.\n")
sink()

message("[Package4b] done.")
