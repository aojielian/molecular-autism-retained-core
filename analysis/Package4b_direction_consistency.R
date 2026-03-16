#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
})

DEFAULT_PROJECT_DIR <- normalizePath(Sys.getenv("MA_PROJECT_DIR", "."), mustWork = FALSE)
option_list <- list(
  make_option("--indir", type = "character", default = file.path(DEFAULT_PROJECT_DIR, "results", "Package4_CrossCohortValidation_Velmeshev_v2")),
  make_option("--outdir", type = "character", default = file.path(DEFAULT_PROJECT_DIR, "results", "Package4b_DirectionConsistencySummary_v2_manual"))
)
opt <- parse_args(OptionParser(option_list = option_list))
INDIR  <- opt$indir
OUTDIR <- opt$outdir
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

f_cell  <- file.path(INDIR, "12_validation_cell_scores.tsv.gz")
f_prop  <- file.path(INDIR, "15_validation_sample_score_defined_proportions.tsv")
f_shift <- file.path(INDIR, "16_validation_sample_score_defined_shift.tsv")

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
  if (!file.exists(f)) stop("File not found: ", f)
  if (grepl("\\.csv(\\.gz)?$", f, perl = TRUE)) {
    fread(f)
  } else {
    fread(f, sep = "\t")
  }
}

standardize_sample_dx <- function(dt, file_tag = "") {
  setDT(dt)

  old0 <- colnames(dt)

  # 第一步：先把原始重复列名变唯一，避免 setnames(old=...) 因重复报错
  old_unique <- make.unique(old0, sep = "_origdup")
  setnames(dt, old = seq_along(old0), new = old_unique)

  # 第二步：规范化列名；如果规范化后再次重复，再次强制唯一
  new0 <- norm_name(old_unique)
  new_unique <- make.unique(new0, sep = "_dup")
  setnames(dt, old = seq_along(new_unique), new = new_unique)

  dx_col <- find_col(colnames(dt), c(
    "^diagnosis$",
    "^dx$",
    "^group$",
    "^status$",
    "^condition$"
  ))

  sample_col <- find_col(colnames(dt), c(
    "^sample$",
    "^sample_id$",
    "^sampleid$",
    "^donor$",
    "^donor_id$",
    "^donorid$",
    "^subject$",
    "^subject_id$",
    "^subjectid$",
    "^individual$",
    "^individual_id$",
    "^iid$",
    "^orig_ident$"
  ))

  if (is.na(dx_col)) {
    stop(
      "No diagnosis column found in ", file_tag,
      ". Available columns: ", paste(colnames(dt), collapse = ", ")
    )
  }

  if (is.na(sample_col)) {
    stop(
      "No sample column found in ", file_tag,
      ". Available columns: ", paste(colnames(dt), collapse = ", ")
    )
  }

  dt[, dx2 := recode_dx(get(dx_col))]
  dt <- dt[!is.na(dx2)]
  dt[, sample2 := as.character(get(sample_col))]

  list(dt = dt, dx_col = dx_col, sample_col = sample_col)
}

summarize_metric <- function(dt, metric_col, metric_label, expected_direction, source_file) {
  if (is.na(metric_col) || !(metric_col %in% colnames(dt))) {
    return(data.table(
      metric_label = metric_label,
      expected_direction = expected_direction,
      observed_direction = NA_character_,
      direction_match = NA,
      n_asd = NA_integer_,
      n_control = NA_integer_,
      lm_beta_asd = NA_real_,
      lm_p = NA_real_,
      wilcox_p = NA_real_,
      mean_asd = NA_real_,
      mean_control = NA_real_,
      median_asd = NA_real_,
      median_control = NA_real_,
      aligned_effect = NA_real_,
      source_file_basename = basename(source_file),
      source_column = ifelse(is.na(metric_col), NA_character_, metric_col),
      evidence_note = "metric_not_found"
    ))
  }

  sub <- copy(dt)[, .(sample2, dx2, value = suppressWarnings(as.numeric(get(metric_col))))]
  sub <- sub[!is.na(value)]

  if (nrow(sub) == 0) {
    return(data.table(
      metric_label = metric_label,
      expected_direction = expected_direction,
      observed_direction = NA_character_,
      direction_match = NA,
      n_asd = NA_integer_,
      n_control = NA_integer_,
      lm_beta_asd = NA_real_,
      lm_p = NA_real_,
      wilcox_p = NA_real_,
      mean_asd = NA_real_,
      mean_control = NA_real_,
      median_asd = NA_real_,
      median_control = NA_real_,
      aligned_effect = NA_real_,
      source_file_basename = basename(source_file),
      source_column = metric_col,
      evidence_note = "metric_not_found"
    ))
  }

  # donor/sample 级聚合
  sub <- sub[, .(value = mean(value, na.rm = TRUE)), by = .(sample2, dx2)]

  asd_vals  <- sub[dx2 == "ASD", value]
  ctrl_vals <- sub[dx2 == "Control", value]

  mean_asd <- mean(asd_vals, na.rm = TRUE)
  mean_control <- mean(ctrl_vals, na.rm = TRUE)
  median_asd <- median(asd_vals, na.rm = TRUE)
  median_control <- median(ctrl_vals, na.rm = TRUE)

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

  p_w <- tryCatch(
    wilcox.test(value ~ dx2, data = sub, exact = FALSE)$p.value,
    error = function(e) NA_real_
  )

  primary_effect <- ifelse(is.na(beta), mean_asd - mean_control, beta)

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

  evidence_note <- if (is.na(primary_effect)) {
    "metric_not_found"
  } else if (direction_match && !is.na(p_w) && p_w < 0.05) {
    "direction_match_and_wilcox_p_lt_0.05"
  } else if (direction_match) {
    "direction_match_but_not_significant"
  } else {
    "direction_mismatch"
  }

  data.table(
    metric_label = metric_label,
    expected_direction = expected_direction,
    observed_direction = observed_direction,
    direction_match = direction_match,
    n_asd = length(asd_vals),
    n_control = length(ctrl_vals),
    lm_beta_asd = beta,
    lm_p = p_lm,
    wilcox_p = p_w,
    mean_asd = mean_asd,
    mean_control = mean_control,
    median_asd = median_asd,
    median_control = median_control,
    aligned_effect = aligned_effect,
    source_file_basename = basename(source_file),
    source_column = metric_col,
    evidence_note = evidence_note
  )
}

# -----------------------------
# 1) 读入三个核心文件
# -----------------------------
cell0  <- read_dt(f_cell)
prop0  <- read_dt(f_prop)
shift0 <- read_dt(f_shift)

cellx  <- standardize_sample_dx(cell0,  "12_validation_cell_scores")
propx  <- standardize_sample_dx(prop0,  "15_validation_sample_score_defined_proportions")
shiftx <- standardize_sample_dx(shift0, "16_validation_sample_score_defined_shift")

# 输出 debug 列名
writeLines(colnames(cellx$dt),  con = file.path(OUTDIR, "debug_colnames_12_validation_cell_scores.txt"))
writeLines(colnames(propx$dt),  con = file.path(OUTDIR, "debug_colnames_15_validation_sample_score_defined_proportions.txt"))
writeLines(colnames(shiftx$dt), con = file.path(OUTDIR, "debug_colnames_16_validation_sample_score_defined_shift.txt"))

cell  <- cellx$dt
prop  <- propx$dt
shift <- shiftx$dt

# -----------------------------
# 2) 统一 comparable sample set
#    先取 score-defined 两张表的交集，再与 cell_scores 取交集
# -----------------------------
shared_samples <- intersect(unique(prop$sample2), unique(shift$sample2))
shared_samples <- intersect(shared_samples, unique(cell$sample2))

if (length(shared_samples) == 0) {
  stop("No shared comparable samples found across 12/15/16 tables.")
}

prop  <- prop[sample2 %in% shared_samples]
shift <- shift[sample2 %in% shared_samples]
cell  <- cell[sample2 %in% shared_samples]

shared_dt <- rbindlist(list(
  unique(cell[,  .(sample2, dx2)]),
  unique(prop[,  .(sample2, dx2)]),
  unique(shift[, .(sample2, dx2)])
), fill = TRUE)
shared_dt <- unique(shared_dt)
setorder(shared_dt, dx2, sample2)
fwrite(shared_dt, file.path(OUTDIR, "00_shared_comparable_samples.tsv"), sep = "\t")

sample_count_dt <- shared_dt[, .N, by = dx2]
fwrite(sample_count_dt, file.path(OUTDIR, "00b_shared_sample_counts.tsv"), sep = "\t")

# -----------------------------
# 3) 识别列名
# -----------------------------
cell_cols  <- colnames(cell)
prop_cols  <- colnames(prop)
shift_cols <- colnames(shift)

col_candidate_minus_reference <- find_col(cell_cols, c(
  "^candidate_minus_reference$",
  "^all_microglia_candidate_minus_reference$"
))

col_candidate_state_score <- find_col(cell_cols, c(
  "^candidate_state_score$",
  "^all_microglia_candidate_state_score$",
  "^candidate_score$"
))

col_reference_like_score <- find_col(cell_cols, c(
  "^reference_like_score$",
  "^all_microglia_reference_like_score$",
  "^reference_score$"
))

col_prop_candidate <- find_col(prop_cols, c(
  "^prop_candidate$",
  "^candidate_proportion$",
  "^candidate_prop$",
  "^score_defined_candidate_proportion$"
))

col_prop_reference <- find_col(prop_cols, c(
  "^prop_reference$",
  "^reference_proportion$",
  "^reference_prop$",
  "^score_defined_reference_proportion$"
))

col_shift <- find_col(shift_cols, c(
  "^score_defined_shift$",
  "^candidate_vs_reference_shift$",
  "^score_defined_candidate_vs_reference_shift$"
))

metric_detect_dt <- data.table(
  metric = c(
    "candidate_minus_reference",
    "candidate_state_score",
    "reference_like_score",
    "prop_candidate",
    "prop_reference",
    "score_defined_shift"
  ),
  detected_column = c(
    col_candidate_minus_reference,
    col_candidate_state_score,
    col_reference_like_score,
    col_prop_candidate,
    col_prop_reference,
    col_shift
  )
)
fwrite(metric_detect_dt, file.path(OUTDIR, "01_detected_metric_columns.tsv"), sep = "\t")

# -----------------------------
# 4) 现算 log2(candidate/reference)
# -----------------------------
eps <- 1e-6
if (!is.na(col_prop_candidate) && !is.na(col_prop_reference)) {
  prop[, log2_candidate_reference_ratio := log2((get(col_prop_candidate) + eps) / (get(col_prop_reference) + eps))]
  col_log2_ratio <- "log2_candidate_reference_ratio"
} else {
  col_log2_ratio <- NA_character_
}

# -----------------------------
# 5) 汇总各指标
# -----------------------------
res_list <- list(
  summarize_metric(
    cell,
    col_candidate_minus_reference,
    "All microglia: candidate - reference",
    "ASD > Control",
    f_cell
  ),

  summarize_metric(
    cell,
    col_candidate_state_score,
    "All microglia: candidate-state score",
    "ASD > Control",
    f_cell
  ),

  summarize_metric(
    cell,
    col_reference_like_score,
    "All microglia: reference-like score",
    "ASD < Control",
    f_cell
  ),

  summarize_metric(
    prop,
    col_prop_candidate,
    "Score-defined cells: candidate proportion",
    "ASD > Control",
    f_prop
  ),

  summarize_metric(
    prop,
    col_prop_reference,
    "Score-defined cells: reference proportion",
    "ASD < Control",
    f_prop
  ),

  summarize_metric(
    prop,
    col_log2_ratio,
    "Score-defined cells: log2(candidate/reference)",
    "ASD > Control",
    f_prop
  ),

  summarize_metric(
    shift,
    col_shift,
    "Score-defined cells: candidate vs reference shift",
    "ASD > Control",
    f_shift
  )
)

res <- rbindlist(res_list, fill = TRUE)
res[, lm_fdr := p.adjust(lm_p, method = "BH")]
res[, wilcox_fdr := p.adjust(wilcox_p, method = "BH")]

pretty <- res[, .(
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

fwrite(pretty, file.path(OUTDIR, "03_direction_consistency_summary_pretty.tsv"), sep = "\t")

manu <- copy(pretty)
manu[, expected_short := ifelse(expected_direction == "ASD > Control", ">", "<")]
manu[, observed_short := ifelse(
  observed_direction == "ASD > Control", ">",
  ifelse(observed_direction == "ASD < Control", "<", "=")
)]
manu[, direction_call := ifelse(
  is.na(direction_match), "missing",
  ifelse(direction_match, "match", "mismatch")
)]
manu2 <- manu[, .(
  metric = metric_label,
  expected = expected_short,
  observed = observed_short,
  direction_call,
  beta_asd = signif(lm_beta_asd, 4),
  wilcox_p = signif(wilcox_p, 4),
  n_asd,
  n_control
)]
fwrite(manu2, file.path(OUTDIR, "04_direction_consistency_summary_manuscript.tsv"), sep = "\t")

# -----------------------------
# 6) summary plot
# -----------------------------
plot_dt <- pretty[!is.na(aligned_effect)]

if (nrow(plot_dt) > 0) {
  ord <- rev(c(
    "All microglia: candidate - reference",
    "All microglia: candidate-state score",
    "All microglia: reference-like score",
    "Score-defined cells: candidate proportion",
    "Score-defined cells: reference proportion",
    "Score-defined cells: log2(candidate/reference)",
    "Score-defined cells: candidate vs reference shift"
  ))
  plot_dt[, metric_label := factor(metric_label, levels = ord)]
  plot_dt[, sig_class := ifelse(!is.na(wilcox_p) & wilcox_p < 0.05, "Wilcox p < 0.05", "Wilcox p >= 0.05")]
  plot_dt[, match_class := ifelse(direction_match, "Direction matched", "Direction mismatched")]
  plot_dt[, label_text := ifelse(
    is.na(wilcox_p),
    observed_direction,
    paste0(observed_direction, " | p=", signif(wilcox_p, 3))
  )]

  xmax <- max(abs(plot_dt$aligned_effect), na.rm = TRUE)
  if (!is.finite(xmax) || xmax == 0) xmax <- 1

  p <- ggplot(plot_dt, aes(x = aligned_effect, y = metric_label)) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.5, colour = "grey40") +
    geom_segment(aes(x = 0, xend = aligned_effect, yend = metric_label), linewidth = 0.7, colour = "grey60") +
    geom_point(aes(shape = sig_class, fill = match_class), size = 3.8, colour = "black", stroke = 0.4) +
    geom_text(aes(label = label_text, hjust = ifelse(aligned_effect >= 0, -0.05, 1.05)), size = 3.1) +
    scale_shape_manual(values = c("Wilcox p < 0.05" = 21, "Wilcox p >= 0.05" = 24)) +
    scale_fill_manual(values = c("Direction matched" = "#4DAF4A", "Direction mismatched" = "#E41A1C")) +
    coord_cartesian(xlim = c(-1.2 * xmax, 1.5 * xmax), clip = "off") +
    labs(
      x = "Expected-direction aligned effect (positive = concordant)",
      y = NULL,
      title = "Package 4b v2: direction-consistency summary on the shared comparable sample set",
      subtitle = "All metrics restricted to the same score-comparable validation samples"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.title = element_blank(),
      plot.margin = margin(10, 95, 10, 10)
    )

  ggsave(file.path(OUTDIR, "05_direction_consistency_summary_plot.pdf"), p, width = 11.5, height = 5.2)
  ggsave(file.path(OUTDIR, "05_direction_consistency_summary_plot.png"), p, width = 11.5, height = 5.2, dpi = 320)
}

# -----------------------------
# 7) run notes
# -----------------------------
sink(file.path(OUTDIR, "06_run_notes.txt"))
cat("Package 4b v2 manual summary\n")
cat("============================\n")
cat("Input files:\n")
cat("- ", f_cell, "\n", sep = "")
cat("- ", f_prop, "\n", sep = "")
cat("- ", f_shift, "\n", sep = "")
cat("\nShared comparable samples were defined as the intersection of sample IDs across 12/15/16.\n")
cat("All metrics were forced onto the same sample set.\n")
cat("log2(candidate/reference) was computed directly from prop_candidate and prop_reference.\n\n")
cat("Detected columns:\n")
print(metric_detect_dt)
cat("\nShared sample counts:\n")
print(sample_count_dt)
sink()

message("Package 4b v2 manual done.")
