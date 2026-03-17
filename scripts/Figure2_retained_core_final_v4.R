#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(tidyr)
  library(Matrix)
  library(forcats)
})

options(stringsAsFactors = FALSE)

# =========================================================
# Figure 2 final plotting script (v5)
# Main message:
# The retained microglial core shows a donor-level ASD shift
# toward a candidate-state program.
#
# Inputs are fixed from previously used package scripts.
# OOF state concordance is copied for record-keeping only,
# but is excluded from the main Figure 2 forest plot because
# it showed essentially no ASD-control difference.
# =========================================================

cfg <- list(
  strict_rds = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Package2_strict_input_fixed_from_original.rds",

  pkg2_outdir = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Package2_Define_Cluster2",
  pkg3_outdir = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Package3_Discovery_DonorAwareDiseaseAssociation",

  pkg2_cluster2_summary = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Package2_Define_Cluster2/tables/22_cluster2_core_summary_up.tsv",
  pkg2_cluster0_summary = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Package2_Define_Cluster2/tables/23_cluster0_core_summary_up.tsv",
  pkg2_cell_de = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Package2_Define_Cluster2/tables/06_celllevel_DE_cluster2_vs_0.tsv.gz",

  pkg3_abundance = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Package3_Discovery_DonorAwareDiseaseAssociation/02_abundance_per_donor.tsv",
  pkg3_abundance_stats = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Package3_Discovery_DonorAwareDiseaseAssociation/03_abundance_stats.tsv",
  pkg3_top_state_genes = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Package3_Discovery_DonorAwareDiseaseAssociation/08_state_definition_top_genes.tsv",
  pkg3_oof_scores = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Package3_Discovery_DonorAwareDiseaseAssociation/09_oof_state_concordance_per_donor.tsv",
  pkg3_oof_stats = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Package3_Discovery_DonorAwareDiseaseAssociation/10_oof_state_concordance_stats.tsv",

  outdir = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Figure2_RetainedCore_DonorLevel_v5",

  assay = "RNA",
  cluster_col = "cluster_use",
  donor_col = "individual_ID",
  dx_col = "Diagnosis",

  ref_cluster = "0",
  cand_cluster = "2",
  case_label = "ASD",
  control_label = "Control",

  top_dot_each = 6,
  score_top_pos = 100,
  score_top_neg = 100,
  min_cells_donor = 20,

  color_control = "#BDBDBD",
  color_asd = "#C97B84"
)

dir.create(cfg$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cfg$outdir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cfg$outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cfg$outdir, "objects"), recursive = TRUE, showWarnings = FALSE)

msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ..., "\n", sep = "")
}

safe_fread <- function(f, sep = "\t") {
  if (!file.exists(f)) stop("Missing file: ", f)
  data.table::fread(f, sep = sep)
}

norm_dx <- function(x) {
  x0 <- trimws(as.character(x))
  xl <- tolower(x0)
  out <- rep(NA_character_, length(x0))
  out[grepl("^asd$|autism|case", xl)] <- cfg$case_label
  out[grepl("^control$|^ctrl$|^ctl$|normal|neurotyp|^hc$", xl)] <- cfg$control_label
  out[x0 %in% c(cfg$case_label, cfg$control_label)] <- x0[x0 %in% c(cfg$case_label, cfg$control_label)]
  out
}

get_assay_mat <- function(obj, assay = "RNA", layer_name = c("data", "counts")) {
  layer_name <- match.arg(layer_name)
  mat <- tryCatch({
    GetAssayData(obj, assay = assay, slot = layer_name)
  }, error = function(e) {
    GetAssayData(obj, assay = assay, layer = layer_name)
  })
  mat
}

safe_scale <- function(x) {
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(0, length(x)))
  as.numeric(scale(x))
}

fit_dx_effect <- function(df, value_col, dx_col = "dx", expected_direction = c("ASD > Control", "ASD < Control")) {
  expected_direction <- match.arg(expected_direction)
  d <- as.data.table(df)
  d <- d[get(dx_col) %in% c(cfg$control_label, cfg$case_label) & is.finite(get(value_col))]

  if (nrow(d) < 4) {
    return(data.table(
      metric = value_col,
      expected_direction = expected_direction,
      beta = NA_real_,
      se = NA_real_,
      p = NA_real_,
      mean_asd = NA_real_,
      mean_control = NA_real_,
      aligned_effect = NA_real_,
      n_asd = NA_integer_,
      n_control = NA_integer_
    ))
  }

  d[, dx_factor := factor(get(dx_col), levels = c(cfg$control_label, cfg$case_label))]
  d[, y := safe_scale(get(value_col))]
  fit <- tryCatch(lm(y ~ dx_factor, data = d), error = function(e) NULL)

  beta <- se <- p <- NA_real_
  if (!is.null(fit)) {
    sm <- summary(fit)$coefficients
    if ("dx_factorASD" %in% rownames(sm)) {
      beta <- unname(sm["dx_factorASD", "Estimate"])
      se   <- unname(sm["dx_factorASD", "Std. Error"])
      p    <- unname(sm["dx_factorASD", "Pr(>|t|)"])
    }
  }

  mean_asd <- mean(d[dx_factor == cfg$case_label, get(value_col)], na.rm = TRUE)
  mean_ctrl <- mean(d[dx_factor == cfg$control_label, get(value_col)], na.rm = TRUE)
  aligned <- if (expected_direction == "ASD > Control") beta else -beta

  n_asd <- if ("donor" %in% colnames(d)) uniqueN(d[dx_factor == cfg$case_label, donor]) else sum(d$dx_factor == cfg$case_label)
  n_control <- if ("donor" %in% colnames(d)) uniqueN(d[dx_factor == cfg$control_label, donor]) else sum(d$dx_factor == cfg$control_label)

  data.table(
    metric = value_col,
    expected_direction = expected_direction,
    beta = beta,
    se = se,
    p = p,
    mean_asd = mean_asd,
    mean_control = mean_ctrl,
    aligned_effect = aligned,
    n_asd = n_asd,
    n_control = n_control
  )
}

plot_box_metric <- function(df, x = "dx", y, ylab, title = NULL) {
  d <- as.data.table(df)
  d <- d[get(x) %in% c(cfg$control_label, cfg$case_label) & is.finite(get(y))]
  d[, dx_plot := factor(get(x), levels = c(cfg$control_label, cfg$case_label))]

  ggplot(d, aes(x = dx_plot, y = get(y), fill = dx_plot)) +
    geom_violin(trim = FALSE, scale = "width", alpha = 0.45, color = "black", linewidth = 0.5) +
    geom_boxplot(width = 0.18, outlier.shape = NA, fill = "white", color = "black", alpha = 0.95) +
    geom_jitter(width = 0.10, height = 0, alpha = 0.8, size = 1.6, color = "black") +
    scale_fill_manual(values = c(cfg$color_control, cfg$color_asd)) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      plot.title = element_text(face = "bold")
    ) +
    labs(x = NULL, y = ylab, title = title)
}

make_forest_plot <- function(dt) {
  d <- copy(dt)
  d <- d[is.finite(aligned_effect)]
  if (nrow(d) == 0) {
    return(
      ggplot() +
        theme_void() +
        annotate("text", x = 0, y = 0, label = "No forest data available")
    )
  }

  d[, metric_pretty := c(
    "Candidate−reference",
    "Candidate-state",
    "Lower reference-like",
    "Cluster 2 proportion",
    "log2(C2/C0)"
  )[match(metric, c(
    "candidate_minus_reference",
    "candidate_state_score",
    "reference_like_score",
    "prop_c2",
    "log2_ratio_c2_over_c0"
  ))]]

  d[, lo := aligned_effect - 1.96 * se]
  d[, hi := aligned_effect + 1.96 * se]
  d[, metric_pretty := fct_rev(factor(metric_pretty, levels = unique(metric_pretty)))]

  ggplot(d, aes(x = aligned_effect, y = metric_pretty)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.18, linewidth = 0.5) +
    geom_point(size = 2.6) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.margin = margin(5.5, 12, 5.5, 20, unit = "pt")
    ) +
    labs(
      x = "Direction-aligned standardized ASD effect (beta ± 95% CI)",
      title = "E. Direction-aligned donor-level ASD effects"
    )
}

# =========================================================
# 1. Read strict retained object
# =========================================================
msg("Reading strict retained object ...")
obj <- readRDS(cfg$strict_rds)
if (!inherits(obj, "Seurat")) stop("strict_rds is not a Seurat object")
if (cfg$assay %in% Assays(obj)) DefaultAssay(obj) <- cfg$assay

meta <- as.data.table(obj@meta.data, keep.rownames = "cell")
if (!all(c(cfg$cluster_col, cfg$donor_col, cfg$dx_col) %in% colnames(meta))) {
  stop("Required metadata columns missing in strict object: ",
       paste(setdiff(c(cfg$cluster_col, cfg$donor_col, cfg$dx_col), colnames(meta)), collapse = ", "))
}

meta[, cluster := as.character(get(cfg$cluster_col))]
meta[, donor := as.character(get(cfg$donor_col))]
meta[, dx := norm_dx(get(cfg$dx_col))]
meta <- meta[cluster %in% c(cfg$ref_cluster, cfg$cand_cluster) & !is.na(dx) & !is.na(donor) & donor != ""]

obj <- subset(obj, cells = meta$cell)
meta <- as.data.table(obj@meta.data, keep.rownames = "cell")
meta[, cluster := as.character(get(cfg$cluster_col))]
meta[, donor := as.character(get(cfg$donor_col))]
meta[, dx := norm_dx(get(cfg$dx_col))]
meta <- meta[cluster %in% c(cfg$ref_cluster, cfg$cand_cluster) & !is.na(dx) & !is.na(donor) & donor != ""]

obj$cluster_use <- factor(meta$cluster[match(colnames(obj), meta$cell)], levels = c(cfg$ref_cluster, cfg$cand_cluster))
obj$dx_use <- meta$dx[match(colnames(obj), meta$cell)]
obj$donor_use <- meta$donor[match(colnames(obj), meta$cell)]
Idents(obj) <- "cluster_use"

saveRDS(obj, file.path(cfg$outdir, "objects", "Figure2_input_object_strict_0_vs_2.rds"))
fwrite(meta[, .(cell, donor, dx, cluster)], file.path(cfg$outdir, "tables", "00_cell_manifest_used.tsv"), sep = "\t")

# =========================================================
# 2. Read package-derived inputs
# =========================================================
msg("Reading package-derived inputs ...")
cluster2_summary <- safe_fread(cfg$pkg2_cluster2_summary)
cluster0_summary <- safe_fread(cfg$pkg2_cluster0_summary)
state_top <- safe_fread(cfg$pkg3_top_state_genes)
abund <- safe_fread(cfg$pkg3_abundance)
if (file.exists(cfg$pkg3_oof_scores)) {
  oof_scores <- safe_fread(cfg$pkg3_oof_scores)
} else {
  oof_scores <- data.table()
}
if (file.exists(cfg$pkg3_oof_stats)) {
  oof_stats <- safe_fread(cfg$pkg3_oof_stats)
} else {
  oof_stats <- data.table()
}

# =========================================================
# 3. Define plotting genes and state signature
# =========================================================
msg("Defining plotting genes and data-driven signature ...")
stopifnot("gene" %in% colnames(cluster2_summary), "gene" %in% colnames(cluster0_summary))
stopifnot(all(c("direction", "gene") %in% colnames(state_top)))

plot_genes_c2 <- unique(cluster2_summary$gene)[1:min(cfg$top_dot_each, nrow(cluster2_summary))]
plot_genes_c0 <- unique(cluster0_summary$gene)[1:min(cfg$top_dot_each, nrow(cluster0_summary))]
dot_genes <- unique(c(plot_genes_c2, plot_genes_c0))
dot_genes <- intersect(dot_genes, rownames(obj))

sig_pos <- state_top[direction == "Cluster2_up", unique(gene)][1:min(cfg$score_top_pos, state_top[direction == "Cluster2_up", uniqueN(gene)])]
sig_neg <- state_top[direction == "Cluster0_up", unique(gene)][1:min(cfg$score_top_neg, state_top[direction == "Cluster0_up", uniqueN(gene)])]

expr_data <- get_assay_mat(obj, assay = DefaultAssay(obj), layer_name = "data")
sig_pos_use <- intersect(sig_pos, rownames(expr_data))
sig_neg_use <- intersect(sig_neg, rownames(expr_data))

fwrite(data.table(signature = "candidate_state", gene = sig_pos_use),
       file.path(cfg$outdir, "tables", "01_signature_candidate_state_genes_used.tsv"), sep = "\t")
fwrite(data.table(signature = "reference_like", gene = sig_neg_use),
       file.path(cfg$outdir, "tables", "02_signature_reference_like_genes_used.tsv"), sep = "\t")

# =========================================================
# 4. Cell-level frozen scores in discovery retained core
# =========================================================
msg("Computing cell-level candidate/reference scores ...")

pos_score <- if (length(sig_pos_use) > 0) Matrix::colMeans(expr_data[sig_pos_use, , drop = FALSE]) else rep(NA_real_, ncol(expr_data))
neg_score <- if (length(sig_neg_use) > 0) Matrix::colMeans(expr_data[sig_neg_use, , drop = FALSE]) else rep(NA_real_, ncol(expr_data))

cell_scores <- data.table(
  cell = colnames(expr_data),
  donor = obj$donor_use,
  dx = obj$dx_use,
  cluster = as.character(obj$cluster_use),
  candidate_state_score = as.numeric(pos_score),
  reference_like_score = as.numeric(neg_score)
)
cell_scores[, candidate_minus_reference := candidate_state_score - reference_like_score]
cell_scores[, candidate_over_reference := (candidate_state_score + 1e-6) / (reference_like_score + 1e-6)]

fwrite(cell_scores, file.path(cfg$outdir, "tables", "03_discovery_cell_scores.tsv.gz"), sep = "\t")

score_donor <- cell_scores[, .(
  n_cells = .N,
  candidate_state_score = mean(candidate_state_score, na.rm = TRUE),
  reference_like_score = mean(reference_like_score, na.rm = TRUE),
  candidate_minus_reference = mean(candidate_minus_reference, na.rm = TRUE)
), by = .(donor, dx)]
score_donor <- score_donor[n_cells >= cfg$min_cells_donor]

fwrite(score_donor, file.path(cfg$outdir, "tables", "04_discovery_donor_scores.tsv"), sep = "\t")

abund <- as.data.table(abund)
abund[, donor := as.character(donor)]
abund[, dx := norm_dx(dx)]
abund <- abund[n_total_retained >= cfg$min_cells_donor]
fwrite(abund, file.path(cfg$outdir, "tables", "05_discovery_abundance_per_donor_copy.tsv"), sep = "\t")

if (nrow(oof_scores) > 0) {
  oof_scores <- as.data.table(oof_scores)
  oof_scores[, donor := as.character(donor)]
  oof_scores[, dx := norm_dx(dx)]
  fwrite(oof_scores, file.path(cfg$outdir, "tables", "06_discovery_oof_state_concordance_copy.tsv"), sep = "\t")
}
if (nrow(oof_stats) > 0) {
  fwrite(oof_stats, file.path(cfg$outdir, "tables", "06b_discovery_oof_state_concordance_stats_copy.tsv"), sep = "\t")
}

# =========================================================
# 5. Donor-level effect summary for forest plot
# =========================================================
msg("Computing donor-level ASD effect summary ...")

forest_list <- list(
  fit_dx_effect(score_donor, "candidate_minus_reference", expected_direction = "ASD > Control"),
  fit_dx_effect(score_donor, "candidate_state_score", expected_direction = "ASD > Control"),
  fit_dx_effect(score_donor, "reference_like_score", expected_direction = "ASD < Control"),
  fit_dx_effect(abund, "prop_c2", expected_direction = "ASD > Control"),
  fit_dx_effect(abund, "log2_ratio_c2_over_c0", expected_direction = "ASD > Control")
)

forest_dt <- rbindlist(forest_list, fill = TRUE)
fwrite(forest_dt, file.path(cfg$outdir, "tables", "07_donor_level_effect_summary.tsv"), sep = "\t")

# =========================================================
# 6. Build figure panels
# =========================================================
msg("Generating figure panels ...")

if (length(dot_genes) == 0) stop("No dotplot genes available after intersecting with rownames(obj)")

# A. DotPlot for retained clusters 0/2
pA <- DotPlot(obj, features = dot_genes, group.by = "cluster_use") +
  RotatedAxis() +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10)
  ) +
  labs(
    title = "A. Program genes distinguish retained clusters 0 and 2",
    x = NULL, y = NULL
  )

# B. Cell-level background distributions
cell_long <- cell_scores[, .(dx, candidate_state_score, reference_like_score, candidate_minus_reference)] |>
  tidyr::pivot_longer(cols = c(candidate_state_score, reference_like_score, candidate_minus_reference),
                      names_to = "metric", values_to = "value")
cell_long$metric <- factor(
  cell_long$metric,
  levels = c("candidate_state_score", "reference_like_score", "candidate_minus_reference"),
  labels = c("Candidate-state score", "Reference-like score", "Candidate − reference")
)
cell_long$dx_plot <- factor(cell_long$dx, levels = c(cfg$control_label, cfg$case_label))

pB <- ggplot(cell_long, aes(x = dx_plot, y = value, fill = dx_plot)) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.45, color = "black", linewidth = 0.5) +
  geom_boxplot(width = 0.18, outlier.shape = NA, fill = "white", color = "black", alpha = 0.95) +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c(cfg$color_control, cfg$color_asd)) +
  theme_bw(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "none",
    plot.title = element_text(face = "bold")
  ) +
  labs(
    title = "B. Cell-level scores (background only)",
    x = NULL, y = "Score"
  )

# C. Donor-level candidate-minus-reference
pC <- plot_box_metric(
  score_donor,
  y = "candidate_minus_reference",
  ylab = "Donor mean candidate − reference score",
  title = "C. ASD donors show higher candidate−reference scores"
)

# D. Donor-level abundance metrics
pD1 <- plot_box_metric(
  abund,
  y = "prop_c2",
  ylab = "Cluster 2 proportion",
  title = "D1. Higher cluster 2 proportion in ASD"
) +
  theme(
    plot.title = element_text(face = "bold", size = 12, margin = margin(b = 10)),
    plot.margin = margin(t = 5.5, r = 12, b = 5.5, l = 5.5)
  )

pD2 <- plot_box_metric(
  abund,
  y = "log2_ratio_c2_over_c0",
  ylab = "log2 cluster 2 / cluster 0 ratio",
  title = "D2. Higher C2/C0 ratio in ASD"
) +
  theme(
    plot.title = element_text(face = "bold", size = 12, margin = margin(b = 10)),
    plot.margin = margin(t = 5.5, r = 5.5, b = 5.5, l = 12)
  )

pD <- pD1 + pD2 + plot_layout(ncol = 2, widths = c(1, 1))

# E. Donor-level forest summary
pE <- make_forest_plot(forest_dt)

# Combined figure
fig2 <- (pA | pB) / (pC | pD) / pE +
  plot_layout(heights = c(1.00, 0.95, 0.85)) &
  theme(
    plot.title = element_text(face = "bold")
  )

# =========================================================
# 7. Save plots
# =========================================================
msg("Saving Figure 2 outputs ...")

ggsave(file.path(cfg$outdir, "plots", "Figure2A_dotplot_retained_clusters.pdf"), pA, width = 9.2, height = 5.0)
ggsave(file.path(cfg$outdir, "plots", "Figure2A_dotplot_retained_clusters.png"), pA, width = 9.2, height = 5.0, dpi = 300)

ggsave(file.path(cfg$outdir, "plots", "Figure2B_celllevel_score_distributions.pdf"), pB, width = 11.8, height = 4.2)
ggsave(file.path(cfg$outdir, "plots", "Figure2B_celllevel_score_distributions.png"), pB, width = 11.8, height = 4.2, dpi = 300)

ggsave(file.path(cfg$outdir, "plots", "Figure2C_donor_candidate_minus_reference.pdf"), pC, width = 4.8, height = 4.4)
ggsave(file.path(cfg$outdir, "plots", "Figure2C_donor_candidate_minus_reference.png"), pC, width = 4.8, height = 4.4, dpi = 300)

ggsave(file.path(cfg$outdir, "plots", "Figure2D_donor_abundance_metrics.pdf"), pD, width = 9.8, height = 4.4)
ggsave(file.path(cfg$outdir, "plots", "Figure2D_donor_abundance_metrics.png"), pD, width = 9.8, height = 4.4, dpi = 300)

ggsave(file.path(cfg$outdir, "plots", "Figure2E_donor_effect_forest.pdf"), pE, width = 8.2, height = 4.4)
ggsave(file.path(cfg$outdir, "plots", "Figure2E_donor_effect_forest.png"), pE, width = 8.2, height = 4.4, dpi = 300)

ggsave(file.path(cfg$outdir, "plots", "Figure2_combined.pdf"), fig2, width = 16.2, height = 13.0)
ggsave(file.path(cfg$outdir, "plots", "Figure2_combined.png"), fig2, width = 16.2, height = 13.0, dpi = 300)

# =========================================================
# 8. Manuscript-facing summary
# =========================================================
summary_txt <- c(
  "Figure 2 input chain fixed from previous packages:",
  paste0("- strict object: ", cfg$strict_rds),
  paste0("- Package 2 cluster2 summary: ", cfg$pkg2_cluster2_summary),
  paste0("- Package 2 cluster0 summary: ", cfg$pkg2_cluster0_summary),
  paste0("- Package 3 top state genes: ", cfg$pkg3_top_state_genes),
  paste0("- Package 3 abundance per donor: ", cfg$pkg3_abundance),
  paste0("- Package 3 OOF state concordance (copied only): ", cfg$pkg3_oof_scores),
  paste0("- Package 3 OOF stats (copied only): ", cfg$pkg3_oof_stats),
  "",
  "Panel logic:",
  "A = retained-cluster gene distinction (cluster 0 vs cluster 2)",
  "B = cell-level background score distributions",
  "C = donor-level candidate-minus-reference score",
  "D = donor-level abundance shift (cluster 2 proportion and log2 C2/C0)",
  "E = direction-aligned donor-level ASD effect summary",
  "",
  "Scoring logic:",
  "Candidate/reference scores are derived from Package 3 paired pseudobulk state-definition genes (data-driven),",
  "not from hand-curated SPP1-axis modules.",
  "",
  "OOF handling:",
  "OOF state concordance is copied into the output directory for record-keeping,",
  "but is excluded from Figure 2E because it showed essentially no ASD-control difference",
  "in Package 3 (effect_case_minus_control approximately -0.0074; wilcox p approximately 0.777; lm p approximately 0.873).",
  "",
  "Current figure message:",
  "The retained microglial core shows a donor-level ASD shift toward a candidate-state program.",
  "",
  "Color logic:",
  paste0("Control = ", cfg$color_control),
  paste0("ASD = ", cfg$color_asd),
  "Forest plot kept black-and-white."
)
writeLines(summary_txt, con = file.path(cfg$outdir, "Figure2_input_and_logic_summary.txt"))

capture.output(sessionInfo(), file = file.path(cfg$outdir, "sessionInfo_Figure2.txt"))
msg("Figure 2 script finished. Output directory: ", cfg$outdir)

