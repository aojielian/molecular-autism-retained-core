#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(forcats)
  library(grid)
})

PROJECT_DIR <- "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision"
DIR_P4  <- file.path(PROJECT_DIR, "Package4_CrossCohortValidation_Velmeshev_v2")
DIR_P4B <- file.path(PROJECT_DIR, "Package4b_DirectionConsistencySummary_v2_manual")
DIR_P5  <- file.path(PROJECT_DIR, "Package5_DirectionalConcordance_GlobalSummary")
DIR_P7  <- file.path(PROJECT_DIR, "Package7_MinimalRobustnessAnalysis")
OUTDIR  <- file.path(PROJECT_DIR, "Figure5_Velmeshev_Validation_v5_4panel_final")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

read_dt <- function(path) {
  if (!file.exists(path)) stop("Missing input: ", path)
  if (grepl("\\.gz$", path)) {
    fread(cmd = paste("gzip -dc", shQuote(path)), sep = "\t", header = TRUE, data.table = TRUE)
  } else {
    fread(path, sep = "\t", header = TRUE, data.table = TRUE)
  }
}

f_map          <- file.path(DIR_P4,  "11_validation_subcluster_mapping_descriptive.tsv")
f_counts       <- file.path(DIR_P4,  "03_validation_microglia_dx_counts.tsv")
f_pretty       <- file.path(DIR_P4B, "03_direction_consistency_summary_pretty.tsv")
f_shared       <- file.path(DIR_P4B, "00b_shared_sample_counts.tsv")
f_global       <- file.path(DIR_P5,  "02_global_direction_tests.tsv")
f_stability    <- file.path(DIR_P7,  "02_score_defined_direction_stability.tsv")
f_loo          <- file.path(DIR_P7,  "03_leave_one_metric_out_global_summary.tsv")

map_dt    <- read_dt(f_map)
count_dt  <- read_dt(f_counts)
pretty_dt <- read_dt(f_pretty)
shared_dt <- read_dt(f_shared)
global_dt <- read_dt(f_global)
stab_dt   <- read_dt(f_stability)
loo_dt    <- read_dt(f_loo)

fmt_num <- function(x, digits = 3) {
  if (length(x) == 0) return("NA")
  x <- as.numeric(x)[1]
  if (!is.finite(x)) return("NA")
  formatC(x, digits = digits, format = "fg", flag = "#")
}

fmt_p <- function(x) {
  vapply(x, function(xx) {
    if (!is.finite(xx)) return("NA")
    if (xx < 1e-4) return(format(xx, scientific = TRUE, digits = 2))
    formatC(xx, digits = 3, format = "fg", flag = "#")
  }, character(1))
}

get_global <- function(name) {
  val <- global_dt[test_name == name, value]
  if (length(val) == 0) return(NA_real_)
  as.numeric(val[1])
}

short_metric <- function(x) {
  x <- as.character(x)
  x[x == "All microglia: candidate - reference"] <- "Candidate-reference"
  x[x == "All microglia: candidate-state score"] <- "Candidate-state score"
  x[x == "All microglia: reference-like score"] <- "Lower reference-like score"
  x[x == "Score-defined cells: candidate proportion"] <- "Candidate proportion"
  x[x == "Score-defined cells: reference proportion"] <- "Lower reference proportion"
  x[x == "Score-defined cells: log2(candidate/reference)"] <- "log2(candidate/reference)"
  x[x == "Score-defined cells: candidate vs reference shift"] <- "Candidate-reference shift"
  x
}

theme_clean <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold")
    )
}

# A. mapping heatmap
map_dt[, validation_cluster := as.character(validation_cluster)]
setorder(map_dt, corr_delta_c2_minus_c0)
map_dt[, validation_cluster_f := factor(validation_cluster, levels = rev(validation_cluster))]

heat_dt <- rbindlist(list(
  map_dt[, .(validation_cluster_f, cluster = "Discovery cluster 0", corr = corr_to_disc_cluster0)],
  map_dt[, .(validation_cluster_f, cluster = "Discovery cluster 2", corr = corr_to_disc_cluster2)]
))

pA <- ggplot(heat_dt, aes(x = cluster, y = validation_cluster_f, fill = corr)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.3f", corr)), size = 3.0) +
  scale_fill_gradient(low = "#F0F0F0", high = "#4C78A8",
                      limits = c(min(heat_dt$corr, na.rm = TRUE), max(heat_dt$corr, na.rm = TRUE))) +
  theme_clean(11) +
  theme(
    axis.title = element_blank(),
    legend.position = "right",
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold")
  ) +
  labs(title = "A. Validation subclusters correlate with both discovery clusters", fill = "Spearman r")

# B. delta plot
threshold_delta <- 0.05
max_delta <- max(map_dt$corr_delta_c2_minus_c0, na.rm = TRUE)
max_cluster <- map_dt[which.max(corr_delta_c2_minus_c0), validation_cluster][1]
count_dt[, dx := as.character(dx)]
n_asd_cells <- count_dt[dx == "ASD", N][1]
n_ctrl_cells <- count_dt[dx == "Control", N][1]
n_total_cells <- sum(count_dt$N, na.rm = TRUE)

ann_B <- paste0(
  "Validation microglia: ", comma(n_total_cells), " cells\n",
  "ASD ", comma(n_asd_cells), " | Control ", comma(n_ctrl_cells), "\n",
  "Max delta = ", fmt_num(max_delta, 5)
)

pB <- ggplot(map_dt, aes(x = corr_delta_c2_minus_c0, y = validation_cluster_f)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.5, colour = "grey45") +
  geom_vline(xintercept = threshold_delta, linetype = 3, linewidth = 0.6, colour = "grey55") +
  geom_segment(aes(x = 0, xend = corr_delta_c2_minus_c0, yend = validation_cluster_f), linewidth = 0.7, colour = "grey65") +
  geom_point(size = 2.7, shape = 21, fill = "white", colour = "black", stroke = 0.5) +
  annotate(
    "label",
    x = 0.026,
    y = 2.0,
    label = ann_B,
    hjust = 0.5,
    vjust = 1,
    size = 3.0,
    label.r = unit(0.12, "lines"),
    fill = "white"
  ) +
  theme_clean(11) +
  theme(axis.title.y = element_blank()) +
  labs(
    title = "B. Differential preference for discovery cluster 2 is weak",
    x = "Correlation delta (discovery cluster 2 minus cluster 0)",
    y = NULL
  )

# C. seven metrics aligned effect + aggregate stats
pretty_dt[, metric_label_short := short_metric(metric_label)]
pretty_dt[, nominal_p := fifelse(!is.na(lm_p), lm_p, wilcox_p)]
pretty_dt[, metric_group := ifelse(grepl("^All microglia", metric_label), "All microglia", "Score-defined")]
ord <- rev(c(
  "Candidate-reference",
  "Candidate-state score",
  "Lower reference-like score",
  "Candidate proportion",
  "Lower reference proportion",
  "log2(candidate/reference)",
  "Candidate-reference shift"
))
pretty_dt[, metric_label_short := factor(metric_label_short, levels = ord)]

xmax_d <- max(abs(pretty_dt$aligned_effect), na.rm = TRUE)
if (!is.finite(xmax_d) || xmax_d == 0) xmax_d <- 1

n_total_metrics <- get_global("n_total_metrics")
n_matched       <- get_global("n_direction_matched")
binom_one       <- get_global("binomial_test_one_sided_p")
stouffer_one    <- get_global("stouffer_one_sided_p")
shared_asd      <- get_global("n_asd_strict_shared")
shared_ctrl     <- get_global("n_control_strict_shared")

ann_C <- paste0(
  "Matched metrics: ", n_matched, "/", n_total_metrics,
  "; strict shared samples: ASD ", shared_asd, ", Control ", shared_ctrl,
  "; one-sided binomial p = ", fmt_p(binom_one),
  "; one-sided Stouffer p = ", fmt_p(stouffer_one)
)

pC <- ggplot(pretty_dt, aes(x = aligned_effect, y = metric_label_short)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.5, colour = "grey45") +
  geom_segment(aes(x = 0, xend = aligned_effect, yend = metric_label_short), linewidth = 0.75, colour = "grey65") +
  geom_point(aes(fill = metric_group), shape = 21, size = 3.7, colour = "black", stroke = 0.4) +
  scale_fill_manual(values = c("All microglia" = "#9ecae1", "Score-defined" = "#f4b6c2")) +
  coord_cartesian(xlim = c(-0.10 * xmax_d, 1.25 * xmax_d), clip = "off") +
  theme_clean(11) +
  theme(
    axis.title.y = element_blank(),
    legend.position = "top",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, colour = "black", margin = margin(b = 6))
  ) +
  labs(
    title = "C. Seven validation metrics show concordant score-based support",
    subtitle = ann_C,
    x = "Expected-direction aligned effect (positive = concordant)",
    y = NULL
  )

# D. leave-one-metric-out stability
# D2. leave-one-metric-out lollipop
loo_sub <- copy(loo_dt)[omitted_metric != "none"]
loo_sub[, omitted_short := short_metric(omitted_metric)]
loo_sub[, omitted_short := factor(omitted_short, levels = rev(short_metric(pretty_dt$metric_label)))]
loo_sub[, stouffer_score := -log10(pmax(stouffer_one_sided_p, 1e-300))]
base_stouffer_p <- loo_dt[omitted_metric == "none", stouffer_one_sided_p][1]
base_binom_p    <- loo_dt[omitted_metric == "none", binomial_one_sided_p][1]
base_score      <- -log10(pmax(base_stouffer_p, 1e-300))
loo_sub[, matched_fraction := n_direction_matched / n_total_metrics]

xmax_d2 <- max(c(loo_sub$stouffer_score, base_score), na.rm = TRUE)

pD2 <- ggplot(loo_sub, aes(x = stouffer_score, y = omitted_short)) +
  geom_vline(xintercept = base_score, linetype = 3, linewidth = 0.6, colour = "grey45") +
  geom_segment(aes(x = 0, xend = stouffer_score, yend = omitted_short), linewidth = 0.8, colour = "grey70") +
  geom_point(size = 3.1, shape = 21, fill = "white", colour = "black", stroke = 0.5) +
  theme_clean(10.5) +
  theme(
    axis.title.y = element_blank(),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_blank()
  ) +
  labs(
    title = "D. Leave-one-metric-out stability",
    subtitle = NULL,
    x = expression(-log[10]("one-sided Stouffer p")),
    y = NULL
  )

pD <- pD2

upper <- pA | pB
lower <- pC | pD
fig5 <- (upper / lower) + plot_layout(heights = c(0.92, 1.18), widths = c(1, 1))

ggsave(file.path(OUTDIR, "Figure5_combined.pdf"), fig5, width = 15.8, height = 10.8)
ggsave(file.path(OUTDIR, "Figure5_combined.png"), fig5, width = 15.8, height = 10.8, dpi = 300)

message("Done. Outputs written to: ", OUTDIR)
