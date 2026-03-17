#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(scales)
  library(rlang)
})

# =========================================================
# Figure 3
# Discovery-stage disease association
# =========================================================

cfg <- list(
  indir = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Package3_Discovery_DonorAwareDiseaseAssociation",
  outdir = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Figure3_Discovery_DiseaseAssociation_v3",
  control_color = "#BDBDBD",
  asd_color = "#C97B84"
)

dir.create(cfg$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cfg$outdir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cfg$outdir, "tables"), recursive = TRUE, showWarnings = FALSE)

safe_read <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  if (grepl("\\.gz$", path)) {
    fread(cmd = paste("gzip -dc", shQuote(path)), sep = "\t", header = TRUE, data.table = TRUE)
  } else {
    fread(path, sep = "\t", header = TRUE, data.table = TRUE)
  }
}

std_dx <- function(x) {
  x <- as.character(x)
  out <- ifelse(x %in% c("ASD", "Case", "case", "Autism"), "ASD", "Control")
  factor(out, levels = c("Control", "ASD"))
}

p_fmt <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 1e-4) return(formatC(p, format = "e", digits = 2))
  sprintf("%.3f", p)
}

theme_fig <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      legend.title = element_text(face = "bold"),
      strip.background = element_rect(fill = "white", colour = "black")
    )
}

find_dx_col <- function(dt) {
  cand <- c("dx", "Diagnosis", "diagnosis", "group", "Group")
  hit <- cand[cand %in% colnames(dt)]
  if (length(hit) == 0) stop("Cannot find diagnosis column. Available columns: ", paste(colnames(dt), collapse = ", "))
  hit[1]
}

standardize_dx_col <- function(dt) {
  dx_col <- find_dx_col(dt)
  dt[, dx := std_dx(get(dx_col))]
  dt
}

# -----------------------------
# Input files
# -----------------------------
abund_file <- file.path(cfg$indir, "02_abundance_per_donor.tsv")
abund_stats_file <- file.path(cfg$indir, "03_abundance_stats.tsv")
delta_file <- file.path(cfg$indir, "18_pseudobulk_donorDelta_cluster2_minus_cluster0_ASD_vs_Control.tsv.gz")
donor_pairs_file <- file.path(cfg$indir, "18_pseudobulk_donorDelta_cluster2_minus_cluster0_ASD_vs_Control.donor_pairs.tsv")

abund <- safe_read(abund_file)
abund_stats <- safe_read(abund_stats_file)
delta <- safe_read(delta_file)
donor_pairs <- safe_read(donor_pairs_file)

message("[read] abund columns: ", paste(colnames(abund), collapse = ", "))
message("[read] donor_pairs columns: ", paste(colnames(donor_pairs), collapse = ", "))

abund <- standardize_dx_col(abund)
donor_pairs <- standardize_dx_col(donor_pairs)

# -----------------------------
# Stats extraction
# -----------------------------
get_stat_row <- function(metric_name) {
  x <- abund_stats[metric == metric_name]
  if (nrow(x) == 0) stop("Metric not found in 03_abundance_stats.tsv: ", metric_name)
  x[1]
}

stat_prop <- get_stat_row("prop_c2_within_retained_core")
stat_ratio <- get_stat_row("log2_ratio_c2_over_c0")

# -----------------------------
# Panel A / B
# -----------------------------
plot_box_metric <- function(dt, y, ylab, title, stat_row) {
  ann <- paste0(
    "ASD-Control beta = ", sprintf("%.3f", stat_row$effect_case_minus_control),
    "\nlm P = ", p_fmt(stat_row$lm_p)
  )

  y_rng <- range(dt[[y]], na.rm = TRUE)
  y_pos <- y_rng[2] + 0.04 * diff(y_rng)

  ggplot(dt, aes(x = dx, y = .data[[y]], fill = dx)) +
    geom_violin(width = 0.9, alpha = 0.55, colour = NA, trim = FALSE) +
    geom_boxplot(width = 0.20, outlier.shape = NA, fill = "white", colour = "black") +
    geom_jitter(width = 0.08, height = 0, size = 2.1, alpha = 0.90, colour = "black") +
    scale_fill_manual(values = c("Control" = cfg$control_color, "ASD" = cfg$asd_color)) +
    theme_fig(12) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 12, margin = margin(b = 8))
    ) +
    labs(title = title, x = NULL, y = ylab) +
    annotate(
      "label",
      x = 1.5,
      y = y_pos,
      label = ann,
      size = 3.25,
      hjust = 0.5,
      vjust = 1.1,
      fill = "white"
    )
}

pA <- plot_box_metric(
  abund,
  y = "prop_c2",
  ylab = "Cluster 2 proportion",
  title = "A. Donor-level cluster 2 proportion",
  stat_row = stat_prop
)

pB <- plot_box_metric(
  abund,
  y = "log2_ratio_c2_over_c0",
  ylab = "log2(C2/C0)",
  title = "B. Donor-level log2(C2/C0) ratio",
  stat_row = stat_ratio
)

# -----------------------------
# Panel C volcano
# -----------------------------
setorder(delta, PValue)
delta[, neglog10P := -log10(pmax(PValue, 1e-300))]
delta[, sig_cat := "Background"]
delta[FDR < 0.10, sig_cat := "FDR<0.10"]

n_fdr01 <- delta[FDR < 0.10, .N]
min_fdr <- min(delta$FDR, na.rm = TRUE)
n_pairs <- nrow(donor_pairs)

n_label <- 10L
if (n_fdr01 == 0L) {
  delta[1:n_label, sig_cat := "Top nominal"]
  lab_dt <- delta[1:n_label]
} else {
  lab_dt <- delta[FDR < 0.10][order(PValue)][1:min(n_label, .N)]
}

volcano_cols <- c(
  "Background" = "#CFCFCF",
  "Top nominal" = cfg$asd_color,
  "FDR<0.10" = "#7A1F2B"
)

ann_c <- paste0(
  n_pairs, " paired donors\n",
  n_fdr01, " genes with FDR < 0.10\n",
  "Minimum FDR = ", sprintf("%.3f", min_fdr)
)

x_lim <- quantile(abs(delta$logFC), probs = 0.995, na.rm = TRUE)
x_lim <- max(1.5, x_lim)
y_lim <- quantile(delta$neglog10P, probs = 0.995, na.rm = TRUE)
y_lim <- max(4.5, y_lim)

pC <- ggplot(delta, aes(x = logFC, y = neglog10P)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = 3, colour = "grey65") +
  geom_point(aes(color = sig_cat), size = 1.0, alpha = 0.75) +
  geom_text_repel(
    data = lab_dt,
    aes(label = gene),
    size = 3.0,
    max.overlaps = Inf,
    box.padding = 0.35,
    point.padding = 0.20,
    segment.size = 0.25,
    min.segment.length = 0
  ) +
  scale_color_manual(values = volcano_cols, breaks = c("FDR<0.10", "Top nominal", "Background")) +
  coord_cartesian(xlim = c(-x_lim, x_lim), ylim = c(0, y_lim)) +
  theme_fig(12) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 12, margin = margin(b = 8))
  ) +
  labs(
    title = "C. Donor-paired delta model shows no robust ASD-enhanced internal shift",
    x = expression(log[2]*" fold-change"~"([Cluster 2 - Cluster 0]~ASD~vs~Control)"),
    y = expression(-log[10](italic(P)))
  ) +
  annotate(
    "label",
    x = -0.96 * x_lim,
    y = 0.96 * y_lim,
    label = ann_c,
    hjust = 0,
    vjust = 1,
    size = 3.25,
    fill = "white"
  )

# -----------------------------
# Panel D summary panel
# -----------------------------
summary_dt <- data.table(
  endpoint = c(
    "Cluster 2 proportion",
    "log2(C2/C0) ratio",
    "Donor-paired delta model"
  ),
  finding = c(
    paste0(
      ifelse(stat_prop$effect_case_minus_control > 0, "Higher in ASD", "Lower in ASD"),
      "\nASD-Control beta = ", sprintf("%.3f", stat_prop$effect_case_minus_control),
      "; lm P = ", p_fmt(stat_prop$lm_p)
    ),
    paste0(
      ifelse(stat_ratio$effect_case_minus_control > 0, "Higher in ASD", "Lower in ASD"),
      "\nASD-Control beta = ", sprintf("%.3f", stat_ratio$effect_case_minus_control),
      "; lm P = ", p_fmt(stat_ratio$lm_p)
    ),
    paste0(
      n_pairs, " paired donors",
      "\n0 genes with FDR < 0.10",
      "\nNo robust dx-enhanced internal shift"
    )
  ),
  row = c(3, 2, 1),
  fill = c(alpha(cfg$asd_color, 0.10), alpha(cfg$asd_color, 0.10), "white")
)

pD <- ggplot(summary_dt) +
  geom_tile(
    aes(x = 1, y = row, fill = fill),
    width = 1.92,
    height = 0.84,
    color = "grey35",
    linewidth = 0.35,
    show.legend = FALSE
  ) +
  scale_fill_identity() +
  geom_text(
    aes(x = 0.34, y = row, label = endpoint),
    hjust = 0,
    vjust = 0.5,
    size = 3.85,
    fontface = "bold"
  ) +
  geom_text(
    aes(x = 1.08, y = row, label = finding),
    hjust = 0,
    vjust = 0.5,
    size = 3.35,
    lineheight = 1.05
  ) +
  xlim(0, 2.04) +
  ylim(0.5, 3.5) +
  theme_void(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 12, margin = margin(b = 8)),
    plot.margin = margin(6, 8, 6, 8)
  ) +
  labs(title = "D. Discovery-stage summary")

# -----------------------------
# Save supporting tables
# -----------------------------
fwrite(summary_dt, file.path(cfg$outdir, "tables", "Figure3_summary_panel.tsv"), sep = "\t")
fwrite(abund_stats, file.path(cfg$outdir, "tables", "Figure3_abundance_stats_copy.tsv"), sep = "\t")
fwrite(donor_pairs, file.path(cfg$outdir, "tables", "Figure3_donor_pairs_copy.tsv"), sep = "\t")

# -----------------------------
# Save individual plots
# -----------------------------
ggsave(file.path(cfg$outdir, "plots", "Figure3A_prop_c2.pdf"), pA, width = 5.2, height = 4.6)
ggsave(file.path(cfg$outdir, "plots", "Figure3A_prop_c2.png"), pA, width = 5.2, height = 4.6, dpi = 300)

ggsave(file.path(cfg$outdir, "plots", "Figure3B_log2ratio_c2_over_c0.pdf"), pB, width = 5.2, height = 4.6)
ggsave(file.path(cfg$outdir, "plots", "Figure3B_log2ratio_c2_over_c0.png"), pB, width = 5.2, height = 4.6, dpi = 300)

ggsave(file.path(cfg$outdir, "plots", "Figure3C_donorDelta_volcano.pdf"), pC, width = 6.4, height = 5.2)
ggsave(file.path(cfg$outdir, "plots", "Figure3C_donorDelta_volcano.png"), pC, width = 6.4, height = 5.2, dpi = 300)

ggsave(file.path(cfg$outdir, "plots", "Figure3D_summary_panel.pdf"), pD, width = 6.4, height = 5.2)
ggsave(file.path(cfg$outdir, "plots", "Figure3D_summary_panel.png"), pD, width = 6.4, height = 5.2, dpi = 300)

# -----------------------------
# Combine figure
# -----------------------------
fig3 <- (pA | pB) / (pC | pD) +
  plot_layout(widths = c(1, 1), heights = c(0.92, 1.08))

ggsave(file.path(cfg$outdir, "plots", "Figure3_combined.pdf"), fig3, width = 14.6, height = 11.2)
ggsave(file.path(cfg$outdir, "plots", "Figure3_combined.png"), fig3, width = 14.6, height = 11.2, dpi = 300)

# -----------------------------
# Analysis notes
# -----------------------------
notes <- c(
  "Figure 3 = discovery-stage disease association.",
  paste0("Paired donor count for delta model: ", n_pairs),
  paste0("Genes with FDR<0.10 in donor-delta model: ", n_fdr01),
  paste0("Minimum FDR in donor-delta model: ", sprintf("%.6f", min_fdr)),
  "Interpretation target: upward abundance tendency, but no robust diagnosis-enhanced internal state shift."
)
writeLines(notes, con = file.path(cfg$outdir, "Figure3_notes.txt"))

message("Figure 3 finished: ", cfg$outdir)
