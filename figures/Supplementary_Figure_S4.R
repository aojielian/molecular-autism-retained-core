#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(stringr)
  library(scales)
  library(optparse)
})


project_dir_default <- Sys.getenv("MA_PROJECT_DIR", ".")
option_list <- list(
  make_option("--root", type = "character", default = Sys.getenv("MA_DRUG_OUTDIR", file.path(project_dir_default, "results", "Package9_PerturbationPrioritization")), help = "Package9 drug-prioritization directory [default: %default]"),
  make_option("--outdir", type = "character", default = file.path(project_dir_default, "results", "Supplementary_Figures", "S4_Perturbation"), help = "Output directory [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))
root <- opt$root
outdir <- opt$outdir
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) {
  fread(path, sep = "\t", header = TRUE, fill = TRUE, quote = "", data.table = TRUE)
}

short_label <- function(x, n = 45) {
  x <- gsub("\\s+", " ", x)
  x <- trimws(x)
  str_trunc(x, width = n)
}

primary <- read_tsv(file.path(root, "tables", "05_primary_reversal_combined.tsv"))
sens    <- read_tsv(file.path(root, "tables", "10_sensitivity_reversal_combined.tsv"))
buckets <- read_tsv(file.path(root, "postprocessed", "07_primary_mechanism_buckets.tsv"))

for (dt in list(primary, sens)) {
  if (!"perturbagen_norm" %in% names(dt)) {
    dt[, perturbagen_norm := tolower(gsub("\\s+", " ", Term_candDown))]
  }
}

primary[, label := short_label(perturbagen_norm, 40)]
sens[, label := short_label(perturbagen_norm, 40)]

# Top hits
primary_top <- primary[order(-reversal_score)][1:min(12, .N)]
sens_top    <- sens[order(-reversal_score)][1:min(12, .N)]

primary_top[, label := factor(label, levels = rev(label))]
sens_top[, label := factor(label, levels = rev(label))]

# Scatter labels: keep representative, manuscript-relevant candidates
label_hits <- unique(c(
  primary[order(-reversal_score)][1:min(12, .N), perturbagen_norm],
  c("vemurafenib", "plx4032 db05238", "lipopolysaccharide", "etanercept db00005",
    "metformin db00331", "curcumin", "pioglitazone db01132")
))
primary_scatter <- copy(primary)
primary_scatter[, annotate_hit := perturbagen_norm %in% label_hits]

col_primary <- "#C97B84"
col_sens <- "#B98E6C"
col_point <- "#BFBFBF"

pA <- ggplot(primary_top, aes(x = reversal_score, y = label)) +
  geom_col(fill = col_primary, width = 0.72) +
  geom_text(aes(label = sprintf("%.2f", reversal_score)),
            hjust = -0.1, size = 3.2) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(
    title = "A. Primary perturbation-based reversal candidates",
    x = "Reversal score",
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

# Calculate safe annotation position inside panel
x_max <- max(primary_scatter$neglog10_adjP_down, na.rm = TRUE)
y_max <- max(primary_scatter$neglog10_adjP_up, na.rm = TRUE)

pB <- ggplot(primary_scatter, aes(x = neglog10_adjP_down, y = neglog10_adjP_up)) +
  geom_point(color = col_point, alpha = 0.8, size = 1.8) +
  geom_point(
    data = primary_scatter[annotate_hit == TRUE],
    color = col_primary, size = 2.1
  ) +
  geom_text_repel(
    data = primary_scatter[annotate_hit == TRUE],
    aes(label = perturbagen_norm),
    size = 3.1,
    box.padding = 0.25,
    point.padding = 0.18,
    segment.size = 0.25,
    max.overlaps = Inf
  ) +
  labs(
    title = "B. Primary reversal scatter highlights bidirectional candidates",
    x = expression(-log[10]("candidate-program adj. P (drug-down)")),
    y = expression(-log[10]("reference-program adj. P (drug-up)"))
  ) +
  annotate(
    "label",
    x = x_max * 0.98,
    y = y_max * 1.02,
    hjust = 1,
    vjust = 1,
    label = paste0(
      "Top primary hit: ", primary[which.max(reversal_score), perturbagen_norm], "\n",
      "Top dual-axis signal: ", primary[order(-reversal_score)][1, perturbagen_norm], "\n",
      "MAPK/BRAF-related examples: PLX4032, vemurafenib"
    ),
    size = 3.0,
    label.size = 0.25
  ) +
  coord_cartesian(
    xlim = c(0, x_max * 1.05),
    ylim = c(0, y_max * 1.06),
    clip = "off"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

pC <- ggplot(sens_top, aes(x = reversal_score, y = label)) +
  geom_col(fill = col_sens, width = 0.72) +
  geom_text(aes(label = sprintf("%.2f", reversal_score)),
            hjust = -0.1, size = 3.2) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(
    title = "C. Sensitivity pooled-DE reversal candidates",
    x = "Reversal score",
    y = NULL
  ) +
  annotate(
    "label",
    x = Inf, y = -Inf, hjust = 1.02, vjust = -0.1,
    label = paste0(
      "Sensitivity top hit: ", sens[which.max(reversal_score), perturbagen_norm], "\n",
      "Interpret as hypothesis-generating prioritization"
    ),
    size = 3.0,
    label.size = 0.25
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

fig <- (pA / pB / pC) + plot_layout(heights = c(1.0, 1.15, 1.0))

ggsave(file.path(outdir, "Supplementary_Figure_S3.png"), fig, width = 12.5, height = 14.0, dpi = 300)
ggsave(file.path(outdir, "Supplementary_Figure_S3.pdf"), fig, width = 12.5, height = 14.0)

fwrite(primary_top[, .(perturbagen_norm, reversal_score, support_strength,
                       Adjusted_p_value_candDown, Adjusted_p_value_refUp)],
       file.path(outdir, "Supplementary_Figure_S3_primary_top_hits.tsv"), sep = "\t")
fwrite(sens_top[, .(perturbagen_norm, reversal_score, support_strength,
                    Adjusted_p_value_candDown, Adjusted_p_value_refUp)],
       file.path(outdir, "Supplementary_Figure_S3_sensitivity_top_hits.tsv"), sep = "\t")
if (!is.null(buckets) && nrow(buckets) > 0) {
  fwrite(buckets, file.path(outdir, "Supplementary_Figure_S3_primary_mechanism_buckets.tsv"), sep = "\t")
}

message("Wrote Supplementary Figure S3 to: ", outdir)
