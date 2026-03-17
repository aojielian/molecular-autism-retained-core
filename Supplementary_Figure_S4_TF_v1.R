#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(stringr)
})

root <- "/Users/aojie/PROJECTs/molecular_autism/drug/Package10_TF_prioritization_local"
outdir <- file.path(root, "supplementary_figure_S4")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) {
  fread(path, sep = "\t", header = TRUE, fill = TRUE, quote = "", data.table = TRUE)
}

short_tf <- function(x, n = 20) str_trunc(x, width = n)

cand <- read_tsv(file.path(root, "tables", "06_primary_candidate_tf_summary.tsv"))
pool <- read_tsv(file.path(root, "tables", "10_pooled_up_tf_summary.tsv"))
conv <- read_tsv(file.path(root, "tables", "13_candidate_plus_pooledUp_convergent_TFs.tsv"))

cand[, TF := short_tf(TF)]
pool[, TF := short_tf(TF)]
conv[, TF := short_tf(TF)]

cand_top <- cand[order(-aggregate_score)][1:min(12, .N)]
pool_top <- pool[order(-aggregate_score)][1:min(12, .N)]
conv_top <- conv[order(-joint_score)][1:min(15, .N)]

cand_top[, TF := factor(TF, levels = rev(TF))]
pool_top[, TF := factor(TF, levels = rev(TF))]
conv_top[, TF := factor(TF, levels = rev(TF))]

col_cand <- "#C97B84"
col_pool <- "#B98E6C"
col_conv <- "#8FA8C9"

pA <- ggplot(cand_top, aes(x = aggregate_score, y = TF)) +
  geom_col(fill = col_cand, width = 0.72) +
  geom_text(aes(label = n_libraries_supported), hjust = -0.15, size = 3.1) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.10))) +
  labs(
    title = "A. Top TFs for the primary candidate program",
    subtitle = "Numbers at bar ends denote supporting library count",
    x = "Aggregate TF score",
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

pB <- ggplot(pool_top, aes(x = aggregate_score, y = TF)) +
  geom_col(fill = col_pool, width = 0.72) +
  geom_text(aes(label = n_libraries_supported), hjust = -0.15, size = 3.1) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.10))) +
  labs(
    title = "B. Top TFs for pooled ASD-up genes",
    subtitle = "Numbers at bar ends denote supporting library count",
    x = "Aggregate TF score",
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

highlight_tfs <- c("STAT3", "SPI1", "PPARG", "CEBPB", "IRF8", "NR3C1")
conv_top[, highlight := TF %in% highlight_tfs]

pC <- ggplot(conv_top, aes(x = joint_score, y = TF)) +
  geom_col(aes(fill = highlight), width = 0.72, show.legend = FALSE) +
  scale_fill_manual(values = c(`TRUE` = "#D87474", `FALSE` = col_conv)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.10))) +
  labs(
    title = "C. Convergent TFs shared by the candidate program and pooled ASD-up genes",
    subtitle = "Highlighted regulators are emphasized in the manuscript narrative",
    x = "Joint convergence score",
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

fig <- (pA / pB / pC) + plot_layout(heights = c(1.0, 1.0, 1.15))

ggsave(file.path(outdir, "Supplementary_Figure_S4.png"), fig, width = 12.5, height = 14.0, dpi = 300)
ggsave(file.path(outdir, "Supplementary_Figure_S4.pdf"), fig, width = 12.5, height = 14.0)

fwrite(cand_top, file.path(outdir, "Supplementary_Figure_S4_candidate_topTFs.tsv"), sep = "\t")
fwrite(pool_top, file.path(outdir, "Supplementary_Figure_S4_pooledUp_topTFs.tsv"), sep = "\t")
fwrite(conv_top, file.path(outdir, "Supplementary_Figure_S4_convergent_TFs.tsv"), sep = "\t")

message("Wrote Supplementary Figure S4 to: ", outdir)
