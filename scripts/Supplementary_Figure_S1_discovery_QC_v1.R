#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

root <- "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision"
pkg1 <- file.path(root, "Package1_Figure1_DiscoveryAudit")
pkg1b <- file.path(root, "Package1b_ClusterAnnotation_ContaminantAudit")
outdir <- file.path(root, "Supplementary_Figures_v1", "S1_Discovery_QC")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

asd_col <- "#C97B84"
ctrl_col <- "#BDBDBD"

read_tsv <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path)
  if (grepl("\\.gz$", path)) fread(cmd = paste("gzip -dc", shQuote(path)), sep = "\t", header = TRUE)
  else fread(path, sep = "\t", header = TRUE)
}

find_col <- function(nms, pats) {
  for (p in pats) {
    hit <- grep(p, nms, perl = TRUE, ignore.case = TRUE, value = TRUE)
    if (length(hit) > 0) return(hit[1])
  }
  stop("Cannot find matching column. Candidates: ", paste(pats, collapse = "; "),
       "\nAvailable: ", paste(nms, collapse = ", "))
}

# -------- Panel A: discovery microglia UMAP by donor --------
obj <- readRDS(file.path(pkg1, "19_microglia_object_with_audit.rds"))
meta_info <- read_tsv(file.path(pkg1, "Package1_Audit_Metadata.tsv"))
# donor column was recorded in Package1_Audit_Metadata.tsv as item/value
donor_col <- character()
if ("item" %in% names(meta_info) && "value" %in% names(meta_info)) {
  donor_col <- meta_info$value[meta_info$item == "donor_col"]
}
if (length(donor_col) == 0 || is.na(donor_col[1]) || !donor_col[1] %in% colnames(obj@meta.data)) {
  donor_col <- find_col(
    colnames(obj@meta.data),
    c("^donor$", "donor_id", "individual", "subject", "patient", "sample", "orig.ident", "library", "case_id")
  )
}
donor_col <- as.character(donor_col[1])

pA <- DimPlot(
  obj,
  reduction = "umap",
  group.by = donor_col,
  raster = FALSE
) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold")
  ) +
  labs(title = "A. Discovery microglia UMAP by donor")

# -------- Panel B: donor-level proportions across all clusters --------
prop_df <- read_tsv(file.path(pkg1, "06_allcluster_donor_proportions.tsv.gz"))
dx_col <- find_col(names(prop_df), c("^dx$", "^Diagnosis$", "^diagnosis$"))
cluster_col <- find_col(names(prop_df), c("^cluster$"))
prop_col <- find_col(names(prop_df), c("^prop$"))
prop_df[[dx_col]] <- factor(prop_df[[dx_col]], levels = c("Control", "ASD"))

# order clusters numerically if possible
clv <- unique(as.character(prop_df[[cluster_col]]))
num <- suppressWarnings(as.numeric(clv))
if (all(!is.na(num))) clv <- as.character(sort(num)) else clv <- sort(clv)
prop_df[[cluster_col]] <- factor(prop_df[[cluster_col]], levels = clv)

pB <- ggplot(prop_df, aes(x = .data[[dx_col]], y = .data[[prop_col]], fill = .data[[dx_col]])) +
  geom_boxplot(outlier.shape = NA, width = 0.62, alpha = 0.9) +
  geom_jitter(width = 0.12, size = 0.8, alpha = 0.75) +
  facet_wrap(as.formula(paste("~", cluster_col)), scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c(Control = ctrl_col, ASD = asd_col)) +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold", size = 8),
    axis.text.x = element_text(angle = 20, hjust = 1),
    plot.title = element_text(face = "bold")
  ) +
  labs(
    title = "B. Donor-level proportions across all discovery clusters",
    x = NULL,
    y = "Proportion within donor microglia"
  )

# -------- Panel C: audit superclass counts --------
anno_df <- read_tsv(file.path(pkg1b, "tables", "22_cluster_annotation_summary_for_manuscript.tsv"))
super_col <- find_col(names(anno_df), c("^suggested_superclass$"))
count_df <- as.data.table(anno_df)[, .N, by = super_col]
setnames(count_df, super_col, "suggested_superclass")
count_df[, suggested_superclass := factor(
  suggested_superclass,
  levels = c("microglia","mixed_microglia","uncertain","contaminant")
)]
pC <- ggplot(count_df, aes(x = suggested_superclass, y = N, fill = suggested_superclass)) +
  geom_col(width = 0.72) +
  scale_fill_manual(values = c(
    microglia = "#1b9e77",
    mixed_microglia = "#7570b3",
    uncertain = "#999999",
    contaminant = "#d95f02"
  ), drop = FALSE) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 20, hjust = 1),
    plot.title = element_text(face = "bold")
  ) +
  labs(title = "C. Number of clusters by audit superclass", x = NULL, y = "Number of clusters")

# -------- Panel D: lineage score heatmap --------
lineage_long <- read_tsv(file.path(pkg1b, "tables", "03_lineage_scores_long.tsv"))
cluster_col2 <- find_col(names(lineage_long), c("^cluster$"))
lineage_col <- find_col(names(lineage_long), c("^lineage$"))
score_col <- find_col(names(lineage_long), c("^total_score$"))
clv2 <- unique(as.character(lineage_long[[cluster_col2]]))
num2 <- suppressWarnings(as.numeric(clv2))
if (all(!is.na(num2))) clv2 <- as.character(sort(num2)) else clv2 <- sort(clv2)
lineage_long[[cluster_col2]] <- factor(lineage_long[[cluster_col2]], levels = clv2)

pD <- ggplot(lineage_long, aes(x = .data[[lineage_col]], y = .data[[cluster_col2]], fill = .data[[score_col]])) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1),
    plot.title = element_text(face = "bold")
  ) +
  labs(title = "D. Lineage total scores by cluster", x = NULL, y = "Cluster", fill = "Total score")

fig <- (pA | pC) / (pB | pD) +
  plot_layout(widths = c(1.15, 0.85), heights = c(0.9, 1.1))

ggsave(file.path(outdir, "Supplementary_Figure_S1.pdf"), fig, width = 15.5, height = 12)
ggsave(file.path(outdir, "Supplementary_Figure_S1.png"), fig, width = 15.5, height = 12, dpi = 300)

message("Supplementary Figure S1 written to: ", outdir)

