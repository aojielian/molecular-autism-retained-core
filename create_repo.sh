#!/bin/bash
set -euo pipefail

BASE="/Users/aojie/PROJECTs/molecular_autism"
REPO="${BASE}/MA-risk-microglia-ASD"

echo "Creating repo at: ${REPO}"
mkdir -p "${REPO}"

# -----------------------------
# 1. Create directory structure
# -----------------------------
mkdir -p "${REPO}"/{data_dictionary,scripts/scripts_old,results/{figures,supplementary_tables,logs},docs,env}
mkdir -p "${REPO}"/scripts/utils

# -----------------------------
# 2. Create README
# -----------------------------
cat > "${REPO}/README.md" <<'EOF'
# MA-risk-microglia-ASD

This repository contains the analysis code, documentation, and reproducibility materials for the Molecular Autism revision project on ASD-associated SPP1+ lipid-metabolic microglial states.

## Repository purpose
This local repository is used for:
- organizing figure-generation scripts
- documenting data sources and paths
- tracking manuscript-related code revisions
- preparing materials for GitHub/Zenodo release

## Important note
Raw input data are **not stored in this repository**.
The analyses are run on the HPC server using the original paths under `/gpfs/...`.

## Main server-side data locations
- PsychENCODE discovery cohort:
  `/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/PsychENCODE_Science_2024/`
- Velmeshev validation cohort:
  `/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/Velmeshev_scRNA/`
- Bulk RNA-seq cohorts:
  `/gpfs/hpc/home/lijc/lianaoj/autism_bulk_RNA/`

## Planned script modules
- `01_fig1_discovery_clusters.R`
- `02_fig2_validation_donor_level.R`
- `03_fig3_cellchat_bulk.R`
- `04_fig4_metabolism.R`
- `05_fig5_tf_activity.R`
- `06_supp_tables.R`

## Reproducibility
This repository will eventually include:
- session information
- software/package versions
- curated gene-set definitions
- supplementary table generation code
EOF

# -----------------------------
# 3. Create LICENSE placeholder
# -----------------------------
cat > "${REPO}/LICENSE" <<'EOF'
MIT License

Copyright (c) 2026

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction...
EOF

# -----------------------------
# 4. Create .gitignore
# -----------------------------
cat > "${REPO}/.gitignore" <<'EOF'
# macOS
.DS_Store

# R
.Rhistory
.RData
.Rproj.user

# logs
*.log
*.out
*.err

# temporary
tmp/
temp/
temp_f2/

# large results
results/figures/*.pdf
results/figures/*.png
results/supplementary_tables/*.csv
results/supplementary_tables/*.tsv
results/supplementary_tables/*.xlsx
results/logs/*

# compressed / binary
*.rds
*.rda
*.mtx
*.mtx.gz
*.zip
*.gz
*.tar
*.tar.gz

# editor
.vscode/
.idea/
EOF

# -----------------------------
# 5. Create data dictionary
# -----------------------------
cat > "${REPO}/data_dictionary/data_dictionary.md" <<'EOF'
# Data Dictionary

## 1. PsychENCODE discovery cohort
Base path:
`/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/PsychENCODE_Science_2024/`

Files:
- `counts_matrix.mtx.gz`
- `counts_barcodes.tsv.gz`
- `counts_features.tsv.gz`
- `meta.tsv`
- `SeuratObj_RawConstructed.rds`
- `Step3_Microglia/Microglia_Clustered.rds`
- `Step3_Microglia/Dual_Signatures_List.rds`

## 2. Velmeshev validation cohort
Base path:
`/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/Velmeshev_scRNA/`

Files:
- `rawMatrix.zip`
- `meta.tsv`

## 3. Bulk RNA-seq cohorts
Base path:
`/gpfs/hpc/home/lijc/lianaoj/autism_bulk_RNA/`

Files:
- `GSE102741_RNAseq_GRCh38_expr_pheno.rds`
- `GSE64018_RNAseq_countlevel_expr_pheno.rds`

## 4. Reference resources
Base path:
`/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/`

Files:
- `dorothea_hs.rda`
- `safari.genes.Sand1.txt`
- `MEF2C_ASD_denovo.txt`
EOF

# -----------------------------
# 6. Create docs
# -----------------------------
cat > "${REPO}/docs/reproducibility_notes.md" <<'EOF'
# Reproducibility Notes

## Goals
- make analysis scripts version-controlled
- document software versions
- separate discovery, validation, metabolism, and TF modules
- support manuscript revision and point-by-point response

## To add
- exact package versions
- sessionInfo()
- HPC submission scripts
- supplementary table generation workflow
EOF

cat > "${REPO}/docs/software_versions.md" <<'EOF'
# Software and Package Versions

To be filled after collecting exact server environment details.

Suggested fields:
- R version
- Seurat
- ggplot2
- dplyr
- tidyr
- data.table
- Matrix
- CellChat
- clusterProfiler
- org.Hs.eg.db
- decoupleR
- dorothea
- tibble
EOF

cat > "${REPO}/docs/manuscript_change_log.md" <<'EOF'
# Manuscript Change Log

## Planned major changes
- revise title to lower causal strength
- unify homeostatic cluster definition
- add donor-level validation plots
- add all-cluster marker summaries
- moderate eicosanoid/plasma interpretation
- revise TF activity wording to match displayed rankings
EOF

# -----------------------------
# 7. Create env files
# -----------------------------
cat > "${REPO}/env/sessionInfo.txt" <<'EOF'
Session info will be exported from the HPC environment and pasted here.
EOF

cat > "${REPO}/env/renv.lock" <<'EOF'
{}
EOF

# -----------------------------
# 8. Create utility scripts
# -----------------------------
cat > "${REPO}/scripts/utils/plotting_utils.R" <<'EOF'
# plotting_utils.R
# Common plotting helpers for Molecular Autism revision project

theme_ma <- function() {
  ggplot2::theme_classic() +
    ggplot2::theme(
      text = ggplot2::element_text(size = 12),
      axis.title = ggplot2::element_text(size = 12),
      axis.text = ggplot2::element_text(size = 10)
    )
}
EOF

cat > "${REPO}/scripts/utils/scoring_utils.R" <<'EOF'
# scoring_utils.R
# Helper functions for donor-level aggregation and module-score handling

aggregate_donor_scores <- function(meta_df, donor_col, group_col, score_cols) {
  stopifnot(all(c(donor_col, group_col, score_cols) %in% colnames(meta_df)))
  dplyr::as_tibble(meta_df) %>%
    dplyr::mutate(Donor = .data[[donor_col]], Group = .data[[group_col]]) %>%
    dplyr::group_by(Donor, Group) %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(score_cols), ~mean(.x, na.rm = TRUE)),
      n_cells = dplyr::n(),
      .groups = "drop"
    )
}
EOF

# -----------------------------
# 9. Create main R scripts
# -----------------------------
cat > "${REPO}/scripts/00_session_info.R" <<'EOF'
# 00_session_info.R
# Export session info from HPC environment

sink("env/sessionInfo.txt")
sessionInfo()
sink()
EOF

cat > "${REPO}/scripts/01_fig1_discovery_clusters.R" <<'EOF'
# 01_fig1_discovery_clusters.R
# Discovery cohort: microglial atlas, all-cluster markers, proportions

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(scales)
})

sc_dir  <- "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/PsychENCODE_Science_2024/"
out_dir <- "results/figures"
tab_dir <- "results/supplementary_tables"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

mg_clus <- readRDS(file.path(sc_dir, "Step3_Microglia/Microglia_Clustered.rds"))

# Main Figure 1B
my_cols <- setNames(rep("grey85", length(levels(mg_clus$seurat_clusters))), levels(mg_clus$seurat_clusters))
my_cols["2"] <- "#E41A1C"
p1b <- DimPlot(
  mg_clus, reduction = "umap", group.by = "seurat_clusters",
  cols = my_cols, label = TRUE, label.size = 6
) + theme_void() + theme(legend.position = "none")
ggsave(file.path(out_dir, "Fig1B_Microglia_UMAP.pdf"), p1b, width = 6, height = 6)

# Supplementary: all clusters unique colors
all_lvls <- levels(mg_clus$seurat_clusters)
cluster_cols <- setNames(scales::hue_pal()(length(all_lvls)), all_lvls)
p_umap_all <- DimPlot(
  mg_clus, reduction = "umap", group.by = "seurat_clusters",
  cols = cluster_cols, label = TRUE, label.size = 4
) + theme_void()
ggsave(file.path(out_dir, "Supp_MG_UMAP_AllClusters.pdf"), p_umap_all, width = 7, height = 6)

# Main Figure 1C
props <- mg_clus@meta.data %>%
  group_by(individual_ID, Diagnosis, seurat_clusters) %>%
  summarise(n = n(), .groups = "drop") %>%
  complete(individual_ID, seurat_clusters, fill = list(n = 0)) %>%
  group_by(individual_ID) %>%
  mutate(
    Diagnosis = unique(Diagnosis[!is.na(Diagnosis)]),
    Freq = n / sum(n)
  ) %>%
  filter(seurat_clusters == "2" & !is.na(Diagnosis)) %>%
  mutate(Diagnosis = factor(Diagnosis, levels = c("CTL", "ASD")))

p1c <- ggplot(props, aes(x = Diagnosis, y = Freq, fill = Diagnosis)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2) +
  scale_fill_manual(values = c("CTL" = "grey80", "ASD" = "#E41A1C")) +
  theme_classic() +
  labs(y = "Fraction", x = "") +
  theme(legend.position = "none")
ggsave(file.path(out_dir, "Fig1C_Proportion.pdf"), p1c, width = 3.5, height = 4.5)

# Supplementary: all-cluster proportions
props_all <- mg_clus@meta.data %>%
  group_by(individual_ID, Diagnosis, seurat_clusters) %>%
  summarise(n = n(), .groups = "drop") %>%
  complete(individual_ID, seurat_clusters, fill = list(n = 0)) %>%
  group_by(individual_ID) %>%
  mutate(
    Diagnosis = unique(Diagnosis[!is.na(Diagnosis)]),
    Freq = n / sum(n)
  ) %>%
  filter(!is.na(Diagnosis))

write.csv(props_all, file.path(tab_dir, "Supp_AllCluster_Proportions.csv"), row.names = FALSE)

# Main Figure 1D
features <- c("P2RY12", "CX3CR1", "TMEM119", "CSF1R", "SPP1", "GAS6", "APOE", "CD74", "C3", "LPL", "CTSB")
p1d <- DotPlot(
  mg_clus,
  features = intersect(features, rownames(mg_clus)),
  cols = c("lightgrey", "firebrick"),
  dot.scale = 6
) + theme_bw() + labs(x = "", y = "Cluster") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"))
ggsave(file.path(out_dir, "Fig1D_DotPlot.pdf"), p1d, width = 7.5, height = 5.5)

# Main Figure 1E
markers2 <- FindMarkers(mg_clus, ident.1 = "2", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
eg <- bitr(rownames(markers2), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
go_res <- enrichGO(gene = eg$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pvalueCutoff = 0.05, readable = TRUE)
go_df <- as.data.frame(go_res) %>% arrange(p.adjust) %>% head(10) %>% mutate(Description = substr(Description, 1, 45))
p1e <- ggplot(go_df, aes(x = reorder(Description, -log10(p.adjust)), y = -log10(p.adjust))) +
  geom_bar(stat = "identity", fill = "#377EB8") +
  coord_flip() + theme_classic() + labs(y = "-Log10(FDR)", x = "")
ggsave(file.path(out_dir, "Fig1E_GO.pdf"), p1e, width = 6, height = 4)

# Supplementary: all markers
all_markers <- FindAllMarkers(mg_clus, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(all_markers, file.path(tab_dir, "Supp_AllCluster_Markers.csv"), row.names = FALSE)
EOF

cat > "${REPO}/scripts/02_fig2_validation_donor_level.R" <<'EOF'
# 02_fig2_validation_donor_level.R
# Validation cohort: module scores, donor-level aggregation

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(data.table)
  library(Matrix)
})

raw_dir  <- "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/Velmeshev_scRNA/"
disc_dir <- "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/PsychENCODE_Science_2024/Step3_Microglia/"
out_dir  <- "results/figures"
tab_dir  <- "results/supplementary_tables"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)

unzip(file.path(raw_dir, "rawMatrix.zip"), exdir = "temp_f2")
counts <- Matrix::ReadMtx(
  mtx = "temp_f2/matrix.mtx",
  cells = "temp_f2/barcodes.tsv",
  features = "temp_f2/genes.tsv",
  feature.column = 2
)
meta <- fread(file.path(raw_dir, "meta.tsv"), data.table = FALSE)
rownames(meta) <- meta[[1]]
common <- intersect(colnames(counts), rownames(meta))
sc_val <- CreateSeuratObject(counts = counts[, common], meta.data = meta[common, ])
unlink("temp_f2", recursive = TRUE)

clust_col <- grep("cluster", colnames(sc_val@meta.data), value = TRUE, ignore.case = TRUE)[1]
diag_col  <- grep("diagnosis", colnames(sc_val@meta.data), value = TRUE, ignore.case = TRUE)[1]

Idents(sc_val) <- clust_col
mg_val <- subset(sc_val, idents = grep("Microglia|MG", unique(Idents(sc_val)), value = TRUE, ignore.case = TRUE))
mg_val <- NormalizeData(mg_val, verbose = FALSE)

raw_vec <- as.character(mg_val[[diag_col]])
mg_val$Group <- NA
mg_val$Group[raw_vec == "ASD"] <- "ASD"
mg_val$Group[raw_vec == "Control"] <- "Control"
mg_val <- subset(mg_val, subset = Group %in% c("ASD", "Control"))
mg_val$Group <- factor(mg_val$Group, levels = c("Control", "ASD"))

sigs <- readRDS(file.path(disc_dir, "Dual_Signatures_List.rds"))
mg_val <- AddModuleScore(mg_val, features = list(sigs$Cluster2_Up), name = "Risk_Score")

# NOTE:
# This still uses the original object key.
# After manuscript revision, update this name if the homeostatic reference changes from Cluster4 to Clusters0/1.
mg_val <- AddModuleScore(mg_val, features = list(sigs$Cluster4_Down), name = "Homeo_Score")

df <- mg_val@meta.data

# Per-cell plots
p2b <- ggplot(df, aes(x = Group, y = Risk_Score1, fill = Group)) +
  geom_violin(alpha = 0.5, color = NA) +
  geom_boxplot(width = 0.2) +
  scale_fill_manual(values = c("Control" = "grey70", "ASD" = "#E41A1C")) +
  theme_classic() + theme(legend.position = "none") +
  labs(y = "Risk MG Score", x = "")
ggsave(file.path(out_dir, "Fig2B_Risk_PerCell.pdf"), p2b, width = 3.5, height = 4.5)

p2c <- ggplot(df, aes(x = Group, y = Homeo_Score1, fill = Group)) +
  geom_violin(alpha = 0.5, color = NA) +
  geom_boxplot(width = 0.2) +
  scale_fill_manual(values = c("Control" = "grey70", "ASD" = "#377EB8")) +
  theme_classic() + theme(legend.position = "none") +
  labs(y = "Homeostatic Score", x = "")
ggsave(file.path(out_dir, "Fig2C_Homeo_PerCell.pdf"), p2c, width = 3.5, height = 4.5)

# Donor-level aggregation
donor_col <- grep("individual|donor|sample", colnames(mg_val@meta.data), value = TRUE, ignore.case = TRUE)[1]
if (length(donor_col) == 0 || is.na(donor_col)) donor_col <- colnames(mg_val@meta.data)[1]

df_donor <- mg_val@meta.data %>%
  mutate(Donor = .data[[donor_col]]) %>%
  group_by(Donor, Group) %>%
  summarise(
    Risk_Score = mean(Risk_Score1, na.rm = TRUE),
    Homeo_Score = mean(Homeo_Score1, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  )

write.csv(df_donor, file.path(tab_dir, "Validation_DonorLevel_Scores.csv"), row.names = FALSE)

p2b_donor <- ggplot(df_donor, aes(x = Group, y = Risk_Score, fill = Group)) +
  geom_violin(alpha = 0.4, color = NA) +
  geom_boxplot(width = 0.18, outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 1.8) +
  scale_fill_manual(values = c("Control" = "grey70", "ASD" = "#E41A1C")) +
  theme_classic() + theme(legend.position = "none") +
  labs(y = "Donor-level Risk MG Score", x = "")
ggsave(file.path(out_dir, "Fig2B_Risk_DonorLevel.pdf"), p2b_donor, width = 3.8, height = 4.5)

p2c_donor <- ggplot(df_donor, aes(x = Group, y = Homeo_Score, fill = Group)) +
  geom_violin(alpha = 0.4, color = NA) +
  geom_boxplot(width = 0.18, outlier.shape = NA) +
  geom_jitter(width = 0.12, size = 1.8) +
  scale_fill_manual(values = c("Control" = "grey70", "ASD" = "#377EB8")) +
  theme_classic() + theme(legend.position = "none") +
  labs(y = "Donor-level Homeostatic Score", x = "")
ggsave(file.path(out_dir, "Fig2C_Homeo_DonorLevel.pdf"), p2c_donor, width = 3.8, height = 4.5)
EOF

cat > "${REPO}/scripts/03_fig3_cellchat_bulk.R" <<'EOF'
# 03_fig3_cellchat_bulk.R
# Placeholder for CellChat and bulk correlation analyses
message("Fill in finalized Figure 3 workflow here.")
EOF

cat > "${REPO}/scripts/04_fig4_metabolism.R" <<'EOF'
# 04_fig4_metabolism.R
# Placeholder for immunometabolic and eicosanoid analyses
message("Fill in finalized Figure 4 workflow here.")
EOF

cat > "${REPO}/scripts/05_fig5_tf_activity.R" <<'EOF'
# 05_fig5_tf_activity.R
# Placeholder for TF activity inference and plotting
message("Fill in finalized Figure 5 workflow here.")
EOF

cat > "${REPO}/scripts/06_supp_tables.R" <<'EOF'
# 06_supp_tables.R
# Placeholder for supplementary table generation
message("Generate cohort summary, software versions, curated gene sets, and full rankings here.")
EOF

# -----------------------------
# 10. Create helper shell script
# -----------------------------
cat > "${REPO}/init_git.sh" <<'EOF'
#!/bin/bash
set -euo pipefail
git init
git add .
git commit -m "Initialize Molecular Autism revision repository scaffold"
EOF
chmod +x "${REPO}/init_git.sh"

echo "Done."
echo "Repository created at: ${REPO}"
