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
