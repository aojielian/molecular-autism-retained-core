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
