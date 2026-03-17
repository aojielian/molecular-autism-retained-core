#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
})

infile <- "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/PsychENCODE_Science_2024/Step3_Microglia/Microglia_Clustered.rds"
outfile <- "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Package2_strict_input_fixed_from_original.rds"

obj <- readRDS(infile)

cat("Metadata columns:\n")
print(colnames(obj@meta.data))

# 优先使用已有 cluster_use；如果没有，就退回 seurat_clusters
if ("cluster_use" %in% colnames(obj@meta.data) &&
    sum(!is.na(obj@meta.data$cluster_use)) > 0) {
  cluster_col <- "cluster_use"
} else if ("seurat_clusters" %in% colnames(obj@meta.data) &&
           sum(!is.na(obj@meta.data$seurat_clusters)) > 0) {
  cluster_col <- "seurat_clusters"
} else {
  stop("Could not find a usable cluster column in original microglia object.")
}

obj$cluster_use <- as.character(obj@meta.data[[cluster_col]])
obj$cluster_use <- trimws(obj$cluster_use)
obj$cluster_use <- sub("^Cluster", "", obj$cluster_use, ignore.case = TRUE)
obj$cluster_use <- sub("^cl", "", obj$cluster_use, ignore.case = TRUE)
obj$cluster_use <- sub("\\.0$", "", obj$cluster_use)

cat("\nUsing cluster column:", cluster_col, "\n")
print(table(obj$cluster_use, useNA = "ifany"))

keep_clusters <- c("0", "2")
cells_keep <- colnames(obj)[obj$cluster_use %in% keep_clusters]

cat("\nRetained cells:", length(cells_keep), "\n")
if (length(cells_keep) == 0) stop("No cells found for clusters 0 and 2 in original object.")

obj_fix <- subset(obj, cells = cells_keep)
obj_fix$cluster_use <- factor(as.character(obj_fix$cluster_use), levels = keep_clusters)
Idents(obj_fix) <- obj_fix$cluster_use

cat("\nFixed object cluster table:\n")
print(table(obj_fix$cluster_use, useNA = "ifany"))

# 简单检查关键列
check_cols <- intersect(c("individual_ID", "Diagnosis", "annotation", "Subtype", "cluster_use", "seurat_clusters"), colnames(obj_fix@meta.data))
cat("\nPreview of key metadata columns:\n")
print(head(obj_fix@meta.data[, check_cols, drop = FALSE]))

saveRDS(obj_fix, outfile)
cat("\nSaved fixed strict object to:\n", outfile, "\n")
