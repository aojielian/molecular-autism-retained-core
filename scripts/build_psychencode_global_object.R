#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(Seurat)
})

options(stringsAsFactors = FALSE)

cfg <- list(
  indir = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/PsychENCODE_Science_2024",
  out_rds = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/PsychENCODE_Science_2024/PsychENCODE_global_object.rds",
  out_meta_preview = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/PsychENCODE_Science_2024/PsychENCODE_global_object_metadata_preview.tsv"
)

msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ..., "\n", sep = "")
}

# -----------------------------
# 1. read matrix / features / barcodes
# -----------------------------
msg("Reading features ...")
features <- fread(file.path(cfg$indir, "counts_features.tsv.gz"), header = FALSE)
barcodes <- fread(file.path(cfg$indir, "counts_barcodes.tsv.gz"), header = FALSE)

if (ncol(features) < 1) stop("counts_features.tsv.gz has no columns.")
if (ncol(barcodes) < 1) stop("counts_barcodes.tsv.gz has no columns.")

gene_ids <- as.character(features[[1]])
if (ncol(features) >= 2) {
  gene_names <- as.character(features[[2]])
} else {
  gene_names <- gene_ids
}
cell_ids <- as.character(barcodes[[1]])

# gene symbols may contain duplicates
gene_names[is.na(gene_names) | gene_names == ""] <- gene_ids[is.na(gene_names) | gene_names == ""]
gene_names <- make.unique(gene_names)

msg("Reading sparse matrix ...")
mat <- readMM(file.path(cfg$indir, "counts_matrix.mtx.gz"))

if (nrow(mat) != length(gene_names)) {
  stop("Row number of matrix does not match features.")
}
if (ncol(mat) != length(cell_ids)) {
  stop("Column number of matrix does not match barcodes.")
}

rownames(mat) <- gene_names
colnames(mat) <- cell_ids

msg("Matrix dimensions: ", nrow(mat), " genes x ", ncol(mat), " cells")

# -----------------------------
# 2. read metadata
# -----------------------------
msg("Reading metadata ...")
meta <- fread(file.path(cfg$indir, "meta.tsv"))

if (!"Cell_ID" %in% colnames(meta)) {
  stop("meta.tsv does not contain Cell_ID.")
}

meta$Cell_ID <- as.character(meta$Cell_ID)

# align metadata to matrix columns
common_cells <- intersect(colnames(mat), meta$Cell_ID)
msg("Common cells between matrix and meta: ", length(common_cells))

if (length(common_cells) == 0) {
  stop("No overlapping Cell_ID between matrix columns and meta.tsv")
}

mat <- mat[, common_cells, drop = FALSE]
meta <- meta[match(common_cells, meta$Cell_ID), ]

if (!all(meta$Cell_ID == colnames(mat))) {
  stop("Metadata alignment failed.")
}

# set rownames for Seurat metadata
meta_df <- as.data.frame(meta)
rownames(meta_df) <- meta_df$Cell_ID

# -----------------------------
# 3. create Seurat object
# -----------------------------
msg("Creating Seurat object ...")
obj <- CreateSeuratObject(
  counts = mat,
  meta.data = meta_df,
  project = "PsychENCODE_Science_2024"
)

# if meta already contains nCount_RNA / nFeature_RNA, keep them
# Seurat may regenerate these; that's okay

# -----------------------------
# 4. read and attach UMAP if possible
# -----------------------------
umap_file <- file.path(cfg$indir, "UMAP.coords.tsv.gz")
if (file.exists(umap_file)) {
  msg("Reading UMAP coordinates ...")
  umap_dt <- fread(umap_file, header = FALSE)

  added_umap <- FALSE

  # Case A: V1 = Cell_ID, V2/V3 = coordinates
  if (ncol(umap_dt) >= 3) {
    umap_dt[[1]] <- as.character(umap_dt[[1]])
    overlap_n <- sum(umap_dt[[1]] %in% colnames(obj))
    msg("UMAP V1 overlap with object cells: ", overlap_n)

    if (overlap_n > 0) {
      umap_sub <- umap_dt[umap_dt[[1]] %in% colnames(obj), 1:3, with = FALSE]
      colnames(umap_sub) <- c("Cell_ID", "UMAP_1", "UMAP_2")
      umap_sub <- umap_sub[match(colnames(obj), umap_sub$Cell_ID), ]

      keep <- !is.na(umap_sub$UMAP_1) & !is.na(umap_sub$UMAP_2)
      if (sum(keep) > 0) {
        emb <- as.matrix(umap_sub[keep, c("UMAP_1", "UMAP_2")])
        rownames(emb) <- umap_sub$Cell_ID[keep]
        obj[["umap"]] <- CreateDimReducObject(
          embeddings = emb,
          key = "UMAP_",
          assay = DefaultAssay(obj)
        )
        added_umap <- TRUE
        msg("UMAP added using Cell_ID matching.")
      }
    }
  }

  # Case B: same number/order as object cells, but no Cell_ID
  if (!added_umap && nrow(umap_dt) == ncol(obj) && ncol(umap_dt) >= 2) {
    emb <- as.matrix(umap_dt[, 1:2, with = FALSE])
    rownames(emb) <- colnames(obj)
    colnames(emb) <- c("UMAP_1", "UMAP_2")
    obj[["umap"]] <- CreateDimReducObject(
      embeddings = emb,
      key = "UMAP_",
      assay = DefaultAssay(obj)
    )
    added_umap <- TRUE
    msg("UMAP added by row-order matching.")
  }

  if (!added_umap) {
    msg("UMAP could not be confidently attached; object created without UMAP.")
  }
}

# -----------------------------
# 5. save
# -----------------------------
preview_cols <- intersect(
  c("Cell_ID", "individual_ID", "Diagnosis", "annotation", "Subtype", "Brain_Region", "Sex_Chromosome"),
  colnames(obj@meta.data)
)

fwrite(
  as.data.table(obj@meta.data[, preview_cols, drop = FALSE], keep.rownames = "rowname")[1:min(20, ncol(obj)), ],
  cfg$out_meta_preview,
  sep = "\t"
)

saveRDS(obj, cfg$out_rds)
msg("Saved global object to: ", cfg$out_rds)
msg("Done.")
