#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(Matrix)
  library(patchwork)
})

options(stringsAsFactors = FALSE)

cfg <- list(
  microglia_rds  = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/PsychENCODE_Science_2024/Step3_Microglia/Microglia_Clustered.rds",
  outdir         = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Package1b_ClusterAnnotation_ContaminantAudit",
  assay          = "RNA",
  cluster_col    = NULL,
  donor_col      = NULL,
  dx_col         = NULL,
  marker_file    = NA_character_,   # 可填 Package1 的 marker 表；若为空则自动 FindAllMarkers
  top_n_markers  = 50,
  min_pct        = 0.10,
  logfc          = 0.25,
  delta_strict   = 0.35,            # microglia total score 比 top non-mg 高多少才判 strict microglia
  delta_mixed    = 0.15,            # microglia 与 top non-mg 接近时判 mixed
  min_mg_anchor_hits_strict = 2,
  min_nonmg_hits_contam     = 2
)

dir.create(cfg$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cfg$outdir, "plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cfg$outdir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cfg$outdir, "objects"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cfg$outdir, "marker_panels"), recursive = TRUE, showWarnings = FALSE)

msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ..., "\n", sep = "")
}

save_session_info <- function(path) {
  writeLines(capture.output(sessionInfo()), con = path)
}

find_first_col <- function(df, patterns) {
  nm <- names(df)
  for (p in patterns) {
    hit <- nm[grepl(p, nm, ignore.case = TRUE)]
    if (length(hit) > 0) return(hit[1])
  }
  NULL
}

normalize_dx <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  y <- x
  y[grepl("asd|autism|case", x, ignore.case = TRUE)] <- "ASD"
  y[grepl("control|ctrl|ctl|healthy|non[- ]?psy", x, ignore.case = TRUE)] <- "Control"
  y
}

choose_cluster_col <- function(obj, user_col = NULL) {
  meta <- obj@meta.data
  if (!is.null(user_col) && user_col %in% names(meta)) return(user_col)
  if ("seurat_clusters" %in% names(meta)) return("seurat_clusters")
  cand <- find_first_col(meta, c("^cluster$", "cluster", "seurat", "res\\.", "RNA_snn", "wsnn"))
  if (!is.null(cand)) return(cand)
  obj$cluster_auto_from_idents <- as.character(Idents(obj))
  return("cluster_auto_from_idents")
}

choose_donor_col <- function(meta, user_col = NULL) {
  if (!is.null(user_col) && user_col %in% names(meta)) return(user_col)
  find_first_col(meta, c("^donor$", "donor_id", "individual", "subject", "patient", "sample", "orig.ident", "library", "case_id"))
}

choose_dx_col <- function(meta, user_col = NULL) {
  if (!is.null(user_col) && user_col %in% names(meta)) return(user_col)
  find_first_col(meta, c("^dx$", "^group$", "diagnosis", "condition", "status", "phenotype", "disease"))
}

get_fc_col <- function(markers) {
  fc_candidates <- c("avg_log2FC", "avg_logFC", "log2FC", "logFC")
  hit <- fc_candidates[fc_candidates %in% colnames(markers)]
  if (length(hit) == 0) return(NULL)
  hit[1]
}

safe_scale <- function(x) {
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  if (sd(x, na.rm = TRUE) == 0) return(rep(0, length(x)))
  as.numeric(scale(x))
}

row_zscore_by_gene <- function(mat) {
  if (is.null(dim(mat))) return(mat)
  z <- t(apply(mat, 1, function(v) {
    if (all(is.na(v))) return(rep(NA_real_, length(v)))
    if (sd(v, na.rm = TRUE) == 0) return(rep(0, length(v)))
    as.numeric(scale(v))
  }))
  rownames(z) <- rownames(mat)
  colnames(z) <- colnames(mat)
  z
}

calc_wilcox_p <- function(df, value_col, group_col) {
  df <- as.data.frame(df)
  g <- unique(df[[group_col]])
  g <- g[!is.na(g)]
  if (length(g) != 2) return(NA_real_)
  x1 <- df[[value_col]][df[[group_col]] == g[1]]
  x2 <- df[[value_col]][df[[group_col]] == g[2]]
  x1 <- x1[!is.na(x1)]
  x2 <- x2[!is.na(x2)]
  if (length(x1) < 1 || length(x2) < 1) return(NA_real_)
  suppressWarnings(wilcox.test(x1, x2, exact = FALSE)$p.value)
}

mean_diff_asd_minus_ctrl <- function(df, value_col, group_col) {
  df <- as.data.frame(df)
  grp <- as.character(df[[group_col]])
  val <- df[[value_col]]
  if (!all(c("ASD", "Control") %in% unique(grp))) return(NA_real_)
  mean(val[grp == "ASD"], na.rm = TRUE) - mean(val[grp == "Control"], na.rm = TRUE)
}

# ------------------------------------------------------------
# Transparent marker panels for contaminant auditing
# anchor = stronger lineage-defining markers
# support = supportive but less exclusive markers
# ------------------------------------------------------------
panel_def <- list(
  microglia = list(
    anchor  = c("P2RY12","CX3CR1","TMEM119","CSF1R","AIF1","C1QA","C1QB","C1QC","TREM2","TYROBP"),
    support = c("MERTK","FCER1G","LAPTM5","SPP1","GAS6","APOE","CD74","C3","LGALS3","GPNMB","LPL","FABP5","CTSB","CTSD","LIPA","ABCA1","MSR1")
  ),
  t_cell = list(
    anchor  = c("CD3D","CD3E","CD247","TRBC1","TRBC2","THEMIS"),
    support = c("IL7R","LTB","ETS1")
  ),
  b_cell = list(
    anchor  = c("MS4A1","CD79A","CD79B","IGHM","PAX5"),
    support = c("BLK","FCRL1","CD22","BANK1")
  ),
  endothelial = list(
    anchor  = c("CLDN5","FLT1","KDR","VWF","EPAS1"),
    support = c("ESAM","RAMP2","PLVAP","EMCN")
  ),
  neuronal = list(
    anchor  = c("RBFOX3","SNAP25","SYT1","STMN2"),
    support = c("SLC17A7","GAD1","GAD2","DLG4","GRIN2B","CAMK2A")
  ),
  astrocyte = list(
    anchor  = c("GFAP","AQP4","ALDH1L1"),
    support = c("SLC1A2","SLC1A3","GJA1","ALDOC")
  ),
  oligodendrocyte = list(
    anchor  = c("MBP","MOG","PLP1"),
    support = c("MOBP","MAG","CNP")
  ),
  opc = list(
    anchor  = c("PDGFRA","VCAN","PTPRZ1"),
    support = c("CSPG4","OLIG1","OLIG2")
  ),
  pericyte = list(
    anchor  = c("PDGFRB","RGS5"),
    support = c("MCAM","COL1A1","COL1A2","TAGLN")
  ),
  erythroid = list(
    anchor  = c("HBB","HBA1","HBA2"),
    support = c("ALAS2","AHSP")
  )
)

lineage_names <- names(panel_def)
nonmg_lineages <- setdiff(lineage_names, "microglia")

lineage_marker_panel <- unique(c(
  panel_def$microglia$anchor, panel_def$microglia$support,
  panel_def$t_cell$anchor, panel_def$t_cell$support,
  panel_def$b_cell$anchor, panel_def$b_cell$support,
  panel_def$endothelial$anchor, panel_def$endothelial$support,
  panel_def$neuronal$anchor, panel_def$neuronal$support,
  panel_def$astrocyte$anchor, panel_def$astrocyte$support,
  panel_def$oligodendrocyte$anchor, panel_def$oligodendrocyte$support,
  panel_def$opc$anchor, panel_def$opc$support,
  panel_def$pericyte$anchor, panel_def$pericyte$support,
  panel_def$erythroid$anchor, panel_def$erythroid$support
))

# save panels for provenance
panel_export <- rbindlist(lapply(names(panel_def), function(lin) {
  data.frame(
    lineage = lin,
    tier = c(rep("anchor", length(panel_def[[lin]]$anchor)),
             rep("support", length(panel_def[[lin]]$support))),
    gene = c(panel_def[[lin]]$anchor, panel_def[[lin]]$support),
    stringsAsFactors = FALSE
  )
}))
fwrite(panel_export, file.path(cfg$outdir, "marker_panels", "transparent_lineage_marker_panels.tsv"), sep = "\t")

# -----------------------------
# scoring helpers
# -----------------------------
get_panel_scores <- function(mat, genes_anchor, genes_support) {
  anchor_use <- intersect(genes_anchor, rownames(mat))
  support_use <- intersect(genes_support, rownames(mat))

  if (length(anchor_use) > 0) {
    sub_anchor <- mat[anchor_use, , drop = FALSE]
    z_anchor <- row_zscore_by_gene(sub_anchor)
    anchor_score <- colMeans(z_anchor, na.rm = TRUE)
  } else {
    anchor_score <- setNames(rep(NA_real_, ncol(mat)), colnames(mat))
  }

  if (length(support_use) > 0) {
    sub_support <- mat[support_use, , drop = FALSE]
    z_support <- row_zscore_by_gene(sub_support)
    support_score <- colMeans(z_support, na.rm = TRUE)
  } else {
    support_score <- setNames(rep(NA_real_, ncol(mat)), colnames(mat))
  }

  # anchor weighted more strongly than support
  total_score <- 0.75 * anchor_score + 0.25 * support_score

  list(
    anchor_score = anchor_score,
    support_score = support_score,
    total_score = total_score,
    anchor_genes_used = anchor_use,
    support_genes_used = support_use
  )
}

classify_cluster_row <- function(row_df, cfg) {
  mg_total <- as.numeric(row_df$microglia_total_score)
  mg_anchor <- as.numeric(row_df$microglia_anchor_score)
  mg_support <- as.numeric(row_df$microglia_support_score)
  mg_anchor_hits <- as.numeric(row_df$n_microglia_anchor_hits_top50)
  mg_support_hits <- as.numeric(row_df$n_microglia_support_hits_top50)
  mg_hits_total <- mg_anchor_hits + mg_support_hits

  nonmg_scores <- setNames(
    sapply(nonmg_lineages, function(lin) as.numeric(row_df[[paste0(lin, "_total_score")]])),
    nonmg_lineages
  )
  nonmg_scores[is.na(nonmg_scores)] <- -Inf

  nonmg_hits <- setNames(
    sapply(nonmg_lineages, function(lin) {
      h1 <- as.numeric(row_df[[paste0("n_", lin, "_anchor_hits_top50")]])
      h2 <- as.numeric(row_df[[paste0("n_", lin, "_support_hits_top50")]])
      h1[is.na(h1)] <- 0
      h2[is.na(h2)] <- 0
      h1 + h2
    }),
    nonmg_lineages
  )

  top_nonmg_lineage <- names(which.max(nonmg_scores))[1]
  top_nonmg_score <- max(nonmg_scores, na.rm = TRUE)
  top_nonmg_hits <- max(nonmg_hits, na.rm = TRUE)

  activation_delta <- mg_support - mg_anchor
  mg_state <- ifelse(
    is.na(activation_delta), "microglia_unspecified",
    ifelse(activation_delta > 0.35, "microglia_activation_like",
           ifelse(activation_delta < -0.35, "microglia_homeostatic_like",
                  "microglia_intermediate"))
  )

  out <- list(
    suggested_superclass = "uncertain",
    suggested_lineage = "uncertain",
    suggested_label = "uncertain",
    keep_strict = FALSE,
    keep_lenient = FALSE,
    rationale = ""
  )

  # strict microglia
  if (!is.na(mg_total) &&
      mg_total >= top_nonmg_score + cfg$delta_strict &&
      mg_anchor_hits >= cfg$min_mg_anchor_hits_strict) {
    out$suggested_superclass <- "microglia"
    out$suggested_lineage <- "microglia"
    out$suggested_label <- mg_state
    out$keep_strict <- TRUE
    out$keep_lenient <- TRUE
    out$rationale <- paste0(
      "microglia_total_score exceeds top_nonmg_score by >= ", cfg$delta_strict,
      "; microglia anchor hits top50 = ", mg_anchor_hits
    )
    return(out)
  }

  # clear contaminant
  if (!is.na(top_nonmg_score) &&
      top_nonmg_score >= mg_total + cfg$delta_strict &&
      top_nonmg_hits >= cfg$min_nonmg_hits_contam) {
    out$suggested_superclass <- "contaminant"
    out$suggested_lineage <- top_nonmg_lineage
    out$suggested_label <- paste0("contaminant_", top_nonmg_lineage)
    out$keep_strict <- FALSE
    out$keep_lenient <- FALSE
    out$rationale <- paste0(
      "top_nonmg_lineage=", top_nonmg_lineage,
      " dominates; nonmg hits top50 = ", top_nonmg_hits
    )
    return(out)
  }

  # mixed microglia-like
  if (!is.na(mg_total) &&
      mg_total >= top_nonmg_score - cfg$delta_mixed &&
      mg_hits_total >= 1) {
    out$suggested_superclass <- "mixed_microglia"
    out$suggested_lineage <- paste0("microglia_plus_", top_nonmg_lineage)
    out$suggested_label <- paste0("mixed_", mg_state, "_with_", top_nonmg_lineage)
    out$keep_strict <- FALSE
    out$keep_lenient <- TRUE
    out$rationale <- paste0(
      "microglia_total_score close to top_nonmg_score within ", cfg$delta_mixed,
      "; mg hits top50 = ", mg_hits_total
    )
    return(out)
  }

  # uncertain but likely non-microglia leaning
  if (!is.na(top_nonmg_score) && top_nonmg_hits >= 1) {
    out$suggested_superclass <- "uncertain"
    out$suggested_lineage <- top_nonmg_lineage
    out$suggested_label <- paste0("uncertain_", top_nonmg_lineage, "_leaning")
    out$keep_strict <- FALSE
    out$keep_lenient <- FALSE
    out$rationale <- paste0("uncertain; top_nonmg_lineage=", top_nonmg_lineage)
    return(out)
  }

  out$rationale <- "no lineage strongly supported"
  out
}

# ------------------------------------------------------------
# 1. Read object
# ------------------------------------------------------------
msg("Reading Seurat object ...")
obj <- readRDS(cfg$microglia_rds)
if (!inherits(obj, "Seurat")) stop("Input RDS is not a Seurat object.")

if (cfg$assay %in% Assays(obj)) DefaultAssay(obj) <- cfg$assay

cluster_col <- choose_cluster_col(obj, cfg$cluster_col)
obj$cluster_use <- as.character(obj@meta.data[[cluster_col]])
Idents(obj) <- obj$cluster_use

meta <- obj@meta.data
donor_col <- choose_donor_col(meta, cfg$donor_col)
dx_col <- choose_dx_col(meta, cfg$dx_col)

if (is.null(donor_col)) stop("Could not identify donor column.")
if (is.null(dx_col)) stop("Could not identify diagnosis column.")

cluster_levels <- unique(as.character(obj$cluster_use))
num_try <- suppressWarnings(as.numeric(cluster_levels))
if (all(!is.na(num_try))) {
  cluster_levels <- as.character(sort(num_try))
} else {
  cluster_levels <- sort(cluster_levels)
}
obj$cluster_use <- factor(obj$cluster_use, levels = cluster_levels)

audit_meta <- data.frame(
  item = c("n_cells","n_features","default_assay","cluster_col","donor_col","dx_col"),
  value = c(ncol(obj), nrow(obj), DefaultAssay(obj), cluster_col, donor_col, dx_col),
  stringsAsFactors = FALSE
)
fwrite(audit_meta, file.path(cfg$outdir, "tables", "00_package1b_metadata.tsv"), sep = "\t")

# ------------------------------------------------------------
# 2. Load or compute markers
# ------------------------------------------------------------
if (!is.na(cfg$marker_file) && file.exists(cfg$marker_file)) {
  msg("Reading external marker table ...")
  markers <- fread(cfg$marker_file)
} else {
  msg("Running FindAllMarkers ...")
  markers <- FindAllMarkers(
    object = obj,
    only.pos = TRUE,
    min.pct = cfg$min_pct,
    logfc.threshold = cfg$logfc
  )
}

if (!"gene" %in% colnames(markers)) markers$gene <- rownames(markers)
fc_col <- get_fc_col(markers)
if (is.null(fc_col)) stop("Could not identify FC column in marker table.")

markers <- as.data.frame(markers) %>%
  arrange(cluster, desc(.data[[fc_col]]), p_val_adj)

fwrite(markers, file.path(cfg$outdir, "tables", "01_allcluster_markers.tsv.gz"), sep = "\t")

top10_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = .data[[fc_col]], n = 10, with_ties = FALSE) %>%
  ungroup()
fwrite(top10_markers, file.path(cfg$outdir, "tables", "02_allcluster_markers_top10.tsv.gz"), sep = "\t")

# ------------------------------------------------------------
# 3. Average expression by cluster
# ------------------------------------------------------------
msg("Computing AverageExpression by cluster ...")
avg_list <- AverageExpression(
  obj,
  assays = DefaultAssay(obj),
  group.by = "cluster_use",
  slot = "data",
  verbose = FALSE
)
avg_mat <- avg_list[[DefaultAssay(obj)]]
avg_mat <- avg_mat[, cluster_levels, drop = FALSE]

# ------------------------------------------------------------
# 4. Lineage scores
# ------------------------------------------------------------
msg("Computing lineage scores ...")
lineage_score_list <- list()
genes_used_export <- list()

for (lin in lineage_names) {
  res <- get_panel_scores(
    mat = avg_mat,
    genes_anchor = panel_def[[lin]]$anchor,
    genes_support = panel_def[[lin]]$support
  )

  lineage_score_list[[lin]] <- data.frame(
    cluster = cluster_levels,
    lineage = lin,
    anchor_score = as.numeric(res$anchor_score[cluster_levels]),
    support_score = as.numeric(res$support_score[cluster_levels]),
    total_score = as.numeric(res$total_score[cluster_levels]),
    n_anchor_genes_used = length(res$anchor_genes_used),
    n_support_genes_used = length(res$support_genes_used),
    stringsAsFactors = FALSE
  )

  genes_used_export[[lin]] <- data.frame(
    lineage = lin,
    tier = c(rep("anchor_used", length(res$anchor_genes_used)),
             rep("support_used", length(res$support_genes_used))),
    gene = c(res$anchor_genes_used, res$support_genes_used),
    stringsAsFactors = FALSE
  )
}

lineage_score_long <- rbindlist(lineage_score_list, fill = TRUE)
fwrite(lineage_score_long, file.path(cfg$outdir, "tables", "03_lineage_scores_long.tsv"), sep = "\t")

genes_used_df <- rbindlist(genes_used_export, fill = TRUE)
fwrite(genes_used_df, file.path(cfg$outdir, "tables", "04_lineage_genes_used.tsv"), sep = "\t")

lineage_score_wide <- Reduce(function(x, y) merge(x, y, by = "cluster", all = TRUE),
                             lapply(lineage_names, function(lin) {
                               df <- lineage_score_long[lineage_score_long$lineage == lin, c("cluster","anchor_score","support_score","total_score")]
                               colnames(df) <- c(
                                 "cluster",
                                 paste0(lin, "_anchor_score"),
                                 paste0(lin, "_support_score"),
                                 paste0(lin, "_total_score")
                               )
                               df
                             }))
fwrite(lineage_score_wide, file.path(cfg$outdir, "tables", "05_lineage_scores_wide.tsv"), sep = "\t")

# ------------------------------------------------------------
# 5. Marker-hit audit by cluster
# ------------------------------------------------------------
msg("Computing marker-hit audit ...")
cluster_marker_hit_list <- lapply(cluster_levels, function(cl) {
  subm <- markers %>% filter(cluster == cl) %>%
    arrange(desc(.data[[fc_col]]), p_val_adj)

  top50 <- head(unique(subm$gene), cfg$top_n_markers)
  top10 <- head(unique(subm$gene), 10)

  row <- data.frame(
    cluster = cl,
    top10_preview = paste(top10, collapse = ", "),
    stringsAsFactors = FALSE
  )

  for (lin in lineage_names) {
    row[[paste0("n_", lin, "_anchor_hits_top50")]] <- sum(top50 %in% panel_def[[lin]]$anchor)
    row[[paste0("n_", lin, "_support_hits_top50")]] <- sum(top50 %in% panel_def[[lin]]$support)
  }

  row
})
cluster_marker_hits <- rbindlist(cluster_marker_hit_list, fill = TRUE)
fwrite(cluster_marker_hits, file.path(cfg$outdir, "tables", "06_cluster_marker_hits_top50.tsv"), sep = "\t")

# ------------------------------------------------------------
# 6. Cell numbers and donor-level abundance
# ------------------------------------------------------------
msg("Computing donor-level abundance ...")
cell_count_df <- obj@meta.data %>%
  mutate(cluster = as.character(cluster_use)) %>%
  count(cluster, name = "n_cells") %>%
  mutate(frac_cells = n_cells / sum(n_cells))
fwrite(cell_count_df, file.path(cfg$outdir, "tables", "07_cluster_cell_counts.tsv"), sep = "\t")

meta2 <- obj@meta.data
meta2$cell_id <- colnames(obj)
meta2$donor <- as.character(meta2[[donor_col]])
meta2$dx <- normalize_dx(meta2[[dx_col]])
meta2$cluster <- as.character(meta2$cluster_use)

prop_df <- meta2 %>%
  count(donor, dx, cluster, name = "n_cluster") %>%
  left_join(meta2 %>% count(donor, dx, name = "n_total"), by = c("donor","dx")) %>%
  mutate(prop = n_cluster / n_total) %>%
  arrange(cluster, dx, donor)

abund_stats <- lapply(cluster_levels, function(cl) {
  sub <- prop_df %>% filter(cluster == cl)
  data.frame(
    cluster = cl,
    n_donors = length(unique(sub$donor)),
    mean_prop_ASD = mean(sub$prop[sub$dx == "ASD"], na.rm = TRUE),
    mean_prop_Control = mean(sub$prop[sub$dx == "Control"], na.rm = TRUE),
    effect_ASD_minus_Control = mean_diff_asd_minus_ctrl(sub, "prop", "dx"),
    wilcox_p = calc_wilcox_p(sub, "prop", "dx"),
    stringsAsFactors = FALSE
  )
})
abund_stats_df <- rbindlist(abund_stats, fill = TRUE)
abund_stats_df$wilcox_fdr <- p.adjust(abund_stats_df$wilcox_p, method = "BH")
fwrite(prop_df, file.path(cfg$outdir, "tables", "08_cluster_donor_proportions.tsv.gz"), sep = "\t")
fwrite(abund_stats_df, file.path(cfg$outdir, "tables", "09_cluster_donor_abundance_stats.tsv"), sep = "\t")

# ------------------------------------------------------------
# 7. Merge all evidence and classify clusters
# ------------------------------------------------------------
msg("Building annotation table ...")
anno_df <- lineage_score_wide %>%
  left_join(cluster_marker_hits, by = "cluster") %>%
  left_join(cell_count_df, by = "cluster") %>%
  left_join(abund_stats_df, by = "cluster")

anno_list <- lapply(seq_len(nrow(anno_df)), function(i) {
  x <- classify_cluster_row(anno_df[i, , drop = FALSE], cfg)
  data.frame(
    cluster = anno_df$cluster[i],
    suggested_superclass = x$suggested_superclass,
    suggested_lineage = x$suggested_lineage,
    suggested_label = x$suggested_label,
    keep_strict = x$keep_strict,
    keep_lenient = x$keep_lenient,
    rationale = x$rationale,
    stringsAsFactors = FALSE
  )
})
anno_call_df <- rbindlist(anno_list, fill = TRUE)

final_anno <- anno_df %>%
  left_join(anno_call_df, by = "cluster") %>%
  arrange(
    factor(suggested_superclass,
           levels = c("microglia","mixed_microglia","uncertain","contaminant")),
    desc(n_cells)
  )

fwrite(final_anno, file.path(cfg$outdir, "tables", "10_cluster_annotation_contaminant_audit.tsv"), sep = "\t")

# save strict/lenient lists
strict_clusters <- final_anno %>% filter(keep_strict) %>% pull(cluster) %>% as.character()
lenient_clusters <- final_anno %>% filter(keep_lenient) %>% pull(cluster) %>% as.character()
contam_clusters <- final_anno %>% filter(suggested_superclass == "contaminant") %>% pull(cluster) %>% as.character()

writeLines(strict_clusters, file.path(cfg$outdir, "tables", "11_strict_microglia_clusters.txt"))
writeLines(lenient_clusters, file.path(cfg$outdir, "tables", "12_lenient_microglia_clusters.txt"))
writeLines(contam_clusters, file.path(cfg$outdir, "tables", "13_contaminant_clusters.txt"))

# ------------------------------------------------------------
# 8. Add metadata and save objects
# ------------------------------------------------------------
msg("Saving annotated objects ...")
cluster_map <- final_anno %>%
  select(cluster, suggested_superclass, suggested_lineage, suggested_label, keep_strict, keep_lenient)

obj$cluster_use_chr <- as.character(obj$cluster_use)
obj@meta.data <- obj@meta.data %>%
  left_join(cluster_map, by = c("cluster_use_chr" = "cluster"))

saveRDS(obj, file.path(cfg$outdir, "objects", "Package1b_microglia_object_annotated.rds"))

strict_cells <- colnames(obj)[obj$keep_strict %in% TRUE]
lenient_cells <- colnames(obj)[obj$keep_lenient %in% TRUE]

obj_strict <- subset(obj, cells = strict_cells)
obj_lenient <- subset(obj, cells = lenient_cells)

saveRDS(obj_strict, file.path(cfg$outdir, "objects", "Package1b_microglia_object_strict.rds"))
saveRDS(obj_lenient, file.path(cfg$outdir, "objects", "Package1b_microglia_object_lenient.rds"))

# ------------------------------------------------------------
# 9. Recompute donor-level proportions within strict/lenient retained sets
# ------------------------------------------------------------
compute_retained_prop <- function(obj_sub, donor_col, dx_col, out_prefix) {
  md <- obj_sub@meta.data
  md$donor <- as.character(md[[donor_col]])
  md$dx <- normalize_dx(md[[dx_col]])
  md$cluster <- as.character(md$cluster_use)

  prop_df_sub <- md %>%
    count(donor, dx, cluster, name = "n_cluster") %>%
    left_join(md %>% count(donor, dx, name = "n_total"), by = c("donor","dx")) %>%
    mutate(prop = n_cluster / n_total) %>%
    arrange(cluster, dx, donor)

  stats_sub <- lapply(unique(as.character(md$cluster)), function(cl) {
    sub <- prop_df_sub %>% filter(cluster == cl)
    data.frame(
      cluster = cl,
      n_donors = length(unique(sub$donor)),
      mean_prop_ASD = mean(sub$prop[sub$dx == "ASD"], na.rm = TRUE),
      mean_prop_Control = mean(sub$prop[sub$dx == "Control"], na.rm = TRUE),
      effect_ASD_minus_Control = mean_diff_asd_minus_ctrl(sub, "prop", "dx"),
      wilcox_p = calc_wilcox_p(sub, "prop", "dx"),
      stringsAsFactors = FALSE
    )
  })
  stats_sub <- rbindlist(stats_sub, fill = TRUE)
  stats_sub$wilcox_fdr <- p.adjust(stats_sub$wilcox_p, method = "BH")

  fwrite(prop_df_sub, file.path(cfg$outdir, "tables", paste0(out_prefix, "_cluster_donor_proportions.tsv.gz")), sep = "\t")
  fwrite(stats_sub, file.path(cfg$outdir, "tables", paste0(out_prefix, "_cluster_donor_abundance_stats.tsv")), sep = "\t")

  prop_df_sub$cluster <- factor(prop_df_sub$cluster, levels = sort(unique(as.character(prop_df_sub$cluster))))
  p <- ggplot(prop_df_sub, aes(x = dx, y = prop, fill = dx)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.75) +
    geom_jitter(width = 0.12, size = 1.2, alpha = 0.8) +
    facet_wrap(~ cluster, scales = "free_y") +
    theme_classic(base_size = 12) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "none",
      axis.text.x = element_text(angle = 20, hjust = 1)
    ) +
    labs(title = paste0(out_prefix, ": donor-level cluster proportions"), x = NULL, y = "Proportion within retained cells")

  ggsave(file.path(cfg$outdir, "plots", paste0(out_prefix, "_cluster_donor_proportions_boxplot.pdf")), p, width = 12, height = 7)
  ggsave(file.path(cfg$outdir, "plots", paste0(out_prefix, "_cluster_donor_proportions_boxplot.png")), p, width = 12, height = 7, dpi = 300)
}

if (length(strict_cells) > 0) {
  compute_retained_prop(obj_strict, donor_col, dx_col, "14_strict")
}
if (length(lenient_cells) > 0) {
  compute_retained_prop(obj_lenient, donor_col, dx_col, "15_lenient")
}

# ------------------------------------------------------------
# 10. Plots
# ------------------------------------------------------------
msg("Generating plots ...")

# original UMAP by cluster
if ("umap" %in% names(obj@reductions)) {
  pal_cluster <- setNames(grDevices::hcl.colors(length(cluster_levels), "Dynamic"), cluster_levels)

  p_umap_cluster <- DimPlot(
    obj,
    reduction = "umap",
    group.by = "cluster_use",
    label = TRUE,
    repel = TRUE,
    raster = FALSE
  ) +
    scale_color_manual(values = pal_cluster) +
    theme_classic(base_size = 12) +
    labs(title = "Original reclustered object: clusters")

  ggsave(file.path(cfg$outdir, "plots", "16_umap_original_clusters.pdf"), p_umap_cluster, width = 9, height = 7)
  ggsave(file.path(cfg$outdir, "plots", "16_umap_original_clusters.png"), p_umap_cluster, width = 9, height = 7, dpi = 300)

  # UMAP by audit superclass
  audit_levels <- c("microglia","mixed_microglia","uncertain","contaminant")
  obj$suggested_superclass <- factor(obj$suggested_superclass, levels = audit_levels)
  audit_cols <- c(
    microglia = "#1b9e77",
    mixed_microglia = "#7570b3",
    uncertain = "#999999",
    contaminant = "#d95f02"
  )

  p_umap_superclass <- DimPlot(
    obj,
    reduction = "umap",
    group.by = "suggested_superclass",
    label = FALSE,
    raster = FALSE
  ) +
    scale_color_manual(values = audit_cols, drop = FALSE) +
    theme_classic(base_size = 12) +
    labs(title = "Cluster audit superclass")

  ggsave(file.path(cfg$outdir, "plots", "17_umap_audit_superclass.pdf"), p_umap_superclass, width = 8, height = 7)
  ggsave(file.path(cfg$outdir, "plots", "17_umap_audit_superclass.png"), p_umap_superclass, width = 8, height = 7, dpi = 300)

  # UMAP by suggested label
  p_umap_label <- DimPlot(
    obj,
    reduction = "umap",
    group.by = "suggested_label",
    label = FALSE,
    raster = FALSE
  ) +
    theme_classic(base_size = 10) +
    labs(title = "Cluster audit detailed label")

  ggsave(file.path(cfg$outdir, "plots", "18_umap_audit_label.pdf"), p_umap_label, width = 12, height = 10)
  ggsave(file.path(cfg$outdir, "plots", "18_umap_audit_label.png"), p_umap_label, width = 12, height = 10, dpi = 300)
}

# lineage score heatmap-like tile
score_heat <- lineage_score_long %>%
  mutate(
    lineage = factor(lineage, levels = lineage_names),
    cluster = factor(cluster, levels = cluster_levels)
  )

p_score <- ggplot(score_heat, aes(x = lineage, y = cluster, fill = total_score)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0) +
  theme_classic(base_size = 12) +
  labs(title = "Lineage total scores by cluster", x = NULL, y = "Cluster") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

ggsave(file.path(cfg$outdir, "plots", "19_lineage_score_heatmap.pdf"), p_score, width = 10, height = 6.5)
ggsave(file.path(cfg$outdir, "plots", "19_lineage_score_heatmap.png"), p_score, width = 10, height = 6.5, dpi = 300)

# lineage marker dotplot
genes_dot <- intersect(lineage_marker_panel, rownames(obj))
if (length(genes_dot) > 0) {
  p_dot <- DotPlot(obj, features = genes_dot, group.by = "cluster_use") +
    RotatedAxis() +
    theme_classic(base_size = 10) +
    labs(title = "Transparent lineage marker panel across clusters")
  ggsave(file.path(cfg$outdir, "plots", "20_lineage_marker_dotplot.pdf"), p_dot, width = 15, height = 6)
  ggsave(file.path(cfg$outdir, "plots", "20_lineage_marker_dotplot.png"), p_dot, width = 15, height = 6, dpi = 300)
}

# barplot of cluster class
class_plot_df <- final_anno %>%
  count(suggested_superclass) %>%
  mutate(suggested_superclass = factor(suggested_superclass, levels = c("microglia","mixed_microglia","uncertain","contaminant")))

p_class <- ggplot(class_plot_df, aes(x = suggested_superclass, y = n, fill = suggested_superclass)) +
  geom_col() +
  scale_fill_manual(values = c(
    microglia = "#1b9e77",
    mixed_microglia = "#7570b3",
    uncertain = "#999999",
    contaminant = "#d95f02"
  ), drop = FALSE) +
  theme_classic(base_size = 12) +
  labs(title = "Number of clusters by audit superclass", x = NULL, y = "Number of clusters") +
  theme(legend.position = "none")

ggsave(file.path(cfg$outdir, "plots", "21_cluster_count_by_superclass.pdf"), p_class, width = 6.5, height = 4.5)
ggsave(file.path(cfg$outdir, "plots", "21_cluster_count_by_superclass.png"), p_class, width = 6.5, height = 4.5, dpi = 300)

# concise summary table for manuscript / rebuttal
summary_export <- final_anno %>%
  select(
    cluster, n_cells, frac_cells,
    mean_prop_ASD, mean_prop_Control, effect_ASD_minus_Control, wilcox_p, wilcox_fdr,
    microglia_anchor_score, microglia_support_score, microglia_total_score,
    suggested_superclass, suggested_lineage, suggested_label,
    top10_preview, rationale
  )
fwrite(summary_export, file.path(cfg$outdir, "tables", "22_cluster_annotation_summary_for_manuscript.tsv"), sep = "\t")

save_session_info(file.path(cfg$outdir, "sessionInfo_package1b.txt"))

manifest <- data.frame(
  output_file = c(
    "10_cluster_annotation_contaminant_audit.tsv",
    "11_strict_microglia_clusters.txt",
    "12_lenient_microglia_clusters.txt",
    "13_contaminant_clusters.txt",
    "Package1b_microglia_object_annotated.rds",
    "Package1b_microglia_object_strict.rds",
    "Package1b_microglia_object_lenient.rds",
    "17_umap_audit_superclass",
    "19_lineage_score_heatmap",
    "20_lineage_marker_dotplot",
    "22_cluster_annotation_summary_for_manuscript.tsv"
  ),
  description = c(
    "full cluster audit table with suggested class and rationale",
    "strictly retained microglia clusters",
    "leniently retained microglia/mixed clusters",
    "clusters classified as contaminants",
    "annotated original object",
    "strict retained object",
    "lenient retained object",
    "UMAP colored by audit superclass",
    "cluster-by-lineage score heatmap",
    "transparent lineage marker panel dotplot",
    "compact summary table for rebuttal/manuscript"
  ),
  stringsAsFactors = FALSE
)
fwrite(manifest, file.path(cfg$outdir, "PACKAGE1b_manifest.tsv"), sep = "\t")

msg("Package 1b finished.")
