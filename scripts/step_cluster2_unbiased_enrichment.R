#!/usr/bin/env Rscript

.libPaths(c(
  "/gpfs/hpc/home/lijc/lianaoj/R/x86_64-conda-linux-gnu-library/4.2",
  "/gpfs/hpc/home/lijc/lianaoj/mambaforge/envs/scRNA3/lib/R/library"
))

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(msigdbr)
})

options(stringsAsFactors = FALSE)

cfg <- list(
  microglia_rds   = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/PsychENCODE_Science_2024/Step3_Microglia/Microglia_Clustered.rds",
  outdir          = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Cluster2_Unbiased_Enrichment_light",
  assay           = "RNA",
  cluster_col     = NULL,
  target_cluster  = "2",

  ## DE thresholds
  min_pct         = 0.05,
  logfc_threshold = 0.0,
  up_logfc_cutoff = 0.25,
  down_logfc_cutoff = -0.25,
  padj_cutoff     = 0.05,

  ## gene set size filters
  min_set_size    = 10,
  max_set_size    = 1000,

  ## plotting
  top_show        = 20
)

dir.create(cfg$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cfg$outdir, "plots"), recursive = TRUE, showWarnings = FALSE)

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

choose_cluster_col <- function(obj, user_col = NULL) {
  meta <- obj@meta.data
  if (!is.null(user_col) && user_col %in% names(meta)) return(user_col)
  if ("seurat_clusters" %in% names(meta)) return("seurat_clusters")
  cand <- find_first_col(meta, c("^cluster$", "cluster", "seurat", "res\\.", "RNA_snn", "wsnn"))
  if (!is.null(cand)) return(cand)
  obj$cluster_auto_from_idents <- as.character(Idents(obj))
  return("cluster_auto_from_idents")
}

get_fc_col <- function(markers) {
  fc_candidates <- c("avg_log2FC", "avg_logFC", "log2FC", "logFC")
  hit <- fc_candidates[fc_candidates %in% colnames(markers)]
  if (length(hit) == 0) return(NULL)
  hit[1]
}

ensure_pkg_loadable <- function(pkg) {
  ok <- tryCatch({
    suppressPackageStartupMessages(
      library(pkg, character.only = TRUE, quietly = TRUE, warn.conflicts = FALSE)
    )
    TRUE
  }, error = function(e) {
    message("Package load failed: ", pkg, " -> ", conditionMessage(e))
    FALSE
  })
  if (!ok) stop("Package is installed but not loadable: ", pkg)
}

safe_write_table <- function(x, file) {
  if (is.null(x) || nrow(x) == 0) {
    fwrite(data.frame(message = "No significant result"), file, sep = "\t")
  } else {
    fwrite(as.data.frame(x), file, sep = "\t")
  }
}

clean_term_name <- function(x) {
  x <- gsub("^HALLMARK_", "", x)
  x <- gsub("^REACTOME_", "", x)
  x <- gsub("^GOBP_", "", x)
  x <- gsub("^GOBP_", "", x)
  x <- gsub("^GO_", "", x)
  x <- gsub("_", " ", x)
  x
}

plot_ora_bar <- function(df, title_text, out_pdf, top_n = 20) {
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  plot_df <- as.data.frame(df)

  needed <- c("display_name", "p_adj")
  if (!all(needed %in% colnames(plot_df))) return(invisible(NULL))

  plot_df <- plot_df %>%
    filter(!is.na(p_adj), overlap_size > 0) %>%
    arrange(p_adj, desc(overlap_size)) %>%
    slice_head(n = top_n) %>%
    mutate(
      neglog10padj = -log10(p_adj),
      display_name = factor(display_name, levels = rev(display_name))
    )

  if (nrow(plot_df) == 0) return(invisible(NULL))

  p <- ggplot(plot_df, aes(x = neglog10padj, y = display_name)) +
    geom_col() +
    theme_classic(base_size = 12) +
    labs(
      title = title_text,
      x = "-log10(adj.P)",
      y = NULL
    )

  ggsave(out_pdf, p, width = 9, height = max(4.5, 0.28 * nrow(plot_df) + 1.5))
}

plot_rank_bar <- function(df, title_text, out_pdf, top_n = 20) {
  if (is.null(df) || nrow(df) == 0) return(invisible(NULL))
  plot_df <- as.data.frame(df)

  needed <- c("display_name", "NES_like", "p_adj")
  if (!all(needed %in% colnames(plot_df))) return(invisible(NULL))

  plot_df <- plot_df %>%
    arrange(p_adj, desc(abs(NES_like))) %>%
    slice_head(n = top_n) %>%
    mutate(display_name = factor(display_name, levels = rev(display_name)))

  if (nrow(plot_df) == 0) return(invisible(NULL))

  p <- ggplot(plot_df, aes(x = NES_like, y = display_name)) +
    geom_col() +
    theme_classic(base_size = 12) +
    labs(
      title = title_text,
      x = "NES-like rank statistic",
      y = NULL
    )

  ggsave(out_pdf, p, width = 9, height = max(4.5, 0.28 * nrow(plot_df) + 1.5))
}

get_msigdb_all <- function() {
  x <- as.data.frame(msigdbr::msigdbr(species = "Homo sapiens"))

  term_col <- find_first_col(x, c("^gs_name$", "^term$", "geneset"))
  gene_col <- find_first_col(x, c("^gene_symbol$", "human_gene_symbol", "symbol"))
  coll_col <- find_first_col(x, c("^gs_collection$", "^collection$", "^gs_cat$", "^category$"))
  subcoll_col <- find_first_col(x, c("^gs_subcollection$", "^subcollection$", "^gs_subcat$", "^subcategory$"))
  desc_col <- find_first_col(x, c("^gs_description$", "description"))

  if (is.null(term_col) || is.null(gene_col)) {
    stop("Could not identify term or gene column in msigdbr output.")
  }

  list(
    df = x,
    term_col = term_col,
    gene_col = gene_col,
    coll_col = coll_col,
    subcoll_col = subcoll_col,
    desc_col = desc_col
  )
}

extract_gene_sets <- function(msig, type = c("hallmark", "reactome", "go_bp")) {
  type <- match.arg(type)
  df <- msig$df
  term_col <- msig$term_col
  gene_col <- msig$gene_col
  coll_col <- msig$coll_col
  subcoll_col <- msig$subcoll_col
  desc_col <- msig$desc_col

  term_upper <- toupper(df[[term_col]])
  coll_upper <- if (!is.null(coll_col)) toupper(df[[coll_col]]) else rep(NA_character_, nrow(df))
  subcoll_upper <- if (!is.null(subcoll_col)) toupper(df[[subcoll_col]]) else rep(NA_character_, nrow(df))

  keep <- rep(FALSE, nrow(df))

  if (type == "hallmark") {
    keep <- grepl("^HALLMARK_", term_upper)
    if (!all(is.na(coll_upper))) {
      keep <- keep | (coll_upper == "H")
    }
  }

  if (type == "reactome") {
    keep <- grepl("^REACTOME_", term_upper)
    if (!all(is.na(coll_upper))) {
      keep <- keep | (coll_upper == "C2" & grepl("REACTOME", subcoll_upper))
    }
  }

  if (type == "go_bp") {
    keep <- grepl("^GOBP_", term_upper) | grepl("^GOBP_", term_upper)
    if (!all(is.na(coll_upper))) {
      keep <- keep | (coll_upper == "C5" & grepl("GO:BP|GOBP|BP", subcoll_upper))
    }
  }

  out <- df[keep, , drop = FALSE]
  if (nrow(out) == 0) return(NULL)

  out2 <- data.frame(
    term_name = out[[term_col]],
    gene_symbol = out[[gene_col]],
    stringsAsFactors = FALSE
  )

  if (!is.null(desc_col) && desc_col %in% colnames(out)) {
    out2$term_description <- out[[desc_col]]
  } else {
    out2$term_description <- out2$term_name
  }

  out2 <- out2[!is.na(out2$gene_symbol) & out2$gene_symbol != "", , drop = FALSE]
  out2 <- unique(out2)
  out2
}

run_ora <- function(query_genes, universe_genes, term2gene,
                    min_set_size = 10, max_set_size = 1000) {

  query_genes <- unique(intersect(query_genes, universe_genes))
  universe_genes <- unique(universe_genes)

  if (length(query_genes) < 3) return(NULL)

  term2gene <- unique(term2gene[, c("term_name", "gene_symbol", "term_description"), drop = FALSE])
  term2gene <- term2gene[term2gene$gene_symbol %in% universe_genes, , drop = FALSE]

  gs_list <- split(term2gene$gene_symbol, term2gene$term_name)
  desc_map <- unique(term2gene[, c("term_name", "term_description"), drop = FALSE])

  N <- length(universe_genes)
  n <- length(query_genes)

  res_list <- lapply(names(gs_list), function(term) {
    geneset <- unique(gs_list[[term]])
    K <- length(geneset)

    if (K < min_set_size || K > max_set_size) return(NULL)

    overlap <- intersect(query_genes, geneset)
    k <- length(overlap)

    if (k == 0) return(NULL)

    pval <- phyper(q = k - 1, m = K, n = N - K, k = n, lower.tail = FALSE)
    fold_enrichment <- (k / n) / (K / N)

    data.frame(
      term_name = term,
      set_size = K,
      query_size = n,
      universe_size = N,
      overlap_size = k,
      fold_enrichment = fold_enrichment,
      p_value = pval,
      overlap_genes = paste(sort(overlap), collapse = ","),
      stringsAsFactors = FALSE
    )
  })

  res <- bind_rows(res_list)
  if (nrow(res) == 0) return(NULL)

  res <- res %>%
    left_join(desc_map, by = "term_name") %>%
    mutate(
      p_adj = p.adjust(p_value, method = "BH"),
      display_name = clean_term_name(term_name)
    ) %>%
    arrange(p_adj, p_value, desc(overlap_size), desc(fold_enrichment))

  res
}

run_rank_enrichment <- function(ranked_vec, term2gene,
                                min_set_size = 10, max_set_size = 1000) {

  ranked_vec <- ranked_vec[!is.na(ranked_vec)]
  ranked_vec <- ranked_vec[!duplicated(names(ranked_vec))]
  ranked_vec <- sort(ranked_vec, decreasing = TRUE)

  genes <- names(ranked_vec)
  N <- length(genes)
  if (N < 100) return(NULL)

  term2gene <- unique(term2gene[, c("term_name", "gene_symbol", "term_description"), drop = FALSE])
  term2gene <- term2gene[term2gene$gene_symbol %in% genes, , drop = FALSE]

  gs_list <- split(term2gene$gene_symbol, term2gene$term_name)
  desc_map <- unique(term2gene[, c("term_name", "term_description"), drop = FALSE])

  ## rank: smaller mean rank => enriched near top
  rk <- rank(-ranked_vec, ties.method = "average")
  expected_mean_rank <- (N + 1) / 2

  res_list <- lapply(names(gs_list), function(term) {
    gs <- unique(gs_list[[term]])
    gs <- intersect(gs, genes)
    m <- length(gs)

    if (m < min_set_size || m > max_set_size) return(NULL)

    idx <- genes %in% gs
    mean_rank_in <- mean(rk[idx])

    ## approximate z under random sampling without replacement
    var_mean_rank <- ((N - m) * (N + 1)) / (12 * m)
    if (var_mean_rank <= 0) return(NULL)

    z <- (expected_mean_rank - mean_rank_in) / sqrt(var_mean_rank)
    pval <- 2 * pnorm(-abs(z))

    leading_top <- head(genes[idx][order(match(genes[idx], genes))], 30)
    leading_bottom <- head(rev(genes[idx][order(match(genes[idx], genes))]), 30)

    data.frame(
      term_name = term,
      set_size = m,
      mean_rank = mean_rank_in,
      NES_like = z,
      direction = ifelse(z >= 0, "top_enriched", "bottom_enriched"),
      p_value = pval,
      leading_genes_top30 = paste(leading_top, collapse = ","),
      trailing_genes_bottom30 = paste(leading_bottom, collapse = ","),
      stringsAsFactors = FALSE
    )
  })

  res <- bind_rows(res_list)
  if (nrow(res) == 0) return(NULL)

  res <- res %>%
    left_join(desc_map, by = "term_name") %>%
    mutate(
      p_adj = p.adjust(p_value, method = "BH"),
      display_name = clean_term_name(term_name)
    ) %>%
    arrange(p_adj, desc(abs(NES_like)))

  res
}

## --------------------------------------------------
## 0. package check
## --------------------------------------------------
ensure_pkg_loadable("msigdbr")

## --------------------------------------------------
## 1. read object
## --------------------------------------------------
msg("Reading microglia object ...")
obj <- readRDS(cfg$microglia_rds)
if (!inherits(obj, "Seurat")) stop("Input RDS is not a Seurat object.")

if (cfg$assay %in% Assays(obj)) {
  DefaultAssay(obj) <- cfg$assay
}

cluster_col <- choose_cluster_col(obj, cfg$cluster_col)
obj$cluster_use <- as.character(obj@meta.data[[cluster_col]])
Idents(obj) <- obj$cluster_use

audit <- data.frame(
  item = c("n_cells", "n_features", "default_assay", "cluster_col", "target_cluster"),
  value = c(ncol(obj), nrow(obj), DefaultAssay(obj), cluster_col, cfg$target_cluster)
)
fwrite(audit, file.path(cfg$outdir, "Audit.tsv"), sep = "\t")

cluster_levels <- sort(unique(as.character(obj$cluster_use)))
writeLines(cluster_levels, file.path(cfg$outdir, "Available_Clusters.txt"))

if (!(cfg$target_cluster %in% cluster_levels)) {
  stop("Target cluster ", cfg$target_cluster, " not found. Available clusters: ",
       paste(cluster_levels, collapse = ", "))
}

## --------------------------------------------------
## 2. DE: target cluster vs others
## --------------------------------------------------
msg("Running FindMarkers for Cluster ", cfg$target_cluster, " vs all other microglia ...")

markers <- FindMarkers(
  object = obj,
  ident.1 = cfg$target_cluster,
  ident.2 = setdiff(cluster_levels, cfg$target_cluster),
  only.pos = FALSE,
  min.pct = cfg$min_pct,
  logfc.threshold = cfg$logfc_threshold
)

markers <- as.data.frame(markers)
markers$gene <- rownames(markers)
fc_col <- get_fc_col(markers)
if (is.null(fc_col)) stop("Could not identify fold-change column in marker table.")

markers <- markers %>%
  arrange(desc(.data[[fc_col]]), p_val_adj)

fwrite(markers, file.path(cfg$outdir, "Cluster2_vs_Others_Markers.tsv"), sep = "\t")

up_genes <- markers %>%
  filter(.data[[fc_col]] >= cfg$up_logfc_cutoff, p_val_adj < cfg$padj_cutoff) %>%
  pull(gene) %>%
  unique()

down_genes <- markers %>%
  filter(.data[[fc_col]] <= cfg$down_logfc_cutoff, p_val_adj < cfg$padj_cutoff) %>%
  pull(gene) %>%
  unique()

universe_genes <- unique(markers$gene)

writeLines(up_genes, file.path(cfg$outdir, "Cluster2_Up_Genes.txt"))
writeLines(down_genes, file.path(cfg$outdir, "Cluster2_Down_Genes.txt"))
writeLines(universe_genes, file.path(cfg$outdir, "Cluster2_Universe_Genes.txt"))

de_summary <- data.frame(
  metric = c("n_up_genes", "n_down_genes", "n_universe_genes"),
  value = c(length(up_genes), length(down_genes), length(universe_genes))
)
fwrite(de_summary, file.path(cfg$outdir, "DE_Summary.tsv"), sep = "\t")

ranked_tbl <- markers[, c("gene", fc_col)]
ranked_tbl <- ranked_tbl[!duplicated(ranked_tbl$gene), , drop = FALSE]
ranked_vec <- ranked_tbl[[fc_col]]
names(ranked_vec) <- ranked_tbl$gene
ranked_vec <- sort(ranked_vec, decreasing = TRUE)
saveRDS(ranked_vec, file.path(cfg$outdir, "Cluster2_ranked_geneList.rds"))

## --------------------------------------------------
## 3. local MSigDB gene sets
## --------------------------------------------------
msg("Loading MSigDB via msigdbr ...")
msig <- get_msigdb_all()

hallmark_sets <- extract_gene_sets(msig, "hallmark")
reactome_sets <- extract_gene_sets(msig, "reactome")
gobp_sets <- extract_gene_sets(msig, "go_bp")

safe_write_table(hallmark_sets, file.path(cfg$outdir, "MSigDB_Hallmark_term2gene.tsv"))
safe_write_table(reactome_sets, file.path(cfg$outdir, "MSigDB_Reactome_term2gene.tsv"))
safe_write_table(gobp_sets, file.path(cfg$outdir, "MSigDB_GOBP_term2gene.tsv"))

## --------------------------------------------------
## 4. enrichment
## --------------------------------------------------
msg("Running GO BP ORA (local MSigDB C5 GO:BP) ...")
go_up <- run_ora(up_genes, universe_genes, gobp_sets,
                 min_set_size = cfg$min_set_size,
                 max_set_size = cfg$max_set_size)

go_down <- run_ora(down_genes, universe_genes, gobp_sets,
                   min_set_size = cfg$min_set_size,
                   max_set_size = cfg$max_set_size)

safe_write_table(go_up, file.path(cfg$outdir, "Cluster2_GO_BP_Up_localMSigDB.tsv"))
safe_write_table(go_down, file.path(cfg$outdir, "Cluster2_GO_BP_Down_localMSigDB.tsv"))

msg("Running Reactome ORA (local MSigDB Reactome) ...")
react_up <- run_ora(up_genes, universe_genes, reactome_sets,
                    min_set_size = cfg$min_set_size,
                    max_set_size = cfg$max_set_size)

react_down <- run_ora(down_genes, universe_genes, reactome_sets,
                      min_set_size = cfg$min_set_size,
                      max_set_size = cfg$max_set_size)

safe_write_table(react_up, file.path(cfg$outdir, "Cluster2_Reactome_Up_localMSigDB.tsv"))
safe_write_table(react_down, file.path(cfg$outdir, "Cluster2_Reactome_Down_localMSigDB.tsv"))

msg("Running Hallmark rank-based enrichment (local MSigDB Hallmark) ...")
hallmark_rank <- run_rank_enrichment(ranked_vec, hallmark_sets,
                                     min_set_size = cfg$min_set_size,
                                     max_set_size = cfg$max_set_size)

safe_write_table(hallmark_rank, file.path(cfg$outdir, "Cluster2_Hallmark_RankEnrichment.tsv"))

## --------------------------------------------------
## 5. plots
## --------------------------------------------------
msg("Saving plots ...")

plot_ora_bar(
  df = go_up,
  title_text = "Cluster 2 vs Others: GO BP ORA (Up genes, local MSigDB)",
  out_pdf = file.path(cfg$outdir, "plots", "Cluster2_GO_BP_Up_localMSigDB.pdf"),
  top_n = cfg$top_show
)

plot_ora_bar(
  df = go_down,
  title_text = "Cluster 2 vs Others: GO BP ORA (Down genes, local MSigDB)",
  out_pdf = file.path(cfg$outdir, "plots", "Cluster2_GO_BP_Down_localMSigDB.pdf"),
  top_n = cfg$top_show
)

plot_ora_bar(
  df = react_up,
  title_text = "Cluster 2 vs Others: Reactome ORA (Up genes, local MSigDB)",
  out_pdf = file.path(cfg$outdir, "plots", "Cluster2_Reactome_Up_localMSigDB.pdf"),
  top_n = cfg$top_show
)

plot_ora_bar(
  df = react_down,
  title_text = "Cluster 2 vs Others: Reactome ORA (Down genes, local MSigDB)",
  out_pdf = file.path(cfg$outdir, "plots", "Cluster2_Reactome_Down_localMSigDB.pdf"),
  top_n = cfg$top_show
)

plot_rank_bar(
  df = hallmark_rank,
  title_text = "Cluster 2 vs Others: Hallmark rank-based enrichment (local MSigDB)",
  out_pdf = file.path(cfg$outdir, "plots", "Cluster2_Hallmark_RankEnrichment.pdf"),
  top_n = cfg$top_show
)

## --------------------------------------------------
## 6. top summaries
## --------------------------------------------------
top_go_up <- if (!is.null(go_up) && nrow(go_up) > 0) {
  go_up %>% arrange(p_adj, desc(overlap_size), desc(fold_enrichment)) %>% slice_head(n = 15)
} else NULL

top_react_up <- if (!is.null(react_up) && nrow(react_up) > 0) {
  react_up %>% arrange(p_adj, desc(overlap_size), desc(fold_enrichment)) %>% slice_head(n = 15)
} else NULL

top_hallmark <- if (!is.null(hallmark_rank) && nrow(hallmark_rank) > 0) {
  hallmark_rank %>% arrange(p_adj, desc(abs(NES_like))) %>% slice_head(n = 15)
} else NULL

safe_write_table(top_go_up, file.path(cfg$outdir, "Top15_GO_BP_Up_localMSigDB.tsv"))
safe_write_table(top_react_up, file.path(cfg$outdir, "Top15_Reactome_Up_localMSigDB.tsv"))
safe_write_table(top_hallmark, file.path(cfg$outdir, "Top15_Hallmark_RankEnrichment.tsv"))

## --------------------------------------------------
## 7. method note
## --------------------------------------------------
method_note <- c(
  "This lightweight script avoids clusterProfiler / ReactomePA / ggtree / enrichplot.",
  "GO BP and Reactome are performed by over-representation analysis (hypergeometric test) using local MSigDB gene sets from msigdbr.",
  "Hallmark is performed by a rank-based enrichment statistic computed from the full ranked gene list; this is not identical to clusterProfiler::GSEA.",
  "Output is intended for robust exploratory interpretation when package dependency conflicts block the original workflow."
)
writeLines(method_note, file.path(cfg$outdir, "Method_Note_lightweight.txt"))

save_session_info(file.path(cfg$outdir, "sessionInfo_Cluster2_Unbiased_Enrichment_light.txt"))
msg("Done.")
