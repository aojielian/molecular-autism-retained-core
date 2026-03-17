#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

options(stringsAsFactors = FALSE)

cfg <- list(
  sig_rds = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/PsychENCODE_Science_2024/Step3_Microglia/Dual_Signatures_List.rds",
  outdir  = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Signature_Audit"
)

dir.create(cfg$outdir, recursive = TRUE, showWarnings = FALSE)

msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ..., "\n", sep = "")
}

safe_chr <- function(x, n = 10) {
  x <- as.character(x)
  x <- x[!is.na(x)]
  x <- unique(x)
  paste(head(x, n), collapse = ", ")
}

is_gene_like_character <- function(x) {
  if (!is.character(x)) return(FALSE)
  x <- x[!is.na(x)]
  x <- x[nzchar(x)]
  if (length(x) == 0) return(FALSE)

  # 粗略判断：大写字母/数字/连字符/点 形式较多时，像基因名
  frac_gene_like <- mean(grepl("^[A-Za-z0-9._-]+$", x))
  uniq_n <- length(unique(x))

  frac_gene_like > 0.7 && uniq_n >= 5
}

extract_gene_vector_from_df <- function(df) {
  cn <- names(df)
  cn_low <- tolower(cn)

  # 优先使用显式 gene 列
  gene_cols <- cn[cn_low %in% c("gene", "genes", "symbol", "gene_symbol", "feature", "features", "external_gene_name")]
  if (length(gene_cols) > 0) {
    g <- unique(as.character(df[[gene_cols[1]]]))
    g <- g[!is.na(g) & nzchar(g)]
    return(g)
  }

  # 再看是否有行名像基因
  if (!is.null(rownames(df)) && !all(grepl("^\\d+$", rownames(df)))) {
    g <- unique(rownames(df))
    g <- g[!is.na(g) & nzchar(g)]
    if (length(g) >= 5) return(g)
  }

  # 再尝试从字符列中找最像基因列
  char_cols <- cn[sapply(df, function(z) is.character(z) || is.factor(z))]
  if (length(char_cols) > 0) {
    best_col <- NULL
    best_n <- 0
    for (cc in char_cols) {
      vec <- as.character(df[[cc]])
      vec <- vec[!is.na(vec) & nzchar(vec)]
      if (length(vec) >= 5 && is_gene_like_character(vec)) {
        if (length(unique(vec)) > best_n) {
          best_n <- length(unique(vec))
          best_col <- cc
        }
      }
    }
    if (!is.null(best_col)) {
      g <- unique(as.character(df[[best_col]]))
      g <- g[!is.na(g) & nzchar(g)]
      return(g)
    }
  }

  character(0)
}

get_preview <- function(x, n = 8) {
  if (is.null(x)) return("NULL")

  if (is.atomic(x) && !is.list(x)) {
    return(safe_chr(x, n = n))
  }

  if (is.data.frame(x)) {
    cn <- names(x)
    return(paste("cols=", paste(head(cn, 10), collapse = ",")))
  }

  if (is.list(x)) {
    nms <- names(x)
    if (is.null(nms)) {
      return(paste0("list(", length(x), ") unnamed"))
    } else {
      return(paste(head(nms, 10), collapse = ", "))
    }
  }

  paste(class(x), collapse = ",")
}

summarize_object <- function(x, path) {
  cls <- paste(class(x), collapse = ";")
  len <- tryCatch(length(x), error = function(e) NA_integer_)
  nr  <- tryCatch(if (is.matrix(x) || is.data.frame(x)) nrow(x) else NA_integer_, error = function(e) NA_integer_)
  nc  <- tryCatch(if (is.matrix(x) || is.data.frame(x)) ncol(x) else NA_integer_, error = function(e) NA_integer_)
  cn  <- tryCatch(if (is.data.frame(x) || is.matrix(x)) paste(head(colnames(x), 20), collapse = ",") else "", error = function(e) "")
  rn_preview <- tryCatch(if (!is.null(rownames(x))) paste(head(rownames(x), 10), collapse = ",") else "", error = function(e) "")

  data.frame(
    path = path,
    class = cls,
    length = len,
    nrow = nr,
    ncol = nc,
    colnames_preview = cn,
    rownames_preview = rn_preview,
    value_preview = get_preview(x, n = 8),
    stringsAsFactors = FALSE
  )
}

flatten_all_objects <- function(x, path = "root", max_depth = 8) {
  out <- list()
  out[[length(out) + 1]] <- summarize_object(x, path)

  if (max_depth <= 0) return(out)

  # list 递归
  if (is.list(x) && !is.data.frame(x)) {
    nms <- names(x)
    if (is.null(nms)) nms <- paste0("[[", seq_along(x), "]]")
    for (i in seq_along(x)) {
      child_path <- paste0(path, " / ", nms[i])
      out <- c(out, flatten_all_objects(x[[i]], child_path, max_depth - 1))
    }
  }

  # data.frame 内部也递归每列
  if (is.data.frame(x)) {
    for (cc in names(x)) {
      child_path <- paste0(path, " / $", cc)
      out <- c(out, flatten_all_objects(x[[cc]], child_path, max_depth - 1))
    }
  }

  out
}

flatten_gene_sets <- function(x, path = "root", max_depth = 8) {
  res <- list()

  if (max_depth <= 0) return(res)
  if (is.null(x)) return(res)

  if (is.factor(x)) x <- as.character(x)

  # 情况1：直接字符向量
  if (is.character(x)) {
    genes <- unique(x)
    genes <- genes[!is.na(genes) & nzchar(genes)]
    if (length(genes) >= 5 && is_gene_like_character(genes)) {
      res[[path]] <- genes
    }
    return(res)
  }

  # 情况2：data.frame
  if (is.data.frame(x)) {
    genes <- extract_gene_vector_from_df(x)
    if (length(genes) >= 5) {
      res[[path]] <- genes
    }
    for (cc in names(x)) {
      child_path <- paste0(path, " / $", cc)
      res <- c(res, flatten_gene_sets(x[[cc]], child_path, max_depth - 1))
    }
    return(res)
  }

  # 情况3：list
  if (is.list(x)) {
    nms <- names(x)
    if (is.null(nms)) nms <- paste0("[[", seq_along(x), "]]")
    for (i in seq_along(x)) {
      child_path <- paste0(path, " / ", nms[i])
      res <- c(res, flatten_gene_sets(x[[i]], child_path, max_depth - 1))
    }
    return(res)
  }

  res
}

pick_keyword_hits <- function(paths, keywords) {
  keep <- rep(FALSE, length(paths))
  for (kw in keywords) {
    keep <- keep | grepl(kw, paths, ignore.case = TRUE)
  }
  keep
}

# --------------------------------------------------
# 1. 读取对象
# --------------------------------------------------
msg("Reading RDS: ", cfg$sig_rds)
obj <- readRDS(cfg$sig_rds)

# 保存 str() 文本
msg("Saving str() output ...")
capture.output(str(obj, max.level = 6), file = file.path(cfg$outdir, "Dual_Signatures_str.txt"))

# 保存顶层信息
top_info <- data.frame(
  item = c("top_class", "top_length", "top_names_preview"),
  value = c(
    paste(class(obj), collapse = ";"),
    tryCatch(length(obj), error = function(e) NA_integer_),
    if (!is.null(names(obj))) paste(head(names(obj), 30), collapse = ", ") else "NULL"
  ),
  stringsAsFactors = FALSE
)
fwrite(top_info, file.path(cfg$outdir, "TopLevel_Info.tsv"), sep = "\t")

# --------------------------------------------------
# 2. 递归展平所有对象
# --------------------------------------------------
msg("Flattening all nested objects ...")
all_objs <- rbindlist(flatten_all_objects(obj, path = "root", max_depth = 8), fill = TRUE)
fwrite(all_objs, file.path(cfg$outdir, "All_Nested_Objects.tsv"), sep = "\t")

# --------------------------------------------------
# 3. 提取所有候选 gene sets
# --------------------------------------------------
msg("Extracting candidate gene sets ...")
gs <- flatten_gene_sets(obj, path = "root", max_depth = 8)

if (length(gs) == 0) {
  msg("No candidate gene sets detected.")
  fwrite(data.frame(message = "No candidate gene sets detected"), file.path(cfg$outdir, "GeneSet_Summary.tsv"), sep = "\t")
} else {
  gs_summary <- data.frame(
    path = names(gs),
    n_genes = sapply(gs, length),
    preview = sapply(gs, function(x) safe_chr(x, 12)),
    stringsAsFactors = FALSE
  )
  gs_summary <- gs_summary[order(-gs_summary$n_genes, gs_summary$path), ]
  fwrite(gs_summary, file.path(cfg$outdir, "GeneSet_Summary.tsv"), sep = "\t")

  # 每个 gene set 单独导出
  dir.create(file.path(cfg$outdir, "gene_sets"), showWarnings = FALSE)
  for (i in seq_along(gs)) {
    nm <- names(gs)[i]
    genes <- unique(gs[[i]])
    genes <- genes[!is.na(genes) & nzchar(genes)]

    safe_name <- gsub("[^A-Za-z0-9._-]+", "_", nm)
    safe_name <- substr(safe_name, 1, 180)

    writeLines(genes, file.path(cfg$outdir, "gene_sets", paste0(sprintf("%03d", i), "_", safe_name, ".txt")))
  }

  # 关键词筛选
  keywords <- c("risk", "homeo", "homeost", "cluster.?0", "cluster.?1", "cluster.?2", "cluster.?4", "spp1", "microgl", "dam")
  hit_idx <- pick_keyword_hits(gs_summary$path, keywords)

  gs_hits <- gs_summary[hit_idx, , drop = FALSE]
  fwrite(gs_hits, file.path(cfg$outdir, "GeneSet_KeywordHits.tsv"), sep = "\t")

  # 导出一个更适合人工判断的总表
  long_preview <- rbindlist(lapply(seq_along(gs), function(i) {
    data.frame(
      path = names(gs)[i],
      gene = unique(gs[[i]]),
      stringsAsFactors = FALSE
    )
  }), fill = TRUE)

  fwrite(long_preview, file.path(cfg$outdir, "GeneSet_LongPreview.tsv"), sep = "\t")
}

# --------------------------------------------------
# 4. 针对顶层 list，单独输出 names
# --------------------------------------------------
if (is.list(obj) && !is.null(names(obj))) {
  fwrite(
    data.frame(idx = seq_along(obj), name = names(obj), stringsAsFactors = FALSE),
    file.path(cfg$outdir, "TopLevel_Names.tsv"),
    sep = "\t"
  )
}

msg("Done. Output directory: ", cfg$outdir)
