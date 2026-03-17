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
  # discovery
  pe_microglia_rds = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/PsychENCODE_Science_2024/Step3_Microglia/Microglia_Clustered.rds",

  # validation
  vel_raw_input = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/Velmeshev_scRNA/rawMatrix.zip",
  vel_meta_file = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/Velmeshev_scRNA/meta.tsv",

  # enrichment result directory from your lightweight run
  enrich_dir = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Cluster2_Unbiased_Enrichment_light",

  # optional published signatures
  published_sig_dir = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/published_signatures",

  # output
  outdir = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/PublicPathway_PublishedSignature_DonorAware",

  assay = "RNA",

  # PE metadata columns
  donor_col_pe = NULL,
  dx_col_pe = NULL,
  age_col_pe = NULL,
  sex_col_pe = NULL,
  pmi_col_pe = NULL,
  rin_col_pe = NULL,
  batch_col_pe = NULL,

  # Vel metadata columns
  donor_col_vel = NULL,
  dx_col_vel = NULL,
  age_col_vel = NULL,
  sex_col_vel = NULL,
  pmi_col_vel = NULL,
  rin_col_vel = NULL,
  batch_col_vel = NULL,
  barcode_col_vel = NULL,
  celltype_col_vel = NULL,

  # selection
  top_n_go = 5,
  top_n_reactome = 5,
  top_n_hallmark = 5,

  # score parameters
  min_genes_mapped = 5,
  high_quantile = 0.75
)

dir.create(cfg$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(cfg$outdir, "plots"), recursive = TRUE, showWarnings = FALSE)

msg <- function(...) {
  cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", ..., "\n", sep = "")
}

save_session_info <- function(path) {
  writeLines(capture.output(sessionInfo()), con = path)
}

safe_write_lines <- function(x, file, empty_message = "No entries") {
  if (is.null(x) || length(x) == 0) {
    writeLines(empty_message, file)
  } else {
    writeLines(as.character(x), file)
  }
}

safe_write_table <- function(df, file) {
  if (is.null(df) || nrow(df) == 0) {
    fwrite(data.frame(message = "No result"), file, sep = "\t")
  } else {
    fwrite(as.data.frame(df), file, sep = "\t")
  }
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

normalize_sex <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  y <- x
  y[grepl("^m$|male|man", x, ignore.case = TRUE)] <- "Male"
  y[grepl("^f$|female|woman", x, ignore.case = TRUE)] <- "Female"
  y
}

safe_numeric <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

strip_bc <- function(x) {
  x <- as.character(x)
  x <- sub("-1$", "", x)
  x
}

choose_donor_col <- function(meta, user_col = NULL) {
  if (!is.null(user_col) && user_col %in% names(meta)) return(user_col)
  find_first_col(meta, c("^donor$", "donor_id", "individual", "subject", "patient", "sample", "orig.ident", "library", "case_id"))
}

choose_dx_col <- function(meta, user_col = NULL) {
  if (!is.null(user_col) && user_col %in% names(meta)) return(user_col)
  find_first_col(meta, c("^dx$", "^group$", "diagnosis", "condition", "status", "phenotype", "disease"))
}

choose_age_col <- function(meta, user_col = NULL) {
  if (!is.null(user_col) && user_col %in% names(meta)) return(user_col)
  find_first_col(meta, c("^age$", "age"))
}

choose_sex_col <- function(meta, user_col = NULL) {
  if (!is.null(user_col) && user_col %in% names(meta)) return(user_col)
  find_first_col(meta, c("^sex$", "gender"))
}

choose_pmi_col <- function(meta, user_col = NULL) {
  if (!is.null(user_col) && user_col %in% names(meta)) return(user_col)
  find_first_col(meta, c("post.?mortem.?interval", "^pmi$", "PMI"))
}

choose_rin_col <- function(meta, user_col = NULL) {
  if (!is.null(user_col) && user_col %in% names(meta)) return(user_col)
  find_first_col(meta, c("RNA.?Integrity.?Number", "^rin$", "RIN"))
}

choose_batch_col <- function(meta, user_col = NULL) {
  if (!is.null(user_col) && user_col %in% names(meta)) return(user_col)
  find_first_col(meta, c("seq.?batch", "^batch$", "capbatch", "library.?batch", "lane"))
}

detect_barcode_col <- function(meta, cell_ids, user_col = NULL) {
  if (!is.null(user_col) && user_col %in% names(meta)) return(user_col)

  cand <- names(meta)[sapply(meta, function(z) is.character(z) || is.factor(z))]
  if (length(cand) == 0) return(NULL)

  exact_hits <- sapply(cand, function(cc) sum(as.character(meta[[cc]]) %in% cell_ids))
  if (max(exact_hits) > 0) return(names(which.max(exact_hits)))

  stripped_ids <- strip_bc(cell_ids)
  stripped_hits <- sapply(cand, function(cc) sum(strip_bc(as.character(meta[[cc]])) %in% stripped_ids))
  if (max(stripped_hits) > 0) return(names(which.max(stripped_hits)))

  find_first_col(meta, c("barcode", "cell", "barcodes"))
}

detect_celltype_col <- function(meta, user_col = NULL) {
  if (!is.null(user_col) && user_col %in% names(meta)) return(user_col)

  char_cols <- names(meta)[sapply(meta, function(z) is.character(z) || is.factor(z))]
  if (length(char_cols) == 0) return(NULL)

  priority <- char_cols[grepl("cell.?type|annotation|class|subclass|type|label|cluster", char_cols, ignore.case = TRUE)]
  cand <- unique(c(priority, char_cols))

  score <- sapply(cand, function(cc) {
    v <- as.character(meta[[cc]])
    sum(grepl("microgl", v, ignore.case = TRUE), na.rm = TRUE) +
      sum(grepl("^MG$", v, ignore.case = TRUE), na.rm = TRUE)
  })

  if (max(score) == 0) return(NULL)
  names(which.max(score))
}

read_counts_auto <- function(path) {
  if (grepl("\\.zip$", path, ignore.case = TRUE)) {
    exdir <- file.path(tempdir(), paste0("Velmeshev_raw_", format(Sys.time(), "%Y%m%d_%H%M%S")))
    dir.create(exdir, recursive = TRUE, showWarnings = FALSE)
    unzip(path, exdir = exdir)
    path <- exdir
  }

  if (dir.exists(path)) {
    mtx_files <- list.files(path, pattern = "matrix\\.mtx(\\.gz)?$", recursive = TRUE, full.names = TRUE)
    if (length(mtx_files) == 0) stop("No matrix.mtx(.gz) found.")
    candidate_dirs <- unique(dirname(mtx_files))

    pick_dir <- NULL
    for (d in candidate_dirs) {
      ff <- list.files(d)
      has_mtx  <- any(grepl("matrix\\.mtx(\\.gz)?$", ff))
      has_feat <- any(grepl("(features|genes)\\.tsv(\\.gz)?$", ff))
      has_bc   <- any(grepl("barcodes\\.tsv(\\.gz)?$", ff))
      if (has_mtx && has_feat && has_bc) {
        pick_dir <- d
        break
      }
    }
    if (is.null(pick_dir)) pick_dir <- candidate_dirs[1]

    msg("Read10X directory: ", pick_dir)
    x <- Read10X(data.dir = pick_dir)
    if (is.list(x)) {
      if ("Gene Expression" %in% names(x)) {
        x <- x[["Gene Expression"]]
      } else {
        x <- x[[1]]
      }
    }
    return(x)
  }

  stop("Unsupported raw_input format: ", path)
}

align_meta_to_cells <- function(meta, barcode_col, cell_ids) {
  bc <- as.character(meta[[barcode_col]])
  idx <- match(cell_ids, bc)

  match_rate <- mean(!is.na(idx))
  if (match_rate < 0.5) {
    idx <- match(strip_bc(cell_ids), strip_bc(bc))
  }

  meta2 <- meta[idx, , drop = FALSE]
  rownames(meta2) <- cell_ids
  meta2
}

harmonize_genes <- function(genes, features) {
  genes <- unique(as.character(genes))
  genes <- genes[!is.na(genes) & nzchar(genes)]

  exact <- intersect(genes, features)
  if (length(exact) >= 5) return(exact)

  fmap <- setNames(features, toupper(features))
  mapped <- unname(fmap[toupper(genes)])
  mapped <- unique(mapped[!is.na(mapped)])
  mapped
}

calc_wilcox_p <- function(df, value_col, group_col = "dx") {
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

mean_diff_asd_minus_ctrl <- function(df, value_col, group_col = "dx") {
  df <- as.data.frame(df)
  grp <- as.character(df[[group_col]])
  val <- df[[value_col]]
  if (!all(c("ASD", "Control") %in% unique(grp))) return(NA_real_)
  mean(val[grp == "ASD"], na.rm = TRUE) - mean(val[grp == "Control"], na.rm = TRUE)
}

extract_coef <- function(fit, term = "dxASD") {
  if (is.null(fit)) {
    return(data.frame(
      term = term, beta = NA_real_, se = NA_real_, stat = NA_real_, p = NA_real_,
      effect = NA_real_, conf.low = NA_real_, conf.high = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  sm <- summary(fit)
  cf <- coef(sm)
  rn <- rownames(cf)
  if (!(term %in% rn)) {
    alt <- rn[grepl("^dx", rn)]
    if (length(alt) == 0) {
      return(data.frame(
        term = term, beta = NA_real_, se = NA_real_, stat = NA_real_, p = NA_real_,
        effect = NA_real_, conf.low = NA_real_, conf.high = NA_real_,
        stringsAsFactors = FALSE
      ))
    }
    term <- alt[1]
  }

  beta <- cf[term, 1]
  se   <- cf[term, 2]
  stat <- cf[term, 3]
  pval <- cf[term, 4]

  data.frame(
    term = term,
    beta = beta,
    se = se,
    stat = stat,
    p = pval,
    effect = beta,
    conf.low = beta - 1.96 * se,
    conf.high = beta + 1.96 * se,
    stringsAsFactors = FALSE
  )
}

fit_glm_safe <- function(formula_obj, dat, family_obj) {
  tryCatch(
    glm(formula_obj, data = dat, family = family_obj),
    error = function(e) NULL,
    warning = function(w) {
      invokeRestart("muffleWarning")
    }
  )
}

fit_lm_safe <- function(formula_obj, dat) {
  tryCatch(
    lm(formula_obj, data = dat),
    error = function(e) NULL,
    warning = function(w) {
      invokeRestart("muffleWarning")
    }
  )
}

safe_module_score <- function(data_mat, genes_use) {
  if (length(genes_use) == 0) return(rep(NA_real_, ncol(data_mat)))

  sub <- as.matrix(data_mat[genes_use, , drop = FALSE])

  z <- t(apply(sub, 1, function(v) {
    s <- sd(v, na.rm = TRUE)
    if (is.na(s) || s == 0) {
      rep(0, length(v))
    } else {
      (v - mean(v, na.rm = TRUE)) / s
    }
  }))

  colMeans(z, na.rm = TRUE)
}

determine_usable_covars <- function(donor_meta) {
  usable_covars <- c()

  if ("age" %in% colnames(donor_meta)) {
    x <- donor_meta$age
    if (sum(!is.na(x)) >= 6 && length(unique(x[!is.na(x)])) >= 3) usable_covars <- c(usable_covars, "age")
  }
  if ("sex" %in% colnames(donor_meta)) {
    x <- donor_meta$sex
    if (sum(!is.na(x)) >= 6 && length(unique(x[!is.na(x)])) >= 2) usable_covars <- c(usable_covars, "sex")
  }
  if ("pmi" %in% colnames(donor_meta)) {
    x <- donor_meta$pmi
    if (sum(!is.na(x)) >= 6 && length(unique(x[!is.na(x)])) >= 3) usable_covars <- c(usable_covars, "pmi")
  }
  if ("rin" %in% colnames(donor_meta)) {
    x <- donor_meta$rin
    if (sum(!is.na(x)) >= 6 && length(unique(x[!is.na(x)])) >= 3) usable_covars <- c(usable_covars, "rin")
  }
  if ("batch" %in% colnames(donor_meta)) {
    x <- donor_meta$batch
    if (sum(!is.na(x)) >= 6 && length(unique(x[!is.na(x)])) >= 2 && length(unique(x[!is.na(x)])) <= 6) {
      usable_covars <- c(usable_covars, "batch")
    }
  }

  usable_covars
}

read_top_terms <- function(file, top_n = 5, positive_only = FALSE) {
  if (!file.exists(file)) return(character(0))
  x <- fread(file, data.table = FALSE)

  term_col <- NULL
  for (cc in c("term_name", "Description", "ID", "term")) {
    if (cc %in% colnames(x)) {
      term_col <- cc
      break
    }
  }
  if (is.null(term_col)) return(character(0))

  p_col <- NULL
  for (cc in c("p_adj", "p.adjust", "p_value", "pvalue")) {
    if (cc %in% colnames(x)) {
      p_col <- cc
      break
    }
  }
  if (is.null(p_col)) return(character(0))

  if (positive_only) {
    if ("NES_like" %in% colnames(x)) {
      x <- x[x$NES_like > 0, , drop = FALSE]
    } else if ("direction" %in% colnames(x)) {
      x <- x[grepl("top|up|positive", x$direction, ignore.case = TRUE), , drop = FALSE]
    }
  }

  x <- x[order(x[[p_col]], decreasing = FALSE), , drop = FALSE]
  x <- head(x, top_n)
  unique(as.character(x[[term_col]]))
}

get_msigdb_term2gene <- function() {
  x <- as.data.frame(msigdbr::msigdbr(species = "Homo sapiens"))

  term_col <- find_first_col(x, c("^gs_name$", "^term$", "geneset"))
  gene_col <- find_first_col(x, c("^gene_symbol$", "human_gene_symbol", "symbol"))
  coll_col <- find_first_col(x, c("^gs_collection$", "^collection$", "^gs_cat$", "^category$"))
  subcoll_col <- find_first_col(x, c("^gs_subcollection$", "^subcollection$", "^gs_subcat$", "^subcategory$"))
  desc_col <- find_first_col(x, c("^gs_description$", "description"))

  if (is.null(term_col) || is.null(gene_col) || is.null(coll_col)) {
    stop("Could not identify key columns in msigdbr output.")
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

extract_msig_terms <- function(msig, terms_keep = character(0), type = c("hallmark","reactome","go_bp")) {
  type <- match.arg(type)

  df <- msig$df
  term_col <- msig$term_col
  gene_col <- msig$gene_col
  coll_col <- msig$coll_col
  subcoll_col <- msig$subcoll_col
  desc_col <- msig$desc_col

  coll <- toupper(df[[coll_col]])
  subc <- if (!is.null(subcoll_col)) toupper(df[[subcoll_col]]) else rep(NA_character_, nrow(df))

  keep_base <- rep(FALSE, nrow(df))
  if (type == "hallmark") keep_base <- coll == "H"
  if (type == "reactome") keep_base <- coll == "C2" & grepl("REACTOME", subc)
  if (type == "go_bp") keep_base <- coll == "C5" & grepl("GO:BP|GOBP|BP", subc)

  out <- df[keep_base, , drop = FALSE]
  if (nrow(out) == 0) return(NULL)

  if (length(terms_keep) > 0) {
    out <- out[out[[term_col]] %in% terms_keep, , drop = FALSE]
  }

  if (nrow(out) == 0) return(NULL)

  out2 <- data.frame(
    term_name = out[[term_col]],
    gene_symbol = out[[gene_col]],
    stringsAsFactors = FALSE
  )
  out2 <- out2[!is.na(out2$gene_symbol) & out2$gene_symbol != "", , drop = FALSE]

  if (!is.null(desc_col) && desc_col %in% colnames(out)) {
    out2$term_description <- out[[desc_col]]
  } else {
    out2$term_description <- out2$term_name
  }

  unique(out2)
}

load_published_signatures <- function(sig_dir) {
  if (!dir.exists(sig_dir)) return(list())

  files <- list.files(sig_dir, pattern = "\\.(txt|tsv|csv)$", full.names = TRUE)
  if (length(files) == 0) return(list())

  out <- list()

  for (f in files) {
    nm <- sub("\\.[^.]+$", "", basename(f))

    if (grepl("\\.txt$", f, ignore.case = TRUE)) {
      genes <- readLines(f, warn = FALSE)
      genes <- trimws(genes)
      genes <- genes[genes != ""]
    } else {
      x <- fread(f, data.table = FALSE)
      gene_col <- find_first_col(x, c("^gene$", "gene_symbol", "symbol", "genes"))
      if (is.null(gene_col)) next
      genes <- as.character(x[[gene_col]])
    }

    genes <- unique(genes[!is.na(genes) & genes != ""])
    out[[nm]] <- genes
  }

  out
}

analyze_cohort <- function(seu, cohort_name,
                           donor_col, dx_col, age_col = NULL, sex_col = NULL,
                           pmi_col = NULL, rin_col = NULL, batch_col = NULL,
                           gene_sets, outdir, min_genes_mapped = 5, high_quantile = 0.75) {

  msg("Analyzing cohort: ", cohort_name)

  md <- seu@meta.data
  md$cell_id <- colnames(seu)
  md$donor <- as.character(md[[donor_col]])
  md$dx <- normalize_dx(md[[dx_col]])

  if (!is.null(age_col) && age_col %in% colnames(md))   md$age   <- safe_numeric(md[[age_col]])
  if (!is.null(sex_col) && sex_col %in% colnames(md))   md$sex   <- normalize_sex(md[[sex_col]])
  if (!is.null(pmi_col) && pmi_col %in% colnames(md))   md$pmi   <- safe_numeric(md[[pmi_col]])
  if (!is.null(rin_col) && rin_col %in% colnames(md))   md$rin   <- safe_numeric(md[[rin_col]])
  if (!is.null(batch_col) && batch_col %in% colnames(md)) md$batch <- as.character(md[[batch_col]])

  donor_meta <- md %>%
    group_by(donor) %>%
    summarise(
      dx = dplyr::first(dx),
      age = if ("age" %in% colnames(md)) median(age, na.rm = TRUE) else NA_real_,
      sex = if ("sex" %in% colnames(md)) dplyr::first(na.omit(sex)) else NA_character_,
      pmi = if ("pmi" %in% colnames(md)) median(pmi, na.rm = TRUE) else NA_real_,
      rin = if ("rin" %in% colnames(md)) median(rin, na.rm = TRUE) else NA_real_,
      batch = if ("batch" %in% colnames(md)) dplyr::first(na.omit(batch)) else NA_character_,
      .groups = "drop"
    ) %>%
    as.data.frame()

  usable_covars <- determine_usable_covars(donor_meta)
  safe_write_lines(usable_covars,
                   file.path(outdir, paste0(cohort_name, "_UsableCovariates.txt")),
                   empty_message = "No usable covariates")

  data_mat <- GetAssayData(seu, assay = DefaultAssay(seu), slot = "data")

  audit_tbl <- data.frame(
    item = c("cohort", "n_cells", "n_features", "donor_col", "dx_col", "age_col", "sex_col", "pmi_col", "rin_col", "batch_col"),
    value = c(cohort_name, ncol(seu), nrow(seu), donor_col, dx_col,
              ifelse(is.null(age_col), NA, age_col),
              ifelse(is.null(sex_col), NA, sex_col),
              ifelse(is.null(pmi_col), NA, pmi_col),
              ifelse(is.null(rin_col), NA, rin_col),
              ifelse(is.null(batch_col), NA, batch_col))
  )
  fwrite(audit_tbl, file.path(outdir, paste0(cohort_name, "_Audit.tsv")), sep = "\t")

  gs_audit <- list()
  donor_rows <- list()
  stat_rows <- list()

  for (gs_name in names(gene_sets)) {
    genes_in <- gene_sets[[gs_name]]
    genes_use <- harmonize_genes(genes_in, rownames(data_mat))

    gs_audit[[length(gs_audit) + 1]] <- data.frame(
      signature = gs_name,
      n_input = length(unique(genes_in)),
      n_mapped = length(unique(genes_use)),
      mapped_preview = paste(head(unique(genes_use), 15), collapse = ", "),
      stringsAsFactors = FALSE
    )

    if (length(genes_use) < min_genes_mapped) next

    score <- safe_module_score(data_mat, genes_use)

    cell_df <- data.frame(
      donor = md$donor,
      dx = md$dx,
      score = score,
      stringsAsFactors = FALSE
    )

    thr <- as.numeric(quantile(cell_df$score, probs = high_quantile, na.rm = TRUE))
    cell_df$is_high <- as.integer(cell_df$score >= thr)

    donor_df <- cell_df %>%
      group_by(donor, dx) %>%
      summarise(
        n_cells = n(),
        n_high = sum(is_high, na.rm = TRUE),
        high_frac = mean(is_high, na.rm = TRUE),
        score_mean = mean(score, na.rm = TRUE),
        score_median = median(score, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      left_join(donor_meta, by = c("donor", "dx")) %>%
      as.data.frame()

    donor_df$signature <- gs_name
    donor_df$threshold_high <- thr
    donor_rows[[length(donor_rows) + 1]] <- donor_df

    donor_df$dx <- factor(donor_df$dx, levels = c("Control", "ASD"))
    if ("sex" %in% colnames(donor_df)) donor_df$sex <- factor(donor_df$sex)
    if ("batch" %in% colnames(donor_df)) donor_df$batch <- factor(donor_df$batch)

    p_w_mean <- calc_wilcox_p(donor_df, "score_mean", "dx")
    eff_mean <- mean_diff_asd_minus_ctrl(donor_df, "score_mean", "dx")

    p_w_high <- calc_wilcox_p(donor_df, "high_frac", "dx")
    eff_high <- mean_diff_asd_minus_ctrl(donor_df, "high_frac", "dx")

    lm_basic <- fit_lm_safe(score_mean ~ dx, donor_df)
    lm_basic_res <- extract_coef(lm_basic, "dxASD")

    lm_adj_res <- data.frame(
      term = "dxASD", beta = NA_real_, se = NA_real_, stat = NA_real_, p = NA_real_,
      effect = NA_real_, conf.low = NA_real_, conf.high = NA_real_, stringsAsFactors = FALSE
    )
    n_complete_lm <- NA_integer_
    if (length(usable_covars) > 0) {
      need_cols <- unique(c("score_mean", "dx", usable_covars))
      dat_lm <- donor_df[, need_cols, drop = FALSE]
      dat_lm <- dat_lm[complete.cases(dat_lm), , drop = FALSE]
      n_complete_lm <- nrow(dat_lm)
      if (nrow(dat_lm) >= 8 && length(unique(dat_lm$dx)) == 2) {
        f_lm_adj <- as.formula(paste("score_mean ~ dx +", paste(usable_covars, collapse = " + ")))
        lm_adj <- fit_lm_safe(f_lm_adj, dat_lm)
        lm_adj_res <- extract_coef(lm_adj, "dxASD")
      }
    }

    qb_basic <- fit_glm_safe(cbind(n_high, n_cells - n_high) ~ dx, donor_df, quasibinomial())
    qb_basic_res <- extract_coef(qb_basic, "dxASD")
    qb_basic_or <- exp(qb_basic_res$beta)

    qb_adj_res <- data.frame(
      term = "dxASD", beta = NA_real_, se = NA_real_, stat = NA_real_, p = NA_real_,
      effect = NA_real_, conf.low = NA_real_, conf.high = NA_real_, stringsAsFactors = FALSE
    )
    qb_adj_or <- NA_real_
    n_complete_qb <- NA_integer_
    if (length(usable_covars) > 0) {
      need_cols <- unique(c("n_high", "n_cells", "dx", usable_covars))
      dat_qb <- donor_df[, need_cols, drop = FALSE]
      dat_qb <- dat_qb[complete.cases(dat_qb), , drop = FALSE]
      n_complete_qb <- nrow(dat_qb)
      if (nrow(dat_qb) >= 8 && length(unique(dat_qb$dx)) == 2) {
        f_qb_adj <- as.formula(paste("cbind(n_high, n_cells - n_high) ~ dx +", paste(usable_covars, collapse = " + ")))
        qb_adj <- fit_glm_safe(f_qb_adj, dat_qb, quasibinomial())
        qb_adj_res <- extract_coef(qb_adj, "dxASD")
        qb_adj_or <- exp(qb_adj_res$beta)
      }
    }

    stat_rows[[length(stat_rows) + 1]] <- data.frame(
      cohort = cohort_name,
      signature = gs_name,
      n_donors = length(unique(donor_df$donor)),
      n_ASD = length(unique(donor_df$donor[donor_df$dx == "ASD"])),
      n_Control = length(unique(donor_df$donor[donor_df$dx == "Control"])),

      meanScore_ASD = mean(donor_df$score_mean[donor_df$dx == "ASD"], na.rm = TRUE),
      meanScore_Control = mean(donor_df$score_mean[donor_df$dx == "Control"], na.rm = TRUE),
      effect_meanScore_ASD_minus_Control = eff_mean,
      wilcox_p_meanScore = p_w_mean,

      meanHighFrac_ASD = mean(donor_df$high_frac[donor_df$dx == "ASD"], na.rm = TRUE),
      meanHighFrac_Control = mean(donor_df$high_frac[donor_df$dx == "Control"], na.rm = TRUE),
      effect_highFrac_ASD_minus_Control = eff_high,
      wilcox_p_highFrac = p_w_high,

      lm_basic_beta = lm_basic_res$beta,
      lm_basic_p = lm_basic_res$p,
      lm_adj_beta = lm_adj_res$beta,
      lm_adj_p = lm_adj_res$p,
      lm_adj_n_complete = n_complete_lm,

      qb_basic_beta = qb_basic_res$beta,
      qb_basic_OR = qb_basic_or,
      qb_basic_p = qb_basic_res$p,
      qb_adj_beta = qb_adj_res$beta,
      qb_adj_OR = qb_adj_or,
      qb_adj_p = qb_adj_res$p,
      qb_adj_n_complete = n_complete_qb,

      threshold_high = thr,
      covariates = ifelse(length(usable_covars) == 0, "", paste(usable_covars, collapse = "+")),
      stringsAsFactors = FALSE
    )
  }

  gs_audit_df <- bind_rows(gs_audit)
  donor_scores_df <- bind_rows(donor_rows)
  stats_df <- bind_rows(stat_rows)

  if (nrow(stats_df) > 0) {
    stats_df$wilcox_fdr_meanScore <- p.adjust(stats_df$wilcox_p_meanScore, method = "BH")
    stats_df$wilcox_fdr_highFrac  <- p.adjust(stats_df$wilcox_p_highFrac, method = "BH")
    stats_df$qb_basic_fdr <- p.adjust(stats_df$qb_basic_p, method = "BH")
    stats_df$qb_adj_fdr   <- p.adjust(stats_df$qb_adj_p, method = "BH")
  }

  fwrite(gs_audit_df, file.path(outdir, paste0(cohort_name, "_GeneSetAudit.tsv")), sep = "\t")
  fwrite(donor_scores_df, file.path(outdir, paste0(cohort_name, "_DonorScores.tsv")), sep = "\t")
  fwrite(stats_df, file.path(outdir, paste0(cohort_name, "_Stats.tsv")), sep = "\t")

  if (nrow(stats_df) > 0) {
    p_rank_mean <- stats_df %>%
      mutate(neglog10p = -log10(pmax(wilcox_p_meanScore, 1e-300))) %>%
      arrange(wilcox_p_meanScore, desc(effect_meanScore_ASD_minus_Control))
    p_rank_mean$signature <- factor(p_rank_mean$signature, levels = p_rank_mean$signature)

    g1 <- ggplot(p_rank_mean, aes(x = effect_meanScore_ASD_minus_Control, y = signature)) +
      geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
      geom_point(aes(size = neglog10p)) +
      theme_classic(base_size = 12) +
      labs(
        title = paste0(cohort_name, ": donor-level mean score ranking"),
        x = "Mean score difference (ASD - Control)",
        y = "Signature / pathway",
        size = "-log10(p)"
      )
    ggsave(file.path(outdir, "plots", paste0(cohort_name, "_Ranking_MeanScore.pdf")), g1, width = 8.2, height = 6)

    p_rank_high <- stats_df %>%
      mutate(neglog10p = -log10(pmax(qb_basic_p, 1e-300))) %>%
      arrange(qb_basic_p, desc(effect_highFrac_ASD_minus_Control))
    p_rank_high$signature <- factor(p_rank_high$signature, levels = p_rank_high$signature)

    g2 <- ggplot(p_rank_high, aes(x = effect_highFrac_ASD_minus_Control, y = signature)) +
      geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
      geom_point(aes(size = neglog10p)) +
      theme_classic(base_size = 12) +
      labs(
        title = paste0(cohort_name, ": donor-level high-fraction ranking"),
        x = "High-fraction difference (ASD - Control)",
        y = "Signature / pathway",
        size = "-log10(p)"
      )
    ggsave(file.path(outdir, "plots", paste0(cohort_name, "_Ranking_HighFrac.pdf")), g2, width = 8.2, height = 6)
  }

  return(list(
    stats = stats_df,
    donor_scores = donor_scores_df,
    audit = gs_audit_df
  ))
}

# --------------------------------------------------
# 1. load PE object
# --------------------------------------------------
msg("Reading PsychENCODE microglia object ...")
pe <- readRDS(cfg$pe_microglia_rds)
if (!inherits(pe, "Seurat")) stop("PE RDS is not a Seurat object.")
if (cfg$assay %in% Assays(pe)) DefaultAssay(pe) <- cfg$assay
if (ncol(GetAssayData(pe, assay = DefaultAssay(pe), slot = "data")) == 0) {
  pe <- NormalizeData(pe, verbose = FALSE)
}

pe_meta <- pe@meta.data
pe_donor_col <- choose_donor_col(pe_meta, cfg$donor_col_pe)
pe_dx_col    <- choose_dx_col(pe_meta, cfg$dx_col_pe)
pe_age_col   <- choose_age_col(pe_meta, cfg$age_col_pe)
pe_sex_col   <- choose_sex_col(pe_meta, cfg$sex_col_pe)
pe_pmi_col   <- choose_pmi_col(pe_meta, cfg$pmi_col_pe)
pe_rin_col   <- choose_rin_col(pe_meta, cfg$rin_col_pe)
pe_batch_col <- choose_batch_col(pe_meta, cfg$batch_col_pe)

if (is.null(pe_donor_col)) stop("Could not identify PE donor column.")
if (is.null(pe_dx_col)) stop("Could not identify PE diagnosis column.")

# --------------------------------------------------
# 2. load Vel
# --------------------------------------------------
msg("Reading Velmeshev counts and metadata ...")
vel_counts <- read_counts_auto(cfg$vel_raw_input)
vel_meta0 <- fread(cfg$vel_meta_file, data.table = FALSE)
vel_meta0 <- as.data.frame(vel_meta0)

vel_barcode_col <- detect_barcode_col(vel_meta0, colnames(vel_counts), cfg$barcode_col_vel)
if (is.null(vel_barcode_col)) stop("Could not identify Vel barcode column.")

vel_meta1 <- align_meta_to_cells(vel_meta0, vel_barcode_col, colnames(vel_counts))
keep <- !is.na(vel_meta1[[vel_barcode_col]])
vel_meta1 <- vel_meta1[keep, , drop = FALSE]
vel_counts <- vel_counts[, rownames(vel_meta1), drop = FALSE]

vel_celltype_col <- detect_celltype_col(vel_meta1, cfg$celltype_col_vel)
if (is.null(vel_celltype_col)) stop("Could not identify Vel cell-type column containing microglia labels.")

ctv <- as.character(vel_meta1[[vel_celltype_col]])
mg_idx <- grepl("microgl", ctv, ignore.case = TRUE) | grepl("^MG$", ctv, ignore.case = TRUE)
if (sum(mg_idx) == 0) stop("No microglia cells matched in Vel metadata.")

vel_meta <- vel_meta1[mg_idx, , drop = FALSE]
vel_counts_mg <- vel_counts[, rownames(vel_meta), drop = FALSE]

vel <- CreateSeuratObject(
  counts = vel_counts_mg,
  meta.data = vel_meta,
  assay = cfg$assay,
  project = "VelmeshevMG",
  min.cells = 0,
  min.features = 0
)
DefaultAssay(vel) <- cfg$assay
vel <- NormalizeData(vel, verbose = FALSE)

vel_donor_col <- choose_donor_col(vel@meta.data, cfg$donor_col_vel)
vel_dx_col    <- choose_dx_col(vel@meta.data, cfg$dx_col_vel)
vel_age_col   <- choose_age_col(vel@meta.data, cfg$age_col_vel)
vel_sex_col   <- choose_sex_col(vel@meta.data, cfg$sex_col_vel)
vel_pmi_col   <- choose_pmi_col(vel@meta.data, cfg$pmi_col_vel)
vel_rin_col   <- choose_rin_col(vel@meta.data, cfg$rin_col_vel)
vel_batch_col <- choose_batch_col(vel@meta.data, cfg$batch_col_vel)

if (is.null(vel_donor_col)) stop("Could not identify Vel donor column.")
if (is.null(vel_dx_col)) stop("Could not identify Vel diagnosis column.")

# --------------------------------------------------
# 3. build public pathway gene sets from enrichment outputs
# --------------------------------------------------
msg("Loading public pathways from enrichment outputs + msigdbr ...")

top_go_terms <- read_top_terms(
  file.path(cfg$enrich_dir, "Top15_GO_BP_Up_localMSigDB.tsv"),
  top_n = cfg$top_n_go,
  positive_only = FALSE
)

top_reactome_terms <- read_top_terms(
  file.path(cfg$enrich_dir, "Top15_Reactome_Up_localMSigDB.tsv"),
  top_n = cfg$top_n_reactome,
  positive_only = FALSE
)

top_hallmark_terms <- read_top_terms(
  file.path(cfg$enrich_dir, "Top15_Hallmark_RankEnrichment.tsv"),
  top_n = cfg$top_n_hallmark,
  positive_only = TRUE
)

msig <- get_msigdb_term2gene()

go_t2g <- extract_msig_terms(msig, terms_keep = top_go_terms, type = "go_bp")
react_t2g <- extract_msig_terms(msig, terms_keep = top_reactome_terms, type = "reactome")
hall_t2g <- extract_msig_terms(msig, terms_keep = top_hallmark_terms, type = "hallmark")

safe_write_table(go_t2g, file.path(cfg$outdir, "Selected_GO_term2gene.tsv"))
safe_write_table(react_t2g, file.path(cfg$outdir, "Selected_Reactome_term2gene.tsv"))
safe_write_table(hall_t2g, file.path(cfg$outdir, "Selected_Hallmark_term2gene.tsv"))

public_gene_sets <- list()

if (!is.null(go_t2g) && nrow(go_t2g) > 0) {
  tmp <- split(go_t2g$gene_symbol, go_t2g$term_name)
  for (nm in names(tmp)) public_gene_sets[[paste0("GOBP__", nm)]] <- unique(tmp[[nm]])
}
if (!is.null(react_t2g) && nrow(react_t2g) > 0) {
  tmp <- split(react_t2g$gene_symbol, react_t2g$term_name)
  for (nm in names(tmp)) public_gene_sets[[paste0("REACTOME__", nm)]] <- unique(tmp[[nm]])
}
if (!is.null(hall_t2g) && nrow(hall_t2g) > 0) {
  tmp <- split(hall_t2g$gene_symbol, hall_t2g$term_name)
  for (nm in names(tmp)) public_gene_sets[[paste0("HALLMARK__", nm)]] <- unique(tmp[[nm]])
}

safe_write_lines(names(public_gene_sets),
                 file.path(cfg$outdir, "Selected_PublicPathways.txt"),
                 empty_message = "No public pathways selected")

# --------------------------------------------------
# 4. load published signatures
# --------------------------------------------------
msg("Loading published microglia signatures ...")
published_signatures <- load_published_signatures(cfg$published_sig_dir)

if (length(published_signatures) > 0) {
  published_signatures <- setNames(
    published_signatures,
    paste0("PUBSIG__", names(published_signatures))
  )
  safe_write_lines(names(published_signatures),
                   file.path(cfg$outdir, "Loaded_PublishedSignatures.txt"),
                   empty_message = "No published signatures loaded")
} else {
  published_signatures <- list()
  safe_write_lines(NULL,
                   file.path(cfg$outdir, "Loaded_PublishedSignatures.txt"),
                   empty_message = "No published signatures loaded")
}

# --------------------------------------------------
# 5. combined gene sets
# --------------------------------------------------
gene_sets <- c(public_gene_sets, published_signatures)
if (length(gene_sets) == 0) {
  stop("No public pathways or published signatures loaded.")
}

# --------------------------------------------------
# 6. analyze both cohorts
# --------------------------------------------------
pe_res <- analyze_cohort(
  seu = pe,
  cohort_name = "PsychENCODE",
  donor_col = pe_donor_col,
  dx_col = pe_dx_col,
  age_col = pe_age_col,
  sex_col = pe_sex_col,
  pmi_col = pe_pmi_col,
  rin_col = pe_rin_col,
  batch_col = pe_batch_col,
  gene_sets = gene_sets,
  outdir = cfg$outdir,
  min_genes_mapped = cfg$min_genes_mapped,
  high_quantile = cfg$high_quantile
)

vel_res <- analyze_cohort(
  seu = vel,
  cohort_name = "Velmeshev",
  donor_col = vel_donor_col,
  dx_col = vel_dx_col,
  age_col = vel_age_col,
  sex_col = vel_sex_col,
  pmi_col = vel_pmi_col,
  rin_col = vel_rin_col,
  batch_col = vel_batch_col,
  gene_sets = gene_sets,
  outdir = cfg$outdir,
  min_genes_mapped = cfg$min_genes_mapped,
  high_quantile = cfg$high_quantile
)

# --------------------------------------------------
# 7. cross-cohort summary
# --------------------------------------------------
if (nrow(pe_res$stats) > 0 && nrow(vel_res$stats) > 0) {
  cross <- pe_res$stats %>%
    select(signature, effect_meanScore_ASD_minus_Control, wilcox_p_meanScore, effect_highFrac_ASD_minus_Control, qb_basic_p) %>%
    rename_with(~ paste0("PE_", .x), -signature) %>%
    full_join(
      vel_res$stats %>%
        select(signature, effect_meanScore_ASD_minus_Control, wilcox_p_meanScore, effect_highFrac_ASD_minus_Control, qb_basic_p) %>%
        rename_with(~ paste0("Vel_", .x), -signature),
      by = "signature"
    ) %>%
    mutate(
      same_direction_meanScore = sign(PE_effect_meanScore_ASD_minus_Control) == sign(Vel_effect_meanScore_ASD_minus_Control),
      same_direction_highFrac = sign(PE_effect_highFrac_ASD_minus_Control) == sign(Vel_effect_highFrac_ASD_minus_Control)
    )

  fwrite(cross, file.path(cfg$outdir, "CrossCohort_Summary.tsv"), sep = "\t")
} else {
  fwrite(data.frame(message = "One or both cohorts produced no stats"),
         file.path(cfg$outdir, "CrossCohort_Summary.tsv"),
         sep = "\t")
}

save_session_info(file.path(cfg$outdir, "sessionInfo_PublicPathway_PublishedSignature_DonorAware.txt"))
msg("Done.")
