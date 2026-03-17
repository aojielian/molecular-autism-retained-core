#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(dplyr)
  library(ggplot2)
})

options(stringsAsFactors = FALSE)

cfg <- list(
  microglia_rds = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/PsychENCODE_Science_2024/Step3_Microglia/Microglia_Clustered.rds",
  outdir        = "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Figure1_DonorAbundance_GLM",
  assay         = "RNA",
  cluster_col   = NULL,
  donor_col     = NULL,
  dx_col        = NULL,

  # 优先尝试自动识别这些协变量；识别不到就跳过
  age_col       = NULL,
  sex_col       = NULL,
  pmi_col       = NULL,
  rin_col       = NULL,
  batch_col     = NULL,

  target_cluster = "2"
)

dir.create(cfg$outdir, recursive = TRUE, showWarnings = FALSE)

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

extract_coef <- function(fit, term = "dxASD") {
  sm <- summary(fit)
  cf <- coef(sm)
  rn <- rownames(cf)
  if (!(term %in% rn)) {
    alt <- rn[grepl("^dx", rn)]
    if (length(alt) == 0) {
      return(data.frame(
        term = term, beta = NA_real_, se = NA_real_, stat = NA_real_, p = NA_real_,
        OR = NA_real_, conf.low = NA_real_, conf.high = NA_real_,
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
    OR = exp(beta),
    conf.low = exp(beta - 1.96 * se),
    conf.high = exp(beta + 1.96 * se),
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

# --------------------------------------------------
# 1. 读取对象
# --------------------------------------------------
msg("Reading object ...")
obj <- readRDS(cfg$microglia_rds)
if (!inherits(obj, "Seurat")) stop("Input RDS is not a Seurat object.")

if (cfg$assay %in% Assays(obj)) {
  DefaultAssay(obj) <- cfg$assay
}

meta <- obj@meta.data

cluster_col <- choose_cluster_col(obj, cfg$cluster_col)
obj$cluster_use <- as.character(meta[[cluster_col]])
Idents(obj) <- obj$cluster_use

donor_col <- choose_donor_col(meta, cfg$donor_col)
dx_col    <- choose_dx_col(meta, cfg$dx_col)
age_col   <- choose_age_col(meta, cfg$age_col)
sex_col   <- choose_sex_col(meta, cfg$sex_col)
pmi_col   <- choose_pmi_col(meta, cfg$pmi_col)
rin_col   <- choose_rin_col(meta, cfg$rin_col)
batch_col <- choose_batch_col(meta, cfg$batch_col)

if (is.null(donor_col)) stop("Could not identify donor column.")
if (is.null(dx_col)) stop("Could not identify diagnosis column.")

audit <- data.frame(
  item = c("cluster_col", "donor_col", "dx_col", "age_col", "sex_col", "pmi_col", "rin_col", "batch_col",
           "n_cells", "n_features"),
  value = c(cluster_col, donor_col, dx_col,
            ifelse(is.null(age_col), NA, age_col),
            ifelse(is.null(sex_col), NA, sex_col),
            ifelse(is.null(pmi_col), NA, pmi_col),
            ifelse(is.null(rin_col), NA, rin_col),
            ifelse(is.null(batch_col), NA, batch_col),
            ncol(obj), nrow(obj))
)
fwrite(audit, file.path(cfg$outdir, "Audit.tsv"), sep = "\t")

# --------------------------------------------------
# 2. donor-level 聚合
# --------------------------------------------------
msg("Aggregating donor-level counts ...")
meta2 <- obj@meta.data
meta2$cell_id  <- colnames(obj)
meta2$cluster  <- as.character(meta2[[cluster_col]])
meta2$donor    <- as.character(meta2[[donor_col]])
meta2$dx       <- normalize_dx(meta2[[dx_col]])

if (!is.null(age_col))   meta2$age   <- safe_numeric(meta2[[age_col]])
if (!is.null(sex_col))   meta2$sex   <- normalize_sex(meta2[[sex_col]])
if (!is.null(pmi_col))   meta2$pmi   <- safe_numeric(meta2[[pmi_col]])
if (!is.null(rin_col))   meta2$rin   <- safe_numeric(meta2[[rin_col]])
if (!is.null(batch_col)) meta2$batch <- as.character(meta2[[batch_col]])

cluster_levels <- unique(meta2$cluster)
num_try <- suppressWarnings(as.numeric(cluster_levels))
if (all(!is.na(num_try))) {
  cluster_levels <- as.character(sort(num_try))
} else {
  cluster_levels <- sort(cluster_levels)
}

# donor-level metadata（每个 donor 一行）
donor_meta <- meta2 %>%
  group_by(donor) %>%
  summarise(
    dx = dplyr::first(dx),
    age = if ("age" %in% colnames(cur_data())) median(age, na.rm = TRUE) else NA_real_,
    sex = if ("sex" %in% colnames(cur_data())) dplyr::first(na.omit(sex)) else NA_character_,
    pmi = if ("pmi" %in% colnames(cur_data())) median(pmi, na.rm = TRUE) else NA_real_,
    rin = if ("rin" %in% colnames(cur_data())) median(rin, na.rm = TRUE) else NA_real_,
    batch = if ("batch" %in% colnames(cur_data())) dplyr::first(na.omit(batch)) else NA_character_,
    n_total = n(),
    .groups = "drop"
  )

counts_df <- meta2 %>%
  count(donor, dx, cluster, name = "n_cluster") %>%
  left_join(donor_meta %>% select(donor, n_total), by = "donor") %>%
  mutate(
    n_other = n_total - n_cluster,
    prop = n_cluster / n_total
  ) %>%
  arrange(cluster, dx, donor)

fwrite(counts_df, file.path(cfg$outdir, "DonorCluster_Counts.tsv"), sep = "\t")
fwrite(donor_meta, file.path(cfg$outdir, "Donor_Metadata.tsv"), sep = "\t")

# --------------------------------------------------
# 3. 构建协变量
# --------------------------------------------------
msg("Selecting usable covariates ...")
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
  if (sum(!is.na(x)) >= 6 && length(unique(x[!is.na(x)])) >= 2 && length(unique(x[!is.na(x)])) <= 5) {
    usable_covars <- c(usable_covars, "batch")
  }
}

writeLines(usable_covars, file.path(cfg$outdir, "Usable_Covariates.txt"))

# --------------------------------------------------
# 4. 每个 cluster 拟合 binomial / quasibinomial
# --------------------------------------------------
msg("Fitting binomial and quasibinomial models ...")
all_res <- list()

for (cl in cluster_levels) {
  dat <- counts_df %>%
    filter(cluster == cl) %>%
    left_join(donor_meta, by = c("donor", "dx", "n_total"))

  dat <- as.data.frame(dat)
  dat$dx <- factor(dat$dx, levels = c("Control", "ASD"))
  if ("sex" %in% colnames(dat)) dat$sex <- factor(dat$sex)
  if ("batch" %in% colnames(dat)) dat$batch <- factor(dat$batch)

  # basic
  f_basic <- as.formula("cbind(n_cluster, n_other) ~ dx")
  fit_bin_basic   <- fit_glm_safe(f_basic, dat, binomial())
  fit_quasi_basic <- fit_glm_safe(f_basic, dat, quasibinomial())

  res_bin_basic <- if (!is.null(fit_bin_basic)) extract_coef(fit_bin_basic, "dxASD") else
    data.frame(term = "dxASD", beta = NA, se = NA, stat = NA, p = NA, OR = NA, conf.low = NA, conf.high = NA)

  res_quasi_basic <- if (!is.null(fit_quasi_basic)) extract_coef(fit_quasi_basic, "dxASD") else
    data.frame(term = "dxASD", beta = NA, se = NA, stat = NA, p = NA, OR = NA, conf.low = NA, conf.high = NA)

  res_bin_basic$model <- "binomial_basic"
  res_quasi_basic$model <- "quasibinomial_basic"

  # adjusted
  if (length(usable_covars) > 0) {
    f_adj <- as.formula(
      paste("cbind(n_cluster, n_other) ~ dx +", paste(usable_covars, collapse = " + "))
    )

    need_cols <- unique(c("n_cluster", "n_other", "dx", usable_covars))
    dat_adj <- dat[, need_cols, drop = FALSE]
    dat_adj <- dat_adj[complete.cases(dat_adj), , drop = FALSE]

    fit_bin_adj   <- NULL
    fit_quasi_adj <- NULL
    if (nrow(dat_adj) >= 8 && length(unique(dat_adj$dx)) == 2) {
      fit_bin_adj   <- fit_glm_safe(f_adj, dat_adj, binomial())
      fit_quasi_adj <- fit_glm_safe(f_adj, dat_adj, quasibinomial())
    }

    res_bin_adj <- if (!is.null(fit_bin_adj)) extract_coef(fit_bin_adj, "dxASD") else
      data.frame(term = "dxASD", beta = NA, se = NA, stat = NA, p = NA, OR = NA, conf.low = NA, conf.high = NA)

    res_quasi_adj <- if (!is.null(fit_quasi_adj)) extract_coef(fit_quasi_adj, "dxASD") else
      data.frame(term = "dxASD", beta = NA, se = NA, stat = NA, p = NA, OR = NA, conf.low = NA, conf.high = NA)

    res_bin_adj$model <- "binomial_adjusted"
    res_quasi_adj$model <- "quasibinomial_adjusted"
    res_bin_adj$n_complete <- nrow(dat_adj)
    res_quasi_adj$n_complete <- nrow(dat_adj)
  } else {
    res_bin_adj <- data.frame(
      term = "dxASD", beta = NA, se = NA, stat = NA, p = NA, OR = NA, conf.low = NA, conf.high = NA,
      model = "binomial_adjusted", n_complete = NA
    )
    res_quasi_adj <- data.frame(
      term = "dxASD", beta = NA, se = NA, stat = NA, p = NA, OR = NA, conf.low = NA, conf.high = NA,
      model = "quasibinomial_adjusted", n_complete = NA
    )
  }

  for (x in list(res_bin_basic, res_quasi_basic, res_bin_adj, res_quasi_adj)) {
    x$cluster <- cl
    x$mean_prop_ASD <- mean(dat$prop[dat$dx == "ASD"], na.rm = TRUE)
    x$mean_prop_Control <- mean(dat$prop[dat$dx == "Control"], na.rm = TRUE)
    x$effect_prop_ASD_minus_Control <- x$mean_prop_ASD - x$mean_prop_Control
    x$n_donors <- length(unique(dat$donor))
    x$n_ASD <- length(unique(dat$donor[dat$dx == "ASD"]))
    x$n_Control <- length(unique(dat$donor[dat$dx == "Control"]))
    x$covariates <- ifelse(length(usable_covars) == 0, "", paste(usable_covars, collapse = "+"))
    all_res[[length(all_res) + 1]] <- x
  }
}

res_df <- rbindlist(all_res, fill = TRUE)
res_df <- as.data.frame(res_df)

# 每个 model 单独做 FDR
res_df <- res_df %>%
  group_by(model) %>%
  mutate(fdr = p.adjust(p, method = "BH")) %>%
  ungroup()

fwrite(res_df, file.path(cfg$outdir, "Cluster_DonorLevel_GLM_Results.tsv"), sep = "\t")

# primary table：我建议 quasi basic 作为 primary，quasi adjusted 作为 sensitivity
primary_df <- res_df %>%
  filter(model == "quasibinomial_basic") %>%
  arrange(p, desc(effect_prop_ASD_minus_Control))
fwrite(primary_df, file.path(cfg$outdir, "Cluster_DonorLevel_GLM_Primary_quasiBasic.tsv"), sep = "\t")

sens_df <- res_df %>%
  filter(model == "quasibinomial_adjusted") %>%
  arrange(p, desc(effect_prop_ASD_minus_Control))
fwrite(sens_df, file.path(cfg$outdir, "Cluster_DonorLevel_GLM_Sensitivity_quasiAdjusted.tsv"), sep = "\t")

# --------------------------------------------------
# 5. 画图
# --------------------------------------------------
msg("Saving plots ...")

plot_df <- primary_df %>%
  mutate(
    cluster_label = paste0("C", cluster),
    neglog10p = -log10(pmax(p, 1e-300))
  )

plot_df$cluster_label <- factor(plot_df$cluster_label, levels = plot_df$cluster_label[order(plot_df$effect_prop_ASD_minus_Control)])

p_rank <- ggplot(plot_df, aes(x = effect_prop_ASD_minus_Control, y = cluster_label)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
  geom_point(aes(size = neglog10p)) +
  theme_classic(base_size = 12) +
  labs(
    title = "Discovery donor-level abundance model across microglial clusters",
    x = "Proportion difference (ASD - Control)",
    y = "Cluster",
    size = "-log10(p)"
  )

ggsave(file.path(cfg$outdir, "Fig1_DonorAbundance_quasiBasic_Ranking.pdf"), p_rank, width = 7, height = 6)

# target cluster boxplot
target_dat <- counts_df %>%
  filter(cluster == cfg$target_cluster) %>%
  mutate(dx = factor(dx, levels = c("Control", "ASD")))

if (nrow(target_dat) > 0) {
  target_row <- primary_df %>% filter(cluster == cfg$target_cluster)
  sub_lab <- if (nrow(target_row) > 0) {
    paste0(
      "quasibinomial p = ", signif(target_row$p[1], 3),
      "; OR = ", signif(target_row$OR[1], 3)
    )
  } else {
    "quasibinomial result unavailable"
  }

  p_target <- ggplot(target_dat, aes(x = dx, y = prop, fill = dx)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.75) +
    geom_jitter(width = 0.12, size = 2, alpha = 0.9) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 20, hjust = 1)
    ) +
    labs(
      title = paste0("Cluster ", cfg$target_cluster, " donor-level proportion"),
      subtitle = sub_lab,
      x = NULL,
      y = "Proportion within donor microglia"
    )

  ggsave(file.path(cfg$outdir, paste0("Fig1_Cluster", cfg$target_cluster, "_DonorProp_quasiBasic.pdf")),
         p_target, width = 5.6, height = 4.8)
}

save_session_info(file.path(cfg$outdir, "sessionInfo_Fig1_DonorAbundance_GLM.txt"))
msg("Done.")
