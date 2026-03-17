#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

root <- "/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision"
outdir <- file.path(root, "Supplementary_Tables_v1")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---------------- helpers ----------------
safe_read <- function(path) {
  if (!file.exists(path)) return(NULL)
  ext_gz <- grepl("\\.gz$", path)
  dt <- tryCatch({
    if (ext_gz) fread(cmd = paste("gzip -dc", shQuote(path)), sep = "\t", header = TRUE, data.table = TRUE)
    else fread(path, sep = "\t", header = TRUE, data.table = TRUE)
  }, error = function(e) NULL)
  dt
}

safe_write <- function(dt, filename) {
  fwrite(dt, file.path(outdir, filename), sep = "\t")
}

txt_write <- function(lines, filename) {
  writeLines(lines, file.path(outdir, filename))
}

find_first <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

add_source <- function(dt, source_name) {
  if (is.null(dt)) return(NULL)
  dt[, source_file := source_name]
  dt
}

# ---------------- package paths ----------------
pkg1  <- file.path(root, "Package1_Figure1_DiscoveryAudit")
pkg1b <- file.path(root, "Package1b_ClusterAnnotation_ContaminantAudit")
pkg2  <- file.path(root, "Package2_Define_Cluster2")
pkg2b <- file.path(root, "Package2b_ResidualBackgroundAudit")
pkg3  <- file.path(root, "Package3_Discovery_DonorAwareDiseaseAssociation")
pkg4v2 <- file.path(root, "Package4_CrossCohortValidation_Velmeshev_v2")
pkg4b <- file.path(root, "Package4b_DirectionConsistencySummary_v2_manual")
pkg5  <- root
pkg7  <- root
pkg8  <- file.path(root, "Package8_ExternalValidation_Gandal2022")
pkg8b <- file.path(root, "Package8b_Gandal2022_LOO_Sensitivity")

# ---------------- Supplementary Table S1 ----------------
# Cohort summary assembled from discovery / validation / external validation outputs
s1_rows <- list()

# Discovery
disc_audit <- safe_read(file.path(pkg1, "Package1_Audit_Metadata.tsv"))
disc_donor_prop <- safe_read(file.path(pkg1, "06_allcluster_donor_proportions.tsv.gz"))
disc_abund <- safe_read(file.path(pkg1, "07_allcluster_abundance_stats.tsv"))

disc_row <- data.table(
  cohort = "Discovery",
  source_primary = "Package1_Figure1_DiscoveryAudit",
  n_samples = NA_integer_,
  n_donors = NA_integer_,
  n_cells = NA_integer_,
  n_ASD = NA_integer_,
  n_Control = NA_integer_,
  notes = NA_character_
)

if (!is.null(disc_donor_prop)) {
  donor_col <- intersect(c("donor","sample","Subject","subject","individual_ID","iid"), names(disc_donor_prop))
  dx_col <- intersect(c("dx","Diagnosis","diagnosis","group","Group","dx_std"), names(disc_donor_prop))
  if (length(donor_col) > 0) disc_row$n_donors <- uniqueN(disc_donor_prop[[donor_col[1]]])
  if (length(dx_col) > 0) {
    vals <- as.character(disc_donor_prop[[dx_col[1]]])
    disc_row$n_ASD <- sum(vals %in% c("ASD","Case","case","Autism"), na.rm = TRUE)
    disc_row$n_Control <- sum(vals %in% c("Control","CTL","control"), na.rm = TRUE)
  }
}
if (!is.null(disc_audit)) {
  if ("n_cells" %in% names(disc_audit)) disc_row$n_cells <- sum(as.numeric(disc_audit$n_cells), na.rm = TRUE)
  if ("n_samples" %in% names(disc_audit)) disc_row$n_samples <- sum(as.numeric(disc_audit$n_samples), na.rm = TRUE)
}
if (!is.null(disc_abund)) {
  disc_row$notes <- "Discovery counts assembled from Package1 audit outputs."
}
s1_rows[[length(s1_rows)+1]] <- disc_row

# Velmeshev validation
val_global <- safe_read(file.path(pkg4v2, "02_validation_global_dx_counts.tsv"))
val_micro  <- safe_read(file.path(pkg4v2, "03_validation_microglia_dx_counts.tsv"))
val_row <- data.table(
  cohort = "Velmeshev validation",
  source_primary = "Package4_CrossCohortValidation_Velmeshev_v2",
  n_samples = NA_integer_,
  n_donors = NA_integer_,
  n_cells = NA_integer_,
  n_ASD = NA_integer_,
  n_Control = NA_integer_,
  notes = "Microglial-cell counts from validation cohort."
)
if (!is.null(val_global) && all(c("dx","N") %in% names(val_global))) {
  val_row$n_ASD <- val_global[dx == "ASD", sum(N)]
  val_row$n_Control <- val_global[dx == "Control", sum(N)]
}
if (!is.null(val_micro) && all(c("dx","N") %in% names(val_micro))) {
  val_row$n_cells <- sum(val_micro$N, na.rm = TRUE)
}
s1_rows[[length(s1_rows)+1]] <- val_row

# Gandal bulk
gandal_dx <- safe_read(file.path(pkg8, "03_filtered_dx_counts.tsv"))
gandal_scores <- safe_read(file.path(pkg8, "07_sample_level_program_scores.tsv"))
gandal_row <- data.table(
  cohort = "Gandal 2022 bulk cortex",
  source_primary = "Package8_ExternalValidation_Gandal2022",
  n_samples = NA_integer_,
  n_donors = NA_integer_,
  n_cells = NA_integer_,
  n_ASD = NA_integer_,
  n_Control = NA_integer_,
  notes = "Bulk cortex external validation cohort."
)
if (!is.null(gandal_dx) && all(c("dx","N") %in% names(gandal_dx))) {
  gandal_row$n_ASD <- gandal_dx[dx == "ASD", sum(N)]
  gandal_row$n_Control <- gandal_dx[dx == "Control", sum(N)]
}
if (!is.null(gandal_scores)) {
  gandal_row$n_samples <- nrow(gandal_scores)
  subj_col <- intersect(c("Subject","subject","donor","sample"), names(gandal_scores))
  if (length(subj_col) > 0) gandal_row$n_donors <- uniqueN(gandal_scores[[subj_col[1]]])
}
s1_rows[[length(s1_rows)+1]] <- gandal_row

s1 <- rbindlist(s1_rows, fill = TRUE)
safe_write(s1, "Supplementary_Table_S1_cohort_summary.tsv")

# ---------------- Supplementary Table S2 ----------------
# Discovery audit and exclusion rationale
s2_parts <- list(
  add_source(safe_read(file.path(pkg1, "10_marker_support.tsv")), "Package1/10_marker_support.tsv"),
  add_source(safe_read(file.path(pkg1, "11_cluster_prioritization.tsv")), "Package1/11_cluster_prioritization.tsv"),
  add_source(safe_read(file.path(pkg1, "12_top_candidate.tsv")), "Package1/12_top_candidate.tsv"),
  add_source(safe_read(file.path(pkg1, "13_cluster2_summary.tsv")), "Package1/13_cluster2_summary.tsv")
)
s2 <- rbindlist(Filter(Negate(is.null), s2_parts), fill = TRUE, use.names = TRUE)
safe_write(s2, "Supplementary_Table_S2_discovery_audit_and_exclusion_rationale.tsv")

# ---------------- Supplementary Table S3 ----------------
c2 <- add_source(safe_read(file.path(pkg2, "tables", "22_cluster2_core_summary_up.tsv")), "Package2/tables/22_cluster2_core_summary_up.tsv")
c0 <- add_source(safe_read(file.path(pkg2, "tables", "23_cluster0_core_summary_up.tsv")), "Package2/tables/23_cluster0_core_summary_up.tsv")
if (!is.null(c2)) c2[, program := "candidate_cluster2"]
if (!is.null(c0)) c0[, program := "reference_cluster0"]
s3 <- rbindlist(list(c2, c0), fill = TRUE, use.names = TRUE)
safe_write(s3, "Supplementary_Table_S3_frozen_candidate_reference_gene_programs.tsv")

# ---------------- Supplementary Table S4 ----------------
s4 <- safe_read(file.path(pkg3, "13_pseudobulk_all_retained_ASD_vs_Control.tsv.gz"))
safe_write(s4, "Supplementary_Table_S4_pooled_retained_microglia_DEGs.tsv")

# ---------------- Supplementary Table S5 ----------------
s5a <- safe_read(file.path(pkg3, "15_pseudobulk_cluster0_ASD_vs_Control.tsv.gz"))
s5b <- safe_read(file.path(pkg3, "17_pseudobulk_cluster2_ASD_vs_Control.tsv.gz"))
if (!is.null(s5a)) s5a[, cluster_context := "cluster0"]
if (!is.null(s5b)) s5b[, cluster_context := "cluster2"]
s5 <- rbindlist(list(s5a, s5b), fill = TRUE, use.names = TRUE)
safe_write(s5, "Supplementary_Table_S5_cluster0_cluster2_DEGs.tsv")

# ---------------- Supplementary Table S6 ----------------
s6_files <- list(
  add_source(safe_read(file.path(pkg4b, "00_shared_comparable_samples.tsv")), "Package4b/00_shared_comparable_samples.tsv"),
  add_source(safe_read(file.path(pkg4b, "00b_shared_sample_counts.tsv")), "Package4b/00b_shared_sample_counts.tsv"),
  add_source(safe_read(file.path(pkg4b, "03_direction_consistency_summary_pretty.tsv")), "Package4b/03_direction_consistency_summary_pretty.tsv"),
  add_source(safe_read(file.path(pkg4b, "04_direction_consistency_summary_manuscript.tsv")), "Package4b/04_direction_consistency_summary_manuscript.tsv")
)
s6 <- rbindlist(Filter(Negate(is.null), s6_files), fill = TRUE, use.names = TRUE)
safe_write(s6, "Supplementary_Table_S6_full_validation_metric_summary.tsv")

# ---------------- Supplementary Table S7 ----------------
s7_files <- list(
  add_source(safe_read(file.path(pkg5, "Package5_DirectionalConcordance_GlobalSummary", "01_direction_concordance_master.tsv")), "Package5/01_direction_concordance_master.tsv"),
  add_source(safe_read(file.path(pkg5, "Package5_DirectionalConcordance_GlobalSummary", "02_global_direction_tests.tsv")), "Package5/02_global_direction_tests.tsv"),
  add_source(safe_read(file.path(pkg5, "Package5_DirectionalConcordance_GlobalSummary", "03_metrics_ranked_by_nominal_p.tsv")), "Package5/03_metrics_ranked_by_nominal_p.tsv"),
  add_source(safe_read(file.path(pkg5, "Package5_DirectionalConcordance_GlobalSummary", "04_metrics_ranked_by_aligned_effect.tsv")), "Package5/04_metrics_ranked_by_aligned_effect.tsv"),
  add_source(safe_read(file.path(pkg7, "Package7_MinimalRobustnessAnalysis", "01_score_defined_harmonization_sensitivity.tsv")), "Package7/01_score_defined_harmonization_sensitivity.tsv"),
  add_source(safe_read(file.path(pkg7, "Package7_MinimalRobustnessAnalysis", "02_score_defined_direction_stability.tsv")), "Package7/02_score_defined_direction_stability.tsv"),
  add_source(safe_read(file.path(pkg7, "Package7_MinimalRobustnessAnalysis", "03_leave_one_metric_out_global_summary.tsv")), "Package7/03_leave_one_metric_out_global_summary.tsv")
)
s7 <- rbindlist(Filter(Negate(is.null), s7_files), fill = TRUE, use.names = TRUE)
safe_write(s7, "Supplementary_Table_S7_directional_concordance_and_robustness.tsv")

# ---------------- Supplementary Table S8 ----------------
s8_files <- list(
  add_source(safe_read(file.path(pkg8, "04_candidate_genes_used.txt")), "Package8/04_candidate_genes_used.txt"),
  add_source(safe_read(file.path(pkg8, "04_reference_genes_used.txt")), "Package8/04_reference_genes_used.txt"),
  add_source(safe_read(file.path(pkg8, "06_candidate_hit_genes.txt")), "Package8/06_candidate_hit_genes.txt"),
  add_source(safe_read(file.path(pkg8, "06_reference_hit_genes.txt")), "Package8/06_reference_hit_genes.txt"),
  add_source(safe_read(file.path(pkg8, "06_gene_set_overlap.tsv")), "Package8/06_gene_set_overlap.tsv")
)
# text files need special handling
txt_to_dt <- function(path, label) {
  if (!file.exists(path)) return(NULL)
  x <- readLines(path, warn = FALSE)
  data.table(gene = x, source_file = label)
}
s8 <- rbindlist(list(
  txt_to_dt(file.path(pkg8, "04_candidate_genes_used.txt"), "Package8/04_candidate_genes_used.txt"),
  txt_to_dt(file.path(pkg8, "04_reference_genes_used.txt"), "Package8/04_reference_genes_used.txt"),
  txt_to_dt(file.path(pkg8, "06_candidate_hit_genes.txt"), "Package8/06_candidate_hit_genes.txt"),
  txt_to_dt(file.path(pkg8, "06_reference_hit_genes.txt"), "Package8/06_reference_hit_genes.txt"),
  add_source(safe_read(file.path(pkg8, "06_gene_set_overlap.tsv")), "Package8/06_gene_set_overlap.tsv")
), fill = TRUE, use.names = TRUE)
safe_write(s8, "Supplementary_Table_S8_Gandal_gene_mapping_summary.tsv")

# ---------------- Supplementary Table S9 ----------------
s9_files <- list(
  add_source(safe_read(file.path(pkg8b, "04_leave_one_region_out_results.tsv")), "Package8b/04_leave_one_region_out_results.tsv"),
  add_source(safe_read(file.path(pkg8b, "05_leave_one_region_out_vs_baseline.tsv")), "Package8b/05_leave_one_region_out_vs_baseline.tsv"),
  add_source(safe_read(file.path(pkg8b, "06_leave_one_region_out_stability_summary.tsv")), "Package8b/06_leave_one_region_out_stability_summary.tsv")
)
s9 <- rbindlist(Filter(Negate(is.null), s9_files), fill = TRUE, use.names = TRUE)
safe_write(s9, "Supplementary_Table_S9_leave_one_region_out_results.tsv")

# ---------------- Supplementary Table S10 ----------------
# SessionInfo and run parameters summary
session_files <- c(
  file.path(pkg1, "sessionInfo_package1.txt"),
  file.path(pkg1b, "sessionInfo_package1b.txt"),
  file.path(pkg2, "sessionInfo_package2.txt"),
  file.path(pkg2b, "sessionInfo_package2b.txt"),
  file.path(pkg4v2, "31_sessionInfo.txt"),
  file.path(root, "Package5_DirectionalConcordance_GlobalSummary", "06_manuscript_summary.txt"),
  file.path(root, "Package7_MinimalRobustnessAnalysis", "04_package7_summary.txt"),
  file.path(pkg8, "13_sessionInfo.txt"),
  file.path(pkg8b, "11_sessionInfo.txt")
)
session_dt <- rbindlist(lapply(session_files[file.exists(session_files)], function(f) {
  data.table(file = basename(f), line_no = seq_along(readLines(f, warn = FALSE)), text = readLines(f, warn = FALSE))
}), fill = TRUE)
safe_write(session_dt, "Supplementary_Table_S10_session_and_parameter_summary.tsv")

# ---------------- Manifest ----------------
manifest <- data.table(
  supplementary_table = paste0("S", 1:10),
  filename = c(
    "Supplementary_Table_S1_cohort_summary.tsv",
    "Supplementary_Table_S2_discovery_audit_and_exclusion_rationale.tsv",
    "Supplementary_Table_S3_frozen_candidate_reference_gene_programs.tsv",
    "Supplementary_Table_S4_pooled_retained_microglia_DEGs.tsv",
    "Supplementary_Table_S5_cluster0_cluster2_DEGs.tsv",
    "Supplementary_Table_S6_full_validation_metric_summary.tsv",
    "Supplementary_Table_S7_directional_concordance_and_robustness.tsv",
    "Supplementary_Table_S8_Gandal_gene_mapping_summary.tsv",
    "Supplementary_Table_S9_leave_one_region_out_results.tsv",
    "Supplementary_Table_S10_session_and_parameter_summary.tsv"
  )
)
safe_write(manifest, "Supplementary_Tables_manifest.tsv")

message("Supplementary tables written to: ", outdir)


