# Molecular Autism retained-core microglial remodeling analysis

This repository contains the analysis code supporting the main and supplementary results of the manuscript on **ASD-associated cortical microglial remodeling with a reproducible but non-discrete candidate/reference program**.

## Overview

The study was designed to address three linked questions:

1. Whether a conservative, donor-supported **retained cortical microglial core** could be defined from the discovery-stage ASD cortical microglial analysis space.
2. Whether the dominant ASD-associated signal was better interpreted as a **broad candidate/reference remodeling axis** rather than a sharply discrete disease-specific subtype.
3. Whether the same axis was supported across **cross-cohort single-nucleus validation** and **independent bulk-cortex validation**.

The analysis is organized into discovery, cross-cohort validation, external bulk validation, and supplementary mechanistic prioritization modules.

## Repository structure

```text
scripts/
├── analysis/
│   ├── package1_discovery_audit.R
│   ├── package1b_cluster_annotation_contaminant_audit.R
│   ├── package2_define_cluster2.R
│   ├── package2b_residual_background_audit.R
│   ├── Package3_discovery_donor_aware_disease_association.R
│   ├── Package4_Velmeshev_v2_score_based_validation.R
│   ├── Package4b_direction_consistency_manual_v2.R
│   ├── Package5_directional_concordance_global_summary.R
│   ├── Package6_validation_results_freezeout.R
│   ├── Package7_minimal_robustness_analysis.R
│   ├── Step08_Gandal2022_bulk_program_validation_v3.R
│   ├── Step08b_Gandal2022_leave_one_region_out_sensitivity.R
│   ├── drug_reversal_prioritization.py
│   ├── postprocess_drug_prioritization.py
│   └── tf_prioritization.py
└── figures/
    ├── Figure1.R
    ├── Figure2.R
    ├── Figure3.R
    ├── Figure4.R
    ├── Figure5.R
    ├── Figure6.R
    ├── Supplementary_Figure_S2.R
    ├── Supplementary_Figure_S3.R
    ├── Supplementary_Figure_S4.R
    └── Supplementary_Figure_S5.R
```

Additional supporting resources such as curated marker panels, frozen candidate/reference programs, session information, and supplementary-table build scripts can be added under `metadata/`, `environment/`, and `docs/` as needed.

## Analytical modules and manuscript mapping

### Discovery retained-core definition

- `package1_discovery_audit.R`
- `package1b_cluster_annotation_contaminant_audit.R`
- `package2_define_cluster2.R`
- `package2b_residual_background_audit.R`

These scripts define the discovery-stage cortical microglial analysis space, perform all-cluster auditing, identify contaminant/background structure, and freeze the retained cortical microglial core.

**Primary manuscript mapping:**
- Figure 1
- Supplementary Figure S2
- Supplementary Table S2

### Discovery donor-aware disease analyses

- `Package3_discovery_donor_aware_disease_association.R`

This script performs donor-level composition analyses, candidate/reference score summaries, donor-paired internal state testing, and pseudobulk case-control analyses across the retained cortical microglial core.

**Primary manuscript mapping:**
- Figure 2
- Figure 3
- Figure 4
- Supplementary Tables S4-S5

### Cross-cohort validation in Velmeshev

- `Package4_Velmeshev_v2_score_based_validation.R`
- `Package4b_direction_consistency_manual_v2.R`
- `Package5_directional_concordance_global_summary.R`
- `Package6_validation_results_freezeout.R`
- `Package7_minimal_robustness_analysis.R`

These scripts implement validation-stage microglial reclustering, mapped-cluster comparison to discovery retained states, frozen candidate/reference program scoring, seven prespecified donor-level validation metrics, aggregate directional concordance testing, and minimal robustness analyses.

**Primary manuscript mapping:**
- Figure 5
- Supplementary Figure S3
- Supplementary Tables S6-S7

### External bulk-cortex validation in Gandal 2022

- `Step08_Gandal2022_bulk_program_validation_v3.R`
- `Step08b_Gandal2022_leave_one_region_out_sensitivity.R`

These scripts project the frozen candidate/reference programs into the independent Gandal 2022 bulk-cortex dataset, test ASD effects for candidate, reference, and candidate-minus-reference scores, and perform leave-one-region-out sensitivity analyses.

**Primary manuscript mapping:**
- Figure 6
- Supplementary Tables S8-S9

### Supplementary mechanistic prioritization

- `drug_reversal_prioritization.py`
- `postprocess_drug_prioritization.py`
- `tf_prioritization.py`

These scripts generate supplementary mechanistic prioritization layers, including perturbation-based reversal ranking and TF-target enrichment.

**Primary manuscript mapping:**
- Supplementary Figure S4
- Supplementary Figure S5

## Frozen analytical definitions

The following definitions are treated as fixed, manuscript-level analytical specifications rather than ad hoc post hoc choices:

- **Retained cortical microglial core:** discovery clusters **0** and **2**
- **Candidate side of retained-core contrast:** cluster **2**
- **Reference-like side of retained-core contrast:** cluster **0**
- **Validation reclustering resolution:** **0.4**
- **Validation PCs used:** **20**
- **Candidate quantile threshold:** **0.75**
- **Reference quantile threshold:** **0.25**
- **Mapped-cluster correlation delta threshold:** **0.05**
- **Mapped-cluster state delta threshold:** **0.00**
- **Minimum cells per validation sample:** **20**
- **Minimum cells per validation group:** **10**

These frozen definitions are also described in the manuscript Methods and summarized in Supplementary Table S10.

## Input data

This repository is intended for analysis of publicly available, de-identified transcriptomic datasets used in the manuscript:

- Discovery cortical single-cell resource: PsychENCODE/Wamsley
- Cross-cohort validation single-nucleus dataset: Velmeshev ASD cortex
- External bulk-cortex validation dataset: Gandal 2022

Large primary input matrices and raw count objects are not included in this repository unless explicitly permitted by the source datasets. Users should obtain the original resources from the corresponding public sources cited in the manuscript.

## Running the analysis

Most R scripts are designed to run through command-line arguments or environment variables rather than local hard-coded paths.

Example environment setup:

```bash
export MA_PROJECT_DIR=/path/to/project
export MA_MICROGLIA_RDS=/path/to/discovery_microglia_object.rds
export MA_STRICT_RDS=/path/to/strict_retained_core_object.rds
export MA_GLOBAL_RDS=/path/to/global_reference_object.rds
export MA_VELMESH_RDS=/path/to/velmeshev_microglia_object.rds
```

Example command:

```bash
Rscript scripts/analysis/package1_discovery_audit.R \
  --microglia_rds "$MA_MICROGLIA_RDS" \
  --outdir "$MA_PROJECT_DIR/results/Package1_Figure1_DiscoveryAudit"
```

Another example:

```bash
Rscript scripts/analysis/Package4_Velmeshev_v2_score_based_validation.R \
  --strict_rds "$MA_STRICT_RDS" \
  --val_rds "$MA_VELMESH_RDS" \
  --outdir "$MA_PROJECT_DIR/results/Package4_CrossCohortValidation_Velmeshev_v2"
```

Python scripts provide `--help` messages and can be run similarly after adapting input paths.

## Figure generation

The `scripts/figures/` directory contains manuscript-level plotting and figure-assembly scripts corresponding to Figure 1-6 and Supplementary Figure S2-S5.

Recommended workflow:

1. Run analysis scripts to generate intermediate tables and summaries.
2. Run figure scripts using the corresponding outputs.
3. Build supplementary tables separately using the supplementary-table workbook builder.

## Reproducibility notes

- Use donor-aware inference rather than cell-level pseudoreplicated testing for case-control claims.
- Interpret mapped-cluster validation and program-level validation as related but distinct analytical layers.
- Treat supplementary perturbation and TF analyses as **hypothesis-generating**, not as causal or therapeutic proof.

## Software environment

The manuscript reports the main software environment, package versions, and key analysis parameters in Supplementary Table S10. A repository-level environment file and/or session information file should be added under `environment/`.

## Code availability

This repository is intended to support the code availability statement in the manuscript. A versioned GitHub release should ideally be archived through Zenodo so that a DOI can be cited in the final publication.
