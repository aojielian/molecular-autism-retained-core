# Molecular Autism retained-core microglial remodeling analysis

This repository contains the analysis code supporting the main and supplementary results of the manuscript:

**Broad cortical microglial remodeling in autism spectrum disorder converges on a reproducible but non-discrete candidate/reference axis**

The revised manuscript is centered on four linked inferential layers:

1. definition of a **retained cortical microglial core** from the discovery-stage cortical microglial landscape
2. identification of a **non-discrete candidate/reference remodeling axis** within that retained core
3. **cross-cohort validation** in the Velmeshev single-nucleus dataset, with weak mapped-cluster replication but stronger score-based program-level support
4. **external bulk-cortex validation** in the Gandal 2022 cohort

Supplementary perturbation and TF analyses are included as **supportive mechanistic context** and should not be interpreted as causal or therapeutic proof.

---

## Repository structure

```text
scripts/
└── analysis/
    ├── package1_discovery_audit.R
    ├── package1b_cluster_annotation_contaminant_audit.R
    ├── package2_define_cluster2.R
    ├── package2b_residual_background_audit.R
    ├── Package3_discovery_donor_aware_disease_association.R
    ├── Package4_Velmeshev_v2_score_based_validation.R
    ├── Package4b_direction_consistency.R
    ├── Package5_directional_concordance_global_summary.R
    ├── Package6_validation_results_freezeout.R
    ├── Package7_minimal_robustness_analysis.R
    ├── Step08_Gandal2022_bulk_program_validation.R
    ├── Step08b_Gandal2022_leave_one_region_out_sensitivity.R
    ├── drug_reversal_prioritization.py
    ├── postprocess_drug_prioritization.py
    ├── tf_prioritization.py
    └── README_molecular_autism_retained_core.md
```

If figure-generation scripts are included separately in your local project, they can be documented in a companion `figures/` directory, but the main reproducibility logic of the revised manuscript is captured by the analysis scripts listed above.

---

## Analytical modules and manuscript mapping

### 1. Discovery retained-core definition

Scripts:
- `package1_discovery_audit.R`
- `package1b_cluster_annotation_contaminant_audit.R`
- `package2_define_cluster2.R`
- `package2b_residual_background_audit.R`

Purpose:
- define the discovery-stage cortical microglial analysis space
- perform all-cluster auditing
- identify contaminant/background structure
- freeze the retained cortical microglial core

Main manuscript mapping:
- **Figure 1**
- **Supplementary Figure S2**
- **Supplementary Table S2**

---

### 2. Discovery donor-aware disease analyses

Script:
- `Package3_discovery_donor_aware_disease_association.R`

Purpose:
- donor-level composition analyses
- candidate/reference score summaries
- donor-paired internal state testing
- pooled retained-core pseudobulk case-control analyses

Main manuscript mapping:
- **Figure 2**
- **Figure 3**
- **Figure 4**
- **Supplementary Tables S4–S5**

---

### 3. Cross-cohort validation in Velmeshev

Scripts:
- `Package4_Velmeshev_v2_score_based_validation.R`
- `Package4b_direction_consistency.R`
- `Package5_directional_concordance_global_summary.R`
- `Package6_validation_results_freezeout.R`
- `Package7_minimal_robustness_analysis.R`

Purpose:
- validation-stage microglial comparison to discovery retained states
- frozen candidate/reference program scoring
- donor-level validation metrics
- aggregate directional concordance testing
- minimal robustness analyses

Main manuscript mapping:
- **Figure 5**
- **Supplementary Figure S3**
- **Supplementary Tables S6–S7**

---

### 4. External bulk-cortex validation in Gandal 2022

Scripts:
- `Step08_Gandal2022_bulk_program_validation.R`
- `Step08b_Gandal2022_leave_one_region_out_sensitivity.R`

Purpose:
- project discovery-derived candidate/reference programs into the Gandal 2022 bulk-cortex dataset
- test ASD effects for candidate, reference, and candidate-minus-reference scores
- perform leave-one-region-out sensitivity analyses

Main manuscript mapping:
- **Figure 6**
- **Supplementary Tables S8–S9**

---

### 5. Supplementary mechanistic prioritization

Scripts:
- `drug_reversal_prioritization.py`
- `postprocess_drug_prioritization.py`
- `tf_prioritization.py`

Purpose:
- perturbation-based reversal ranking
- TF-target enrichment
- supplementary mechanistic prioritization compatible with the retained-core remodeling framework

Main manuscript mapping:
- **Supplementary Figure S4**
- **Supplementary Figure S5**

---

## Frozen analytical definitions

The following are treated as fixed manuscript-level analytical specifications rather than ad hoc post hoc choices:

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

These frozen definitions are also summarized in the manuscript Methods and Supplementary Table S10.

---

## Input data and Zenodo reproducibility bundle

This repository is intended for analysis of publicly available, de-identified transcriptomic datasets used in the manuscript:

- discovery cortical microglial resource: **PsychENCODE / Wamsley**
- cross-cohort validation dataset: **Velmeshev single-nucleus ASD cortex**
- external bulk-cortex validation dataset: **Gandal 2022**

Large primary input matrices are not distributed in this GitHub repository unless explicitly permitted by the source datasets. For reproducibility testing, frozen intermediate objects and supporting files are provided in the associated Zenodo bundle.

### Zenodo file mapping

| Zenodo file | Description | Typical use |
|---|---|---|
| `Microglia_Clustered.rds` | discovery-stage microglia object | discovery/audit analyses |
| `Package2_strict_input_fixed_from_original.rds` | strict retained-core input object | donor-level retained-core and validation analyses |
| `Velmeshev_Object.rds` | validation microglia object | cross-cohort validation |
| `08_discovery_frozen_state_signature.tsv` | frozen candidate/reference signature table | supporting signature-based analyses |
| `candidate_genes.txt` | candidate program gene list | bulk validation |
| `reference_genes.txt` | reference program gene list | bulk validation |
| `13_pseudobulk_all_retained_ASD_vs_Control.tsv.gz` | pooled retained-core pseudobulk DE results | supporting outputs |
| `Gene_NormalizedExpression_Metadata_wModelMatrix.RData` | Gandal 2022 bulk validation input | bulk validation |
| `input_manifest.tsv` | manifest linking frozen inputs to manuscript analyses | reproducibility guidance |

---

## Minimal reproducibility guide

For the fastest validation-oriented reproducibility check, we recommend starting with:

1. **Cross-cohort validation in the Velmeshev single-nucleus dataset**
2. **External bulk-cortex validation in the Gandal 2022 dataset**

These two workflows reproduce the main validation layers of the revised manuscript without requiring the full discovery-stage pipeline to be rerun from raw inputs.

### Important usage note

The examples below use **absolute paths**.  
Users may run the scripts from any local working directory, provided that the script path and all input/output paths are explicitly specified.

---

## Example 1. Cross-cohort validation (Velmeshev)

```bash
Rscript /absolute/path/to/scripts/analysis/Package4_Velmeshev_v2_score_based_validation.R \
  --disc_rds /absolute/path/to/zenodo_bundle/Package2_strict_input_fixed_from_original.rds \
  --val_rds /absolute/path/to/zenodo_bundle/Velmeshev_Object.rds \
  --outdir /absolute/path/to/output/results/Package4_CrossCohortValidation_Velmeshev_v2
```

This workflow supports:

- **Figure 5**
- **Supplementary Figure S3**
- **Supplementary Tables S6–S7**

Expected outputs typically include:

- donor-level validation summaries
- score-based support summaries
- mapped-cluster comparison summaries
- directional concordance summaries

---

## Example 2. External bulk-cortex validation (Gandal 2022)

```bash
Rscript /absolute/path/to/scripts/analysis/Step08_Gandal2022_bulk_program_validation.R \
  --bulk_rdata /absolute/path/to/zenodo_bundle/Gene_NormalizedExpression_Metadata_wModelMatrix.RData \
  --candidate_genes /absolute/path/to/zenodo_bundle/candidate_genes.txt \
  --reference_genes /absolute/path/to/zenodo_bundle/reference_genes.txt \
  --outdir /absolute/path/to/output/results/Step08_Gandal2022_bulk_program_validation
```

This workflow supports:

- **Figure 6**
- **Supplementary Tables S8–S9**

Expected outputs typically include:

- candidate score results
- reference score results
- candidate-minus-reference composite results
- leave-one-region-out sensitivity summaries

---

## Optional example 3. Discovery-stage audit

If desired, users may also inspect the discovery-stage retained-core framework.

```bash
Rscript /absolute/path/to/scripts/analysis/package1_discovery_audit.R \
  --microglia_rds /absolute/path/to/zenodo_bundle/Microglia_Clustered.rds \
  --outdir /absolute/path/to/output/results/Package1_Figure1_DiscoveryAudit
```

This workflow supports:

- **Figure 1**
- **Supplementary Figure S2**

---

## Interpretation of the minimal test

The revised manuscript is centered on the following reproducible layers:

1. **retained cortical microglial core**
2. **non-discrete candidate/reference axis**
3. **weak mapped-cluster replication but stronger score-based program-level support**
4. **independent bulk-cortex concordance**

Accordingly, reproducibility is best assessed first at the level of:

- donor-level validation
- score-based cross-cohort support
- external bulk-cortex validation

rather than expecting one-to-one recovery of a sharply discrete disease-specific subtype.

---

## If a script does not run

Please check the following first:

- are the absolute file paths correct?
- are you using the Zenodo files listed above?
- are the required R packages installed?
- does your local environment match the environment documentation?
- are output directories writable?

Please also consult:

- `input_manifest.tsv`
- `Supplementary Table S10`
- any environment/session information files included with the repository

---

## Reproducibility notes

- Use **donor-aware inference** rather than cell-level pseudoreplicated testing for case-control claims.
- Interpret **mapped-cluster validation** and **program-level validation** as related but distinct analytical layers.
- Treat **perturbation** and **TF-target** analyses as **hypothesis-generating** supportive context, not as causal or therapeutic proof.

---

## Software environment

The manuscript reports the main software environment, package versions, and key analysis parameters in **Supplementary Table S10**.

If available, repository-level environment/session information files should also be provided alongside the repository, for example:

- `environment/sessionInfo.txt`
- `environment/R_packages_versions.tsv`

---

## Code availability

This repository is intended to support the code availability statement in the manuscript.  
A versioned GitHub snapshot should ideally be archived through Zenodo so that a DOI-linked code/data package can be cited in the final publication.
