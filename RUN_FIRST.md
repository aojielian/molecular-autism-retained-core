# RUN_FIRST.md

## Minimal Reproducibility Guide

This repository contains the code used for the manuscript:

**Broad cortical microglial remodeling in autism spectrum disorder converges on a reproducible but non-discrete candidate/reference axis**

This file is intended to help reviewers and readers test reproducibility quickly, using the smallest set of frozen inputs from the Zenodo bundle.

---

## Repository layout

In the current repository version:

- analysis scripts are located under `analysis/`
- figure-generation scripts are located under `figures/`

Please ignore any older draft references to:

- `scripts/analysis/`
- `scripts/figures/`

Run all commands from the repository root directory.

---

## What to download first

Please download and extract the Zenodo reproducibility bundle before running any scripts.

Example local directory:

```text
/path/to/zenodo_bundle/
```

---

## Zenodo file mapping

The Zenodo archive contains the frozen input objects required for minimal reproducibility testing.

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

## Recommended minimal reproducibility tests

For the fastest validation-oriented reproducibility check, we recommend starting with:

1. **Cross-cohort validation in the Velmeshev single-nucleus dataset**
2. **External bulk-cortex validation in the Gandal 2022 dataset**

These two workflows reproduce the main validation layers of the revised manuscript without requiring the full discovery-stage pipeline to be rerun from raw inputs.

---

## Example 1. Cross-cohort validation (Velmeshev)

```bash
Rscript analysis/Package4_Velmeshev_v2_score_based_validation.R \
  --disc_rds /path/to/zenodo_bundle/Package2_strict_input_fixed_from_original.rds \
  --val_rds /path/to/zenodo_bundle/Velmeshev_Object.rds \
  --outdir results/Package4_CrossCohortValidation_Velmeshev_v2
```

This workflow supports:

- **Figure 5**
- **Supplementary Figure S3**
- **Supplementary Tables S6–S7**

Expected outputs typically include:

- donor-level validation summaries
- score-based support summaries
- mapped-cluster comparison summaries

---

## Example 2. External bulk-cortex validation (Gandal 2022)

```bash
Rscript analysis/Step08_Gandal2022_bulk_program_validation.R \
  --bulk_rdata /path/to/zenodo_bundle/Gene_NormalizedExpression_Metadata_wModelMatrix.RData \
  --candidate_genes /path/to/zenodo_bundle/candidate_genes.txt \
  --reference_genes /path/to/zenodo_bundle/reference_genes.txt \
  --outdir results/Step08_Gandal2022_bulk_program_validation
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

## Optional: discovery-stage audit

If desired, users may also inspect the discovery-stage retained-core framework.

```bash
Rscript analysis/package1_discovery_audit.R \
  --microglia_rds /path/to/zenodo_bundle/Microglia_Clustered.rds \
  --outdir results/Package1_Figure1_DiscoveryAudit
```

This workflow supports:

- **Figure 1**
- **Supplementary Figure S2**

---

## Interpretation of the minimal test

The revised manuscript is centered on the following reproducible layers:

1. **Retained cortical microglial core**
2. **Non-discrete candidate/reference axis**
3. **Weak mapped-cluster replication but stronger score-based program-level support**
4. **Independent bulk-cortex concordance**

Accordingly, reproducibility is best assessed first at the level of:

- donor-level validation
- score-based cross-cohort support
- external bulk-cortex validation

rather than expecting one-to-one recovery of a sharply discrete disease-specific subtype.

---

## If a script does not run

Please check the following first:

- Are the input file paths correct?
- Are you using the Zenodo files listed above?
- Are you running the command from the repository root directory?
- Are the required R packages installed?
- Does your local R environment match the environment documentation?

Please also consult:

- `README.md`
- `input_manifest.tsv`
- `environment/sessionInfo.txt`
- `Supplementary Table S10`

---

## Suggested reviewer workflow

For the fastest validation-oriented check, we suggest the following order:

1. Run **cross-cohort validation**
2. Run **bulk-cortex validation**
3. Optionally inspect the **discovery-stage audit**

This order reproduces the key inferential backbone of the revised manuscript with the least setup burden.
