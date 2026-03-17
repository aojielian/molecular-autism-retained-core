# MA-risk-microglia-ASD

This repository contains the analysis code, documentation, and reproducibility materials for the Molecular Autism revision project on ASD-associated SPP1+ lipid-metabolic microglial states.

## Repository purpose
This local repository is used for:
- organizing figure-generation scripts
- documenting data sources and paths
- tracking manuscript-related code revisions
- preparing materials for GitHub/Zenodo release

## Important note
Raw input data are **not stored in this repository**.
The analyses are run on the HPC server using the original paths under `/gpfs/...`.

## Main server-side data locations
- PsychENCODE discovery cohort:
  `/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/PsychENCODE_Science_2024/`
- Velmeshev validation cohort:
  `/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/Velmeshev_scRNA/`
- Bulk RNA-seq cohorts:
  `/gpfs/hpc/home/lijc/lianaoj/autism_bulk_RNA/`

## Planned script modules
- `01_fig1_discovery_clusters.R`
- `02_fig2_validation_donor_level.R`
- `03_fig3_cellchat_bulk.R`
- `04_fig4_metabolism.R`
- `05_fig5_tf_activity.R`
- `06_supp_tables.R`

## Reproducibility
This repository will eventually include:
- session information
- software/package versions
- curated gene-set definitions
- supplementary table generation code
