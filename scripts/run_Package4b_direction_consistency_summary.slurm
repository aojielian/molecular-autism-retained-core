#!/bin/bash
#SBATCH -J Package4b_direction_consistency
#SBATCH -p CU
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 12:00:00
#SBATCH -o /gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/slurm_logs/%x-%j.out
#SBATCH -e /gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/slurm_logs/%x-%j.err

set -euo pipefail

PROJECT_DIR=/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision
PKG4_DIR=${PROJECT_DIR}/Package4_CrossCohortValidation_Velmeshev_v2
OUTDIR=${PROJECT_DIR}/Package4b_DirectionConsistencySummary
SCRIPT=${PROJECT_DIR}/scripts/Package4b_direction_consistency_summary.R

mkdir -p "${OUTDIR}"
mkdir -p "${PROJECT_DIR}/slurm_logs"


Rscript "${SCRIPT}" \
  --pkg4_dir "${PKG4_DIR}" \
  --outdir "${OUTDIR}"
