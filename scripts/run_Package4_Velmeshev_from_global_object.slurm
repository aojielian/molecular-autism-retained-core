#!/bin/bash
#SBATCH -J MA_Pkg4_VelGlobal
#SBATCH -p CU
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=12
#SBATCH --mem=180G
#SBATCH -t 24:00:00
#SBATCH -o /gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/slurm_logs/%x-%j.out
#SBATCH -e /gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/slurm_logs/%x-%j.err

set -euo pipefail

R_SCRIPT="/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/scripts/Package4_Velmeshev_from_global_object.R"
DISC_RDS="/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Package2_strict_input_fixed_from_original.rds"
VAL_RDS="/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/Velmeshev_scRNA/Velmeshev_Object.rds"
OUTDIR="/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Package4_CrossCohortValidation_Velmeshev"

mkdir -p "${OUTDIR}"


export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export VECLIB_MAXIMUM_THREADS=${SLURM_CPUS_PER_TASK}
export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK}

Rscript "${R_SCRIPT}" \
  --disc_rds "${DISC_RDS}" \
  --val_rds "${VAL_RDS}" \
  --outdir "${OUTDIR}" \
  --disc_assay RNA \
  --val_assay RNA \
  --disc_cluster_col cluster_use \
  --disc_donor_col individual_ID \
  --disc_dx_col Diagnosis \
  --val_global_cluster_col cluster \
  --val_donor_col AUTO \
  --val_dx_col diagnosis \
  --val_microglia_label Microglia \
  --case_label ASD \
  --control_label Control \
  --disc_cluster0 0 \
  --disc_cluster2 2 \
  --min_cells_donor 20 \
  --min_cells_cluster 10 \
  --top_n_pos 100 \
  --top_n_neg 100 \
  --val_npcs 20 \
  --val_resolution 0.4 \
  --mapping_corr_delta_threshold 0.05 \
  --mapping_state_delta_threshold 0.00 \
  --alpha 0.05 \
  --threads ${SLURM_CPUS_PER_TASK}
