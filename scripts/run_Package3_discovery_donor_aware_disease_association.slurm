#!/bin/bash
#SBATCH -J MA_Pkg3_discovery_donor
#SBATCH -p CU
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=12
#SBATCH --mem=140G
#SBATCH -t 24:00:00
#SBATCH -o /gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/slurm_logs/%x-%j.out
#SBATCH -e /gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/slurm_logs/%x-%j.err

set -euo pipefail

R_SCRIPT="/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/scripts/Package3_discovery_donor_aware_disease_association.R"
INPUT_RDS="/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Package2_strict_input_fixed_from_original.rds"
OUTDIR="/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/Package3_Discovery_DonorAwareDiseaseAssociation"

mkdir -p "${OUTDIR}"


export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export VECLIB_MAXIMUM_THREADS=${SLURM_CPUS_PER_TASK}
export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK}

Rscript "${R_SCRIPT}" \
  --rds "${INPUT_RDS}" \
  --outdir "${OUTDIR}" \
  --assay RNA \
  --cluster_col AUTO \
  --donor_col AUTO \
  --dx_col AUTO \
  --case_label ASD \
  --control_label Control \
  --min_cells_donor 20 \
  --min_cells_cluster 5 \
  --threads ${SLURM_CPUS_PER_TASK}
