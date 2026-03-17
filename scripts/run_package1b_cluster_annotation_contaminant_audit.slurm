#!/bin/bash
#SBATCH -J MA_Pkg1b_Audit
#SBATCH -p CU
#SBATCH -N 1
#SBATCH --cpus-per-task=12
#SBATCH --mem=120G
#SBATCH -t 24:00:00
#SBATCH -o /gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/slurm_logs/%x.%j.out
#SBATCH -e /gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/slurm_logs/%x.%j.err

set -euo pipefail

LOGDIR=/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/slurm_logs
SCRIPT=/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/scripts/package1b_cluster_annotation_contaminant_audit.R

mkdir -p "${LOGDIR}"
mkdir -p /gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/scripts

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export VECLIB_MAXIMUM_THREADS=${SLURM_CPUS_PER_TASK}
export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK}


Rscript "${SCRIPT}"
