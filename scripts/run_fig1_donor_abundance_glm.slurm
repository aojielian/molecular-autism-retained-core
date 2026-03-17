#!/bin/bash
#SBATCH -J MA_Fig1_GLM
#SBATCH -p CU
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -t 8:00:00
#SBATCH -o /gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/slurm_logs/%x.%j.out
#SBATCH -e /gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/slurm_logs/%x.%j.err

set -euo pipefail

BASE=/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision
SCRIPT=${BASE}/scripts/step_fig1_donor_abundance_glm.R

mkdir -p ${BASE}/scripts
mkdir -p ${BASE}/slurm_logs
mkdir -p ${BASE}/tmp

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export TMPDIR=${BASE}/tmp/${SLURM_JOB_ID}
mkdir -p ${TMPDIR}


Rscript ${SCRIPT}
