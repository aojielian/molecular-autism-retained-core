#!/bin/bash
#SBATCH -J Build_PE_GlobalRDS
#SBATCH -p CU
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH -t 24:00:00
#SBATCH -o /gpfs/hpc/home/lijc/lianaoj/autism_scRNA/PsychENCODE_Science_2024/slurm_logs/%x.%j.out
#SBATCH -e /gpfs/hpc/home/lijc/lianaoj/autism_scRNA/PsychENCODE_Science_2024/slurm_logs/%x.%j.err

set -euo pipefail

LOGDIR=/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/PsychENCODE_Science_2024/slurm_logs
SCRIPT=/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/scripts/build_psychencode_global_object.R

mkdir -p "${LOGDIR}"

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export VECLIB_MAXIMUM_THREADS=${SLURM_CPUS_PER_TASK}
export NUMEXPR_NUM_THREADS=${SLURM_CPUS_PER_TASK}


Rscript "${SCRIPT}"
