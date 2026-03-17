#!/bin/bash
#SBATCH -J Package4b_direction_consistency_v2
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
SCRIPT=${PROJECT_DIR}/scripts/Package4b_direction_consistency_manual_v2.R


Rscript "${SCRIPT}"
