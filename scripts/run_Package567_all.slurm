#!/bin/bash
#SBATCH -J Package567_runall
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
SCRIPT_DIR=${PROJECT_DIR}/scripts


echo "[Start] Package 5"
Rscript "${SCRIPT_DIR}/Package5_directional_concordance_global_summary.R"

echo "[Start] Package 6"
Rscript "${SCRIPT_DIR}/Package6_validation_results_freezeout.R"

echo "[Start] Package 7"
Rscript "${SCRIPT_DIR}/Package7_minimal_robustness_analysis.R"

echo "[Done] Package 5/6/7 all finished."
