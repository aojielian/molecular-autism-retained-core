#!/bin/bash
#SBATCH -J GandalPkg8_LOO
#SBATCH -p CU
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=32G
#SBATCH -t 08:00:00
#SBATCH -o /gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/slurm_logs/%x_%j.out
#SBATCH -e /gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/slurm_logs/%x_%j.err

set -euo pipefail


ROOT=/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision
INDIR=${ROOT}/Package8_ExternalValidation_Gandal2022
INFILE=${INDIR}/07_sample_level_program_scores.tsv
OUTDIR=${ROOT}/Package8b_Gandal2022_LOO_Sensitivity
SCRIPT=${ROOT}/scripts/Step08b_Gandal2022_leave_one_region_out_sensitivity.R

mkdir -p ${ROOT}/slurm_logs
mkdir -p ${OUTDIR}

echo "[INFO] Input file: ${INFILE}"
echo "[INFO] Output dir: ${OUTDIR}"

Rscript ${SCRIPT} ${INFILE} ${OUTDIR}

echo "[INFO] Finished."
