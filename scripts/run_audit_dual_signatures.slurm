#!/bin/bash
#SBATCH -J auditDualSig
#SBATCH -p CU
#SBATCH -N 1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH -t 4:00:00
#SBATCH -o /gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/slurm_logs/%x.%j.out
#SBATCH -e /gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/slurm_logs/%x.%j.err

set -euo pipefail

BASE=/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision
SCRIPT=${BASE}/scripts/audit_dual_signatures.R

mkdir -p ${BASE}/slurm_logs
mkdir -p ${BASE}/Signature_Audit
mkdir -p ${BASE}/scripts


Rscript ${SCRIPT}
