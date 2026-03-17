#!/bin/bash
#SBATCH -J GandalPkg8
#SBATCH -p CU
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=48G
#SBATCH -t 12:00:00
#SBATCH -o /gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/slurm_logs/%x_%j.out
#SBATCH -e /gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/slurm_logs/%x_%j.err

set -euo pipefail


ROOT=/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision
OUTDIR=${ROOT}/Package8_ExternalValidation_Gandal2022
FROZENDIR=${OUTDIR}/FrozenGeneSets
SIG_TSV=${ROOT}/Package4_CrossCohortValidation_Velmeshev_v2/08_discovery_frozen_state_signature.tsv
RDATA=/gpfs/hpc/home/lijc/lianaoj/AAAY_data/Gandal_2022/Gene_NormalizedExpression_Metadata_wModelMatrix.RData
SCRIPT=${ROOT}/scripts/Step08_Gandal2022_bulk_program_validation_v3.R

mkdir -p ${ROOT}/slurm_logs
mkdir -p ${OUTDIR}
mkdir -p ${FROZENDIR}

echo "[INFO] Export frozen candidate/reference gene sets from ${SIG_TSV}"

awk -F'\t' 'NR>1 && $1=="disc_cluster2_up"{print $3}' ${SIG_TSV} > ${FROZENDIR}/frozen_candidate_genes_top100.txt
awk -F'\t' 'NR>1 && $1=="disc_cluster0_up"{print $3}' ${SIG_TSV} > ${FROZENDIR}/frozen_reference_genes_top100.txt

echo "[INFO] Candidate gene count:"
wc -l ${FROZENDIR}/frozen_candidate_genes_top100.txt

echo "[INFO] Reference gene count:"
wc -l ${FROZENDIR}/frozen_reference_genes_top100.txt

echo "[INFO] Run Gandal bulk program-level validation"
Rscript ${SCRIPT} \
  ${RDATA} \
  ${FROZENDIR}/frozen_candidate_genes_top100.txt \
  ${FROZENDIR}/frozen_reference_genes_top100.txt \
  ${OUTDIR}

echo "[INFO] Finished."
