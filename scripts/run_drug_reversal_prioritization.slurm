#!/bin/bash
#SBATCH -J MA_Pkg9_DrugRev
#SBATCH -p CU
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=16G
#SBATCH -t 08:00:00
#SBATCH -o /gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/slurm_logs/%x_%j.out
#SBATCH -e /gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/slurm_logs/%x_%j.err

set -euo pipefail


ROOT=/gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision
OUTDIR=${ROOT}/Package9_PerturbationPrioritization

mkdir -p ${ROOT}/slurm_logs
mkdir -p ${OUTDIR}

python3 drug_reversal_prioritization.py \
  --signature_tsv ${ROOT}/Package4_CrossCohortValidation_Velmeshev_v2/08_discovery_frozen_state_signature.tsv \
  --pooled_de_tsv ${ROOT}/Package3_Discovery_DonorAwareDiseaseAssociation/13_pseudobulk_all_retained_ASD_vs_Control.tsv.gz \
  --outdir ${OUTDIR} \
  --top_primary 100 \
  --top_sensitivity 150
