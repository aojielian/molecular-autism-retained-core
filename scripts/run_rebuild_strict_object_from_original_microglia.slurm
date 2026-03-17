#!/bin/bash
#SBATCH -J Fix_StrictObj
#SBATCH -p CU
#SBATCH -N 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH -t 08:00:00
#SBATCH -o /gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/slurm_logs/%x.%j.out
#SBATCH -e /gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/slurm_logs/%x.%j.err

set -euo pipefail


Rscript /gpfs/hpc/home/lijc/lianaoj/autism_scRNA/MolecularAutism_Revision/scripts/rebuild_strict_object_from_original_microglia.R
