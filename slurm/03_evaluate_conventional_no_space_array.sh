#!/bin/bash
#SBATCH --job-name=evaluate_conv_nospace
#SBATCH --partition=short
#SBATCH --array=1-100%20
#SBATCH --time=00:45:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=logs/evaluate_conv_nospace_%A_%a.log

# Evaluate the conventional methods on one dataset per array task with
# code/03_conventional/03_evaluate_conventional.R.
# Usage: sbatch --export=ALL,COND=E slurm/03_evaluate_conventional_no_space_array.sh   (COND = I, H, F or E)

set -euo pipefail

export PROJECT_ROOT="${PROJECT_ROOT:-$SLURM_SUBMIT_DIR}"   # submit from the project root
cd "$PROJECT_ROOT"

module purge
module load R/4.4.2-gcc141-base

RESULTS_BATCH="single_node_cond${COND:?set COND to I, H, F or E}_ds_cr_seeds001_100"
SEED_DIR="results/${RESULTS_BATCH}/$(printf 'seed_%03d' "$SLURM_ARRAY_TASK_ID")"
DATA_SUBDIR=$(Rscript --vanilla -e "cat(readRDS('${SEED_DIR}/config.rds')\$data_subdir)")
SEED_TAG=$(printf '%03d' "$SLURM_ARRAY_TASK_ID")

Rscript --vanilla code/03_conventional/03_evaluate_conventional_no_space.R \
    "data/${DATA_SUBDIR}/sim_${SEED_TAG}.rds" \
    "data/${DATA_SUBDIR}/gt_${SEED_TAG}.rds" \
    "${SEED_DIR}/evaluation_no_space"
