#!/bin/bash
#SBATCH --job-name=evaluate_bgam
#SBATCH --partition=short
#SBATCH --array=1-100%50
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --output=logs/evaluate_bgam_%A_%a.log

# Evaluate one fit per array task with code/04_evaluate/04_evaluate_bgam.R.
# Usage: sbatch --export=ALL,COND=4 slurm/04_evaluate_bgam_array.sh   (COND = 1, 2, 3 or 4)

set -euo pipefail

export PROJECT_ROOT="${PROJECT_ROOT:-$SLURM_SUBMIT_DIR}"   # submit from the project root
cd "$PROJECT_ROOT"

source "$PROJECT_ROOT/slurm/env.sh"

RESULTS_BATCH="scenario${COND:?set COND to 1, 2, 3 or 4}"
SEED_DIR="results/${RESULTS_BATCH}/$(printf 'seed_%03d' "$SLURM_ARRAY_TASK_ID")"
DATA_SUBDIR=$(Rscript --vanilla -e "cat(readRDS('${SEED_DIR}/config.rds')\$data_subdir)")
SEED_TAG=$(printf '%03d' "$SLURM_ARRAY_TASK_ID")

Rscript --vanilla code/04_evaluate/04_evaluate_bgam.R \
    "${SEED_DIR}/fit.rds" \
    "data/${DATA_SUBDIR}/sim_${SEED_TAG}.rds" \
    "data/${DATA_SUBDIR}/gt_${SEED_TAG}.rds" \
    "${SEED_DIR}/evaluation"
