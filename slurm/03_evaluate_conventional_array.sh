#!/bin/bash
#SBATCH --job-name=evaluate_conv
#SBATCH --partition=short
#SBATCH --array=1-100%50
#SBATCH --time=00:45:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --output=logs/evaluate_conv_%A_%a.log

# Evaluate the conventional methods on one dataset per array task with
# code/03_conventional/03_evaluate_conventional.R.
# Usage: sbatch --export=ALL,COND=4 slurm/03_evaluate_conventional_array.sh   (COND = 1, 2, 3 or 4)

set -euo pipefail

export PROJECT_ROOT="${PROJECT_ROOT:-$SLURM_SUBMIT_DIR}"   # submit from the project root
cd "$PROJECT_ROOT"

source "$PROJECT_ROOT/slurm/env.sh"

RESULTS_BATCH="scenario${COND:?set COND to 1, 2, 3 or 4}"
SEED_DIR="results/${RESULTS_BATCH}/$(printf 'seed_%03d' "$SLURM_ARRAY_TASK_ID")"
DATA_SUBDIR="scenario${COND}"
SEED_TAG=$(printf '%03d' "$SLURM_ARRAY_TASK_ID")

Rscript --vanilla code/03_conventional/03_evaluate_conventional.R \
    "data/${DATA_SUBDIR}/sim_${SEED_TAG}.rds" \
    "data/${DATA_SUBDIR}/gt_${SEED_TAG}.rds" \
    "${SEED_DIR}/evaluation"
