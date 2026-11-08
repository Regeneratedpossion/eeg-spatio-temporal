#!/bin/bash
#SBATCH --job-name=evaluate_bgam
#SBATCH --partition=short
#SBATCH --array=1-100%20
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=logs/evaluate_bgam_%A_%a.log

# Evaluate one fit per array task with code/04_evaluate/04_evaluate_bgam.R.
# Usage: sbatch --export=ALL,RESULTS_BATCH=<results subdirectory> slurm/04_evaluate_bgam_array.sh

set -euo pipefail

export PROJECT_ROOT="${PROJECT_ROOT:-$SLURM_SUBMIT_DIR}"   # submit from the project root
cd "$PROJECT_ROOT"

module purge
module load R/4.4.2-gcc141-base

SEED_DIR="results/${RESULTS_BATCH}/$(printf 'seed_%03d' "$SLURM_ARRAY_TASK_ID")"
DATA_SUBDIR=$(Rscript --vanilla -e "cat(readRDS('${SEED_DIR}/config.rds')\$data_subdir)")
SEED_TAG=$(printf '%03d' "$SLURM_ARRAY_TASK_ID")

Rscript --vanilla code/04_evaluate/04_evaluate_bgam.R \
    "${SEED_DIR}/fit.rds" \
    "data/${DATA_SUBDIR}/sim_${SEED_TAG}.rds" \
    "data/${DATA_SUBDIR}/gt_${SEED_TAG}.rds" \
    "${SEED_DIR}/evaluation"
