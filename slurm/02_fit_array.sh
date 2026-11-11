#!/bin/bash
# Usage: sbatch --job-name=ds_cr_cond4_fit --output=logs/ds_cr_cond4_fit_%A_%a.log --export=ALL,COND=4 slurm/02_fit_array.sh
#SBATCH --job-name=ds_cr_fit
#SBATCH --partition=ultralong
#SBATCH --array=1-100%50
#SBATCH --constraint=xeon_e52640v4
#SBATCH --time=4-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=logs/ds_cr_fit_%A_%a.log
#SBATCH --error=logs/ds_cr_fit_%A_%a.log

set -euo pipefail

export PROJECT_ROOT="${PROJECT_ROOT:-$SLURM_SUBMIT_DIR}"   # submit from the project root
export CMDSTAN_PATH="${CMDSTAN_PATH:-}"
export TMPDIR="${SLURM_TMPDIR:-/tmp}/ds_cr_seed_${SLURM_ARRAY_TASK_ID}_${SLURM_JOB_ID}"
unset STAN_NUM_THREADS
mkdir -p "$TMPDIR" "$PROJECT_ROOT/logs"

source "$PROJECT_ROOT/slurm/env.sh"
cd "$PROJECT_ROOT"

echo "job_id=$SLURM_JOB_ID seed=$SLURM_ARRAY_TASK_ID basis=ds_cr cond=$COND started=$(date --iso-8601=seconds)"
bash slurm/02_fit_with_retry.sh
echo "finished=$(date --iso-8601=seconds)"
