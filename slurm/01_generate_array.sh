#!/bin/bash
# Usage: sbatch --export=ALL,COND=4 --job-name=gen_cond4 --output=logs/gen_cond4_%A_%a.log slurm/01_generate_array.sh
#SBATCH --job-name=gen_cond
#SBATCH --partition=ultralong
#SBATCH --array=1-100
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --output=logs/gen_cond_%A_%a.log
#SBATCH --error=logs/gen_cond_%A_%a.log

set -euo pipefail

export PROJECT_ROOT="${PROJECT_ROOT:-$SLURM_SUBMIT_DIR}"   # submit from the project root
export TMPDIR="${SLURM_TMPDIR:-/tmp}/gen_${COND}_${SLURM_ARRAY_TASK_ID}_${SLURM_JOB_ID}"
mkdir -p "$TMPDIR" "$PROJECT_ROOT/logs"

module purge
module load R/4.4.2-gcc141-base
cd "$PROJECT_ROOT"

Rscript --vanilla code/01_simulate/01_generate_dataset.R
