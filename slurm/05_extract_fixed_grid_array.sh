#!/bin/bash
#SBATCH --job-name=fixed_grid_extract
#SBATCH --partition=med
#SBATCH --array=1-100
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --output=logs/fixed_grid_extract_%A_%a.log
#SBATCH --error=logs/fixed_grid_extract_%A_%a.log

# Usage: sbatch --export=ALL,COND=E slurm/05_extract_fixed_grid_array.sh   (COND = I, H, F or E)

set -euo pipefail

export PROJECT_ROOT="${PROJECT_ROOT:-$SLURM_SUBMIT_DIR}"   # submit from the project root
export FIXED_GRID_BATCH="single_node_cond${COND:?set COND to I, H, F or E}_ds_cr_seeds001_100"
export TMPDIR="${TMPDIR:-/tmp}/fixed_grid_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
export TMP="$TMPDIR"
export TEMP="$TMPDIR"

mkdir -p "$PROJECT_ROOT/logs" "$TMPDIR"
cd "$PROJECT_ROOT"

module purge
module load R/4.4.2-gcc141-base

Rscript --vanilla code/05_aggregate/05_extract_fixed_grid.R "${SLURM_ARRAY_TASK_ID}"
