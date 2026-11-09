#!/bin/bash
#SBATCH --partition=short
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --job-name=agg_passed
#SBATCH --array=0-3
#SBATCH --output=logs/agg_passed_%A_%a.log
set -euo pipefail
export PROJECT_ROOT="${PROJECT_ROOT:-$SLURM_SUBMIT_DIR}"   # submit from the project root
cd "$PROJECT_ROOT"; module purge; module load R/4.4.2-gcc141-base
C=(1 3 2 4); export COND=${C[$SLURM_ARRAY_TASK_ID]}
Rscript --vanilla code/05_aggregate/05_aggregate_fixed_grid.R
