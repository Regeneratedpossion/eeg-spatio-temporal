#!/bin/bash
#SBATCH --job-name=chain_cov
#SBATCH --partition=short
#SBATCH --array=1-362%150
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --output=logs/chain_cov_%A_%a.log
# One passed fit per array task; the list is resources/passed_fits.txt (condition, seed).
set -euo pipefail
export PROJECT_ROOT="${PROJECT_ROOT:-$SLURM_SUBMIT_DIR}"   # submit from the project root
cd "$PROJECT_ROOT"; module purge; module load R/4.4.2-gcc141-base
read -r COND SEED < <(sed -n "${SLURM_ARRAY_TASK_ID}p" resources/passed_fits.txt)
Rscript --vanilla code/05_aggregate/05_per_chain_coverage.R "$COND" "$SEED"
