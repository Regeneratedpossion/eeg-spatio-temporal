#!/bin/bash
#SBATCH --job-name=cov_z
#SBATCH --partition=short
#SBATCH --array=1-362%150
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --output=logs/cov_z_%A_%a.log
# Posterior mean, SD and 90% interval of g for one passed fit per array task (code/extra/04_coverage_z.R).
# Task n takes line n of results/passed_fits.txt (scenario and seed).
# Usage: sbatch slurm/04_coverage_z_array.sh   (submit from the project root)
set -euo pipefail
export PROJECT_ROOT="${PROJECT_ROOT:-$SLURM_SUBMIT_DIR}"
cd "$PROJECT_ROOT"
source "$PROJECT_ROOT/slurm/env.sh"
read -r SCENARIO SEED < <(sed -n "${SLURM_ARRAY_TASK_ID}p" results/passed_fits.txt)
Rscript --vanilla code/extra/04_coverage_z.R "$SCENARIO" "$SEED"
