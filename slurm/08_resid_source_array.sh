#!/bin/bash
#SBATCH --job-name=resid_source
#SBATCH --partition=short
#SBATCH --array=1-362
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --output=logs/resid_source_%A_%a.log

# Residual-source diagnostic of one passed fit per array task (code/08_ppc/check_residual_source.R).
# Task n takes line n of results/passed_fits.txt (scenario and seed). 5G is enough.
# Usage: sbatch slurm/08_resid_source_array.sh   (submit from the project root)

set -euo pipefail

export PROJECT_ROOT="${PROJECT_ROOT:-$SLURM_SUBMIT_DIR}"
cd "$PROJECT_ROOT"

source "$PROJECT_ROOT/slurm/env.sh"

read -r SCENARIO SEED < <(sed -n "${SLURM_ARRAY_TASK_ID}p" results/passed_fits.txt)
echo "scenario ${SCENARIO} seed ${SEED}"

Rscript --vanilla code/08_ppc/check_residual_source.R "$SCENARIO" "$SEED" 200
