#!/bin/bash
#SBATCH --job-name=ppc
#SBATCH --partition=short
#SBATCH --array=1-362
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --output=logs/ppc_%A_%a.log

# Posterior predictive checks of one passed fit per array task (code/08_ppc/ppc_one_dataset.R).
# Task n takes line n of results/passed_fits.txt (scenario and seed). Peak use is about 5 GB, so 6G is requested; 5G is too little.
# Usage: sbatch slurm/08_ppc_array.sh   (submit from the project root)

set -euo pipefail

export PROJECT_ROOT="${PROJECT_ROOT:-$SLURM_SUBMIT_DIR}"
cd "$PROJECT_ROOT"

source "$PROJECT_ROOT/slurm/env.sh"

read -r SCENARIO SEED < <(sed -n "${SLURM_ARRAY_TASK_ID}p" results/passed_fits.txt)
echo "scenario ${SCENARIO} seed ${SEED}"

Rscript --vanilla code/08_ppc/ppc_one_dataset.R "$SCENARIO" "$SEED" 200
