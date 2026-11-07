#!/bin/bash
#SBATCH --job-name=fixed_grid_aggregate
#SBATCH --partition=med
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=logs/fixed_grid_aggregate_%j.log
#SBATCH --error=logs/fixed_grid_aggregate_%j.log

set -euo pipefail

export PROJECT_ROOT="${PROJECT_ROOT:-$SLURM_SUBMIT_DIR}"   # submit from the project root
cd "$PROJECT_ROOT"

module purge
module load R/4.4.2-gcc141-base

Rscript --vanilla code/41_aggregate_and_plot_fixed_grid.R
