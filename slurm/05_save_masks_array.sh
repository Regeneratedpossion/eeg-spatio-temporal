#!/bin/bash
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/save_masks_%x_%A_%a.log

# One mask file per array task. Submit with the memory and time of the method:
#   sbatch --job-name=masks_bgam_I --mem=10G --time=00:15:00 --array=1-100%100 \
#     --export=ALL,COND=1,WHAT=bgam slurm/05_save_masks_array.sh
set -euo pipefail
export PROJECT_ROOT="${PROJECT_ROOT:-$SLURM_SUBMIT_DIR}"   # submit from the project root
mkdir -p "$PROJECT_ROOT/logs"
cd "$PROJECT_ROOT"
module purge
module load R/4.4.2-gcc141-base
Rscript --vanilla code/05_aggregate/05_save_masks.R "$COND" "$SLURM_ARRAY_TASK_ID" "$WHAT"
