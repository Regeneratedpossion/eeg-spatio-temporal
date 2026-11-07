#!/bin/bash
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/rescorr_%x_%A_%a.log

# Residual correlation of one BGAM fit, one seed per array task.
# Usage (low memory, high concurrency):
#   sbatch --job-name=rc_I --mem=8G --time=00:20:00 --array=1-10 \
#     --export=ALL,COND=I slurm/61_coverage_array.sh
set -euo pipefail

export PROJECT_ROOT=/work/smjifan1/bgam_eeg_project
export R_LIBS_USER=/work/smjifan1/R/4.4.2:/home/smjifan1/R/x86_64-pc-linux-gnu-library/4.4
cd "$PROJECT_ROOT"

module purge
module load R/4.4.2-gcc141-base

TAG=$(printf '%03d' "$SLURM_ARRAY_TASK_ID")
SEED_DIR="results/single_node_cond${COND}_ds_cr_seeds001_100/seed_${TAG}"
DATA_DIR="data/single_node_cond${COND}_seeds001_100"
OUT_DIR="${SEED_DIR}/rescorr_20261015"

test -d "$SEED_DIR" && test -f "${DATA_DIR}/sim_${TAG}.rds"
[ -f "${OUT_DIR}/residual_autocorr.csv" ] && exit 0
Rscript --vanilla code/27_residual_autocorr.R \
    "${SEED_DIR}/fit.rds" "${DATA_DIR}/sim_${TAG}.rds" "$OUT_DIR"
