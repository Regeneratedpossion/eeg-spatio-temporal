#!/bin/bash
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/old_rule_%x_%A_%a.log

# Re-run the old onset/offset rule for one seed per array task.
# Old rule: component of at least 2 cells, then at least 3 electrodes for at least 2 consecutive time points.
# The old scripts live in /work/smjifan1/old_rule_run_20261013 so the current code is not touched.
# Output goes to <seed>/evaluation_old_rule_20261009 and nothing else is overwritten.
# Usage (low memory, high concurrency):
#   sbatch --job-name=oldrule_bgam_F --mem=6G --time=00:20:00 --array=1-100%100 \
#     --export=ALL,COND=F,WHAT=bgam slurm/60_rerun_old_rule_array.sh
#   sbatch --job-name=oldrule_conv_F --mem=3G --time=00:30:00 --array=1-100%100 \
#     --export=ALL,COND=F,WHAT=conv slurm/60_rerun_old_rule_array.sh
set -euo pipefail

export PROJECT_ROOT=/work/smjifan1/old_rule_run_20261013
export R_LIBS_USER=/work/smjifan1/R/4.4.2:/home/smjifan1/R/x86_64-pc-linux-gnu-library/4.4
cd "$PROJECT_ROOT"

module purge
module load R/4.4.2-gcc141-base

TAG=$(printf '%03d' "$SLURM_ARRAY_TASK_ID")
SEED_DIR="results/single_node_cond${COND}_ds_cr_seeds001_100/seed_${TAG}"
DATA_DIR="data/single_node_cond${COND}_seeds001_100"
OUT_DIR="${OUT_DIR:-${SEED_DIR}/evaluation_old_rule_20261009}"

test -d "$SEED_DIR" && test -f "${DATA_DIR}/sim_${TAG}.rds"
mkdir -p "$OUT_DIR"

if [ "$WHAT" = "bgam" ]; then
    [ -f "${OUT_DIR}/evaluation_main.csv" ] && exit 0
    Rscript --vanilla code/23_evaluate_bgam_fit.R \
        "${SEED_DIR}/fit.rds" "${DATA_DIR}/sim_${TAG}.rds" "${DATA_DIR}/gt_${TAG}.rds" "$OUT_DIR"
else
    [ -f "${OUT_DIR}/evaluation_conventional.csv" ] && exit 0
    Rscript --vanilla code/25_evaluate_conventional.R \
        "${DATA_DIR}/sim_${TAG}.rds" "${DATA_DIR}/gt_${TAG}.rds" "$OUT_DIR"
fi
