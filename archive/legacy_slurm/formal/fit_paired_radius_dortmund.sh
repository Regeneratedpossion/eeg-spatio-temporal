#!/bin/bash
#SBATCH --job-name=paired_radius_fit
#SBATCH --partition=ultralong
#SBATCH --time=10-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --array=0-19
#SBATCH -o logs/paired_radius_fit_%A_%a.log
#SBATCH -e logs/paired_radius_fit_%A_%a.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=fjhmrqx@126.com

set -euo pipefail

module load R/4.4.2-gcc141-base

export PROJECT_ROOT=/work/smjifan1/bgam_eeg_project
export R_LIBS_USER=/work/smjifan1/R/4.4.2:/home/smjifan1/R/x86_64-pc-linux-gnu-library/4.4
export TMPDIR=/work/smjifan1/tmp/paired_radius_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
export TMP=$TMPDIR
export TEMP=$TMPDIR
export CMDSTAN_CACHE_DIR=$PROJECT_ROOT/cmdstan_cache/paired_radius_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}
mkdir -p "$PROJECT_ROOT/logs" "$TMPDIR" "$CMDSTAN_CACHE_DIR"

RADIUS_RULES=(half-max)
RADIUS_INDEX=$((SLURM_ARRAY_TASK_ID / 20))
REPLICATE_SEED=$((SLURM_ARRAY_TASK_ID % 20 + 1))
RADIUS_RULE="${RADIUS_RULES[$RADIUS_INDEX]}"

cd "$PROJECT_ROOT"
Rscript code/07_fit_paired_radius_bgam.R "$RADIUS_RULE" "$REPLICATE_SEED"
