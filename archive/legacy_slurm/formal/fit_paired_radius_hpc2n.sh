#!/bin/bash
#SBATCH --account=hpc2n2025-235
#SBATCH --job-name=paired_radius_fit
#SBATCH --array=0-19
#SBATCH --partition=cpu_sky
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=7-00:00:00
#SBATCH --output=logs/hpc2n_paired_radius_%A_%a.log
#SBATCH --error=logs/hpc2n_paired_radius_%A_%a.log
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=fjhmrqx@126.com

set -euo pipefail

module --force purge >/dev/null 2>&1
module load GCC/13.2.0 R/4.4.1

export PROJECT_ROOT=/proj/nobackup/hpc2n2025-235/jiahui/bgam_eeg_project
export R_LIBS_USER=/proj/nobackup/hpc2n2025-235/jiahui/Rlibs/4.4
export TMPDIR="$PROJECT_ROOT/tmp/paired_radius_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
export TMP="$TMPDIR"
export TEMP="$TMPDIR"
export CMDSTAN_CACHE_DIR="$PROJECT_ROOT/cmdstan_cache/paired_radius_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"

mkdir -p "$PROJECT_ROOT/logs" "$TMPDIR" "$CMDSTAN_CACHE_DIR"

RADIUS_RULES=(half-max)
RADIUS_INDEX=$((SLURM_ARRAY_TASK_ID / 20))
REPLICATE_SEED=$((SLURM_ARRAY_TASK_ID % 20 + 1))
RADIUS_RULE="${RADIUS_RULES[$RADIUS_INDEX]}"

cd "$PROJECT_ROOT"
Rscript code/07_fit_paired_radius_bgam.R "$RADIUS_RULE" "$REPLICATE_SEED"
