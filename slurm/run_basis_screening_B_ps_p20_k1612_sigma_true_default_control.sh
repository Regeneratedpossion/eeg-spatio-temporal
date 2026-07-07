#!/bin/bash
#SBATCH --job-name=bgam_Bps_p20_k1612
#SBATCH --partition=ultralong
#SBATCH --time=10-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH -o logs/bgam_Bps_p20_k1612_%j.log
#SBATCH -e logs/bgam_Bps_p20_k1612_%j.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=fjhmrqx@126.com

set -euo pipefail

mkdir -p logs /work/smjifan1/tmp
module load R/4.4.2-gcc141-base

export PROJECT_ROOT=/work/smjifan1/bgam_eeg_project
export TMPDIR=/work/smjifan1/tmp
export TMP=/work/smjifan1/tmp
export TEMP=/work/smjifan1/tmp

export SCREEN_DATA_SUBDIR=montecarlo_p20_gaussian_hardlabel
export SCREEN_SEED=1
export SCREEN_ITER=3000
export SCREEN_WARMUP=1000
export SCREEN_K_XY=16
export SCREEN_K_T=12
export SCREEN_ADAPT_DELTA=0.80
export SCREEN_MAX_TREEDEPTH=10
export SCREEN_OUT_SUBDIR=basis_screening_B_ps_p20_k1612_sigma_true_default_control

# code/18_basis_screening_full_ab_sigma_true_control.R uses task 2 for B_ps.
export INNER_TASK_ID=2

cd "$PROJECT_ROOT"
echo "data=$SCREEN_DATA_SUBDIR model=B_ps k=($SCREEN_K_XY,$SCREEN_K_T) adapt=$SCREEN_ADAPT_DELTA td=$SCREEN_MAX_TREEDEPTH out=$SCREEN_OUT_SUBDIR TMPDIR=$TMPDIR"
SLURM_ARRAY_TASK_ID="$INNER_TASK_ID" Rscript code/18_basis_screening_full_ab_sigma_true_control.R
