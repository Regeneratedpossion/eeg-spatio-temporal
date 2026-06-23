#!/bin/bash
#SBATCH --job-name=bgam_B_ps_sig095
#SBATCH --partition=ultralong
#SBATCH --array=1-2
#SBATCH --time=10-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH -o logs/bgam_B_ps_sig095_%A_%a.log
#SBATCH -e logs/bgam_B_ps_sig095_%A_%a.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=fjhmrqx@126.com

mkdir -p logs
module load R/4.4.2-gcc141-base

export PROJECT_ROOT=/work/smjifan1/bgam_eeg_project
export SCREEN_SEED=1
export SCREEN_ITER=3000
export SCREEN_WARMUP=1000
export SCREEN_ADAPT_DELTA=0.95
export SCREEN_MAX_TREEDEPTH=10
export SCREEN_OUT_SUBDIR=basis_screening_B_ps_sigma_true_adapt095_scaling

# The R script uses SLURM_ARRAY_TASK_ID=2 for B_k1212_ps.
export INNER_TASK_ID=2

case "$SLURM_ARRAY_TASK_ID" in
  1)
    export SCREEN_DATA_SUBDIR=montecarlo_p5_gaussian_hardlabel
    ;;
  2)
    export SCREEN_DATA_SUBDIR=montecarlo_p20_gaussian_hardlabel
    ;;
esac

cd "$PROJECT_ROOT"
echo "data=$SCREEN_DATA_SUBDIR inner_task=$INNER_TASK_ID adapt=$SCREEN_ADAPT_DELTA td=$SCREEN_MAX_TREEDEPTH"
SLURM_ARRAY_TASK_ID="$INNER_TASK_ID" Rscript code/18_basis_screening_full_ab_sigma_true_control.R
