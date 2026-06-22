#!/bin/bash
#SBATCH --job-name=lap_basis_ladder
#SBATCH --partition=med
#SBATCH --array=1-4
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH -o logs/lap_basis_ladder_%A_%a.log
#SBATCH -e logs/lap_basis_ladder_%A_%a.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=fjhmrqx@126.com

mkdir -p logs
module load R/4.4.2-gcc141-base

export PROJECT_ROOT=/work/smjifan1/bgam_eeg_project
export LAPLACE_DATA_SUBDIR=montecarlo_p5_gaussian_hardlabel
export LAPLACE_SEED=1
export LAPLACE_TIME_K=10
export LAPLACE_ITER=1000
export LAPLACE_WARMUP=500
export LAPLACE_ADAPT_DELTA=0.999
export LAPLACE_MAX_TREEDEPTH=13

K_VALUES=(10 15 20 30)
K="${K_VALUES[$((SLURM_ARRAY_TASK_ID - 1))]}"

export LAPLACE_K="$K"
export LAPLACE_TASKS="1:${K}"
export LAPLACE_OUT_SUBDIR="laplacian_basis_brms_ladder/k${K}"

cd "$PROJECT_ROOT"
Rscript code/12_laplacian_basis_brms_fixed.R
