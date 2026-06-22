#!/bin/bash
#SBATCH --job-name=lap_basis_k10_fix
#SBATCH --partition=ultralong
#SBATCH --time=10-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH -o logs/lap_basis_k10_fix_%j.log
#SBATCH -e logs/lap_basis_k10_fix_%j.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=fjhmrqx@126.com

mkdir -p logs
module load R/4.4.2-gcc141-base

export PROJECT_ROOT=/work/smjifan1/bgam_eeg_project
export LAPLACE_DATA_SUBDIR=montecarlo_p5_gaussian_hardlabel
export LAPLACE_SEED=1
export LAPLACE_K=10
export LAPLACE_TASKS="1:10"
export LAPLACE_TIME_K=10
export LAPLACE_ITER=3000
export LAPLACE_WARMUP=1000
export LAPLACE_ADAPT_DELTA=0.999
export LAPLACE_MAX_TREEDEPTH=13

cd "$PROJECT_ROOT"
Rscript code/12_laplacian_basis_brms_fixed.R
