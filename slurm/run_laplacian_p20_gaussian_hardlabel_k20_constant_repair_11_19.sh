#!/bin/bash
#SBATCH --job-name=lap_p20gh_const_rep
#SBATCH --partition=med
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH -o logs/lap_p20gh_const_rep_%j.log
#SBATCH -e logs/lap_p20gh_const_rep_%j.log

mkdir -p logs
module load R/4.4.2-gcc141-base

export PROJECT_ROOT=/work/smjifan1/bgam_eeg_project
export LAPLACE_DATA_SUBDIR=montecarlo_p20_gaussian_hardlabel
export LAPLACE_SEED=1
export LAPLACE_K=20
export LAPLACE_INCLUDE_CONSTANT=true
export LAPLACE_TASKS="c(11,19)"
export LAPLACE_TIME_K=10
export LAPLACE_ITER=3000
export LAPLACE_WARMUP=1000
export LAPLACE_ADAPT_DELTA=0.9999
export LAPLACE_MAX_TREEDEPTH=15
export LAPLACE_OUT_SUBDIR=laplacian_basis_brms_p20_gaussian_hardlabel_constant_repair_11_19/k20

cd "$PROJECT_ROOT"
Rscript code/12_laplacian_basis_brms_fixed.R
