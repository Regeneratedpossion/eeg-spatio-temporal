#!/bin/bash
#SBATCH --job-name=lap_p20gh_const_cont
#SBATCH --partition=med
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH -o logs/lap_p20gh_const_cont_%j.log
#SBATCH -e logs/lap_p20gh_const_cont_%j.log

mkdir -p logs
module load R/4.4.2-gcc141-base

export PROJECT_ROOT=/work/smjifan1/bgam_eeg_project
export LAPLACE_DATA_SUBDIR=montecarlo_p20_gaussian_hardlabel
export LAPLACE_SEED=1
export LAPLACE_K=20
export LAPLACE_TIME_K=10
export LAPLACE_ITER=3000
export LAPLACE_WARMUP=1000
export LAPLACE_ADAPT_TAG=0p999
export LAPLACE_TAG_SUFFIX=constant
export LAPLACE_FIT_SUBDIR=laplacian_basis_brms_p20_gaussian_hardlabel_constant_long/k20

cd "$PROJECT_ROOT"
Rscript code/14_laplacian_continuous_recovery.R
