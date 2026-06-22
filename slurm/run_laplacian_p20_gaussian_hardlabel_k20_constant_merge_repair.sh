#!/bin/bash
#SBATCH --job-name=lap_p20gh_merge_rep
#SBATCH --partition=med
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH -o logs/lap_p20gh_merge_rep_%j.log
#SBATCH -e logs/lap_p20gh_merge_rep_%j.log

mkdir -p logs
module load R/4.4.2-gcc141-base

export PROJECT_ROOT=/work/smjifan1/bgam_eeg_project
export LAPLACE_BASE_SUBDIR=laplacian_basis_brms_p20_gaussian_hardlabel_constant_long/k20
export LAPLACE_REPAIR_SUBDIR=laplacian_basis_brms_p20_gaussian_hardlabel_constant_repair_11_19/k20
export LAPLACE_MERGED_SUBDIR=laplacian_basis_brms_p20_gaussian_hardlabel_constant_merged_repair/k20
export LAPLACE_REPAIR_TASKS=11,19
export LAPLACE_ITER=3000
export LAPLACE_WARMUP=1000
export LAPLACE_BASE_ADAPT_TAG=0p999
export LAPLACE_REPAIR_ADAPT_TAG=0p9999
export LAPLACE_MERGED_ADAPT_TAG=0p9999

cd "$PROJECT_ROOT"
Rscript code/16_merge_laplacian_repair_summary.R
