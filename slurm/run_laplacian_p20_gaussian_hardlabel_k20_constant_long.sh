#!/bin/bash
#SBATCH --job-name=lap_p20gh_k20_const
#SBATCH --partition=med
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH -o logs/lap_p20gh_k20_const_%j.log
#SBATCH -e logs/lap_p20gh_k20_const_%j.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=fjhmrqx@126.com

mkdir -p logs
module load R/4.4.2-gcc141-base

export PROJECT_ROOT=/work/smjifan1/bgam_eeg_project

cd "$PROJECT_ROOT"
Rscript code/12_laplacian_basis_brms_fixed.R
