#!/bin/bash
#SBATCH --job-name=lap_recovery_diag
#SBATCH --partition=ultralong
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH -o logs/lap_recovery_diag_%j.log
#SBATCH -e logs/lap_recovery_diag_%j.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=fjhmrqx@126.com

mkdir -p logs
module load R/4.4.2-gcc141-base

export PROJECT_ROOT=/work/smjifan1/bgam_eeg_project
export LAPLACE_DATA_SUBDIR=montecarlo_p5_gaussian_hardlabel
export LAPLACE_SEED=1
export LAPLACE_K=10
export LAPLACE_ADAPT_TAG=0p999

cd "$PROJECT_ROOT"
Rscript code/13_laplacian_recovery_diagnostics.R
