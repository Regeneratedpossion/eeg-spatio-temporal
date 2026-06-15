#!/bin/bash
#SBATCH --job-name=bgam_diag
#SBATCH --partition=med
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH -o logs/bgam_diag_%j.log
#SBATCH -e logs/bgam_diag_%j.log

# Run from project root: sbatch slurm/run_diagnostics.sh
mkdir -p logs

module load R/4.4.2-gcc141-base

export PROJECT_ROOT=/work/smjifan1/bgam_eeg_project
export BRMS_MEEG=$PROJECT_ROOT/external_brms_meeg/code
cd "$PROJECT_ROOT"

Rscript code/02_diagnostics.R
