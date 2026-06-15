#!/bin/bash
#SBATCH --job-name=bgam_fit
#SBATCH --partition=ultralong
#SBATCH --time=4-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH -o logs/bgam_fit_%j.log
#SBATCH -e logs/bgam_fit_%j.log

# Run from project root: sbatch slurm/run_fit.sh
mkdir -p logs

module load R/4.4.2-gcc141-base

export PROJECT_ROOT=/work/smjifan1/bgam_eeg_project
export BRMS_MEEG=$PROJECT_ROOT/external_brms_meeg/code
cd "$PROJECT_ROOT"

Rscript code/01_fit_bgam.R
