#!/bin/bash
#SBATCH --job-name=bgam_basis_full_ab
#SBATCH --partition=ultralong
#SBATCH --array=1-2
#SBATCH --time=10-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH -o logs/bgam_basis_full_ab_%A_%a.log
#SBATCH -e logs/bgam_basis_full_ab_%A_%a.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=fjhmrqx@126.com

mkdir -p logs
module load R/4.4.2-gcc141-base

export PROJECT_ROOT=/work/smjifan1/bgam_eeg_project

cd "$PROJECT_ROOT"
Rscript code/11_basis_screening.R
