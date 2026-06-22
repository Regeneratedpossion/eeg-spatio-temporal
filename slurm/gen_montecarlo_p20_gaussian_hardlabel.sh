#!/bin/bash
#SBATCH --job-name=gen_p20_gauss_hl
#SBATCH --partition=short
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH -o logs/gen_p20_gauss_hl_%j.log
#SBATCH -e logs/gen_p20_gauss_hl_%j.log

mkdir -p logs
module load R/4.4.2-gcc141-base

export PROJECT_ROOT=/work/smjifan1/bgam_eeg_project
export MC_N_PPT=20
export MC_SIGNAL_MODE=continuous_label
export MC_DATA_SUBDIR=montecarlo_p20_gaussian_hardlabel

cd "$PROJECT_ROOT"
Rscript code/00_generate_montecarlo.R 1
