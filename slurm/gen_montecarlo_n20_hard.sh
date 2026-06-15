#!/bin/bash
#SBATCH --job-name=gen_mc_n20_hard
#SBATCH --partition=ultralong
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH -o logs/gen_mc_n20_hard_%j.log
#SBATCH -e logs/gen_mc_n20_hard_%j.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=fjhmrqx@126.com

mkdir -p logs
module load R/4.4.2-gcc141-base

export PROJECT_ROOT=/work/smjifan1/bgam_eeg_project
export MC_DATA_SUBDIR=montecarlo_n20_hard

cd "$PROJECT_ROOT"
Rscript code/00_generate_montecarlo.R 30
