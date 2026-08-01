#!/bin/bash
#SBATCH --job-name=paired_radius_generate
#SBATCH --partition=ultralong
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH -o logs/paired_radius_generate_%j.log
#SBATCH -e logs/paired_radius_generate_%j.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=fjhmrqx@126.com

set -euo pipefail

mkdir -p logs
module load R/4.4.2-gcc141-base

export PROJECT_ROOT=/work/smjifan1/bgam_eeg_project
export BRMS_MEEG=$PROJECT_ROOT/external_brms_meeg/code
export MC_N_PPT=20
export MC_SCENARIO=high

N_SIM="${N_SIM:-20}"

cd "$PROJECT_ROOT"
Rscript code/02_generate_montecarlo.R "$N_SIM"
