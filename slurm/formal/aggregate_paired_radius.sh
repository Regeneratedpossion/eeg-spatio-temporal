#!/bin/bash
#SBATCH --job-name=paired_radius_aggregate
#SBATCH --partition=med
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -o logs/paired_radius_aggregate_%j.log
#SBATCH -e logs/paired_radius_aggregate_%j.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=fjhmrqx@126.com

set -euo pipefail

module load R/4.4.2-gcc141-base
export PROJECT_ROOT=/work/smjifan1/bgam_eeg_project
export R_LIBS_USER=/work/smjifan1/R/4.4.2:/home/smjifan1/R/x86_64-pc-linux-gnu-library/4.4
cd "$PROJECT_ROOT"
Rscript code/10_aggregate_paired_radius_bgam.R
