#!/bin/bash
#SBATCH --job-name=bgam_mc_p5_gauss_hl_sem
#SBATCH --partition=ultralong
#SBATCH --array=1-1
#SBATCH --time=10-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH -o logs/bgam_mc_p5_gauss_hl_sem_%A_%a.log
#SBATCH -e logs/bgam_mc_p5_gauss_hl_sem_%A_%a.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=fjhmrqx@126.com

mkdir -p logs
module load R/4.4.2-gcc141-base

export PROJECT_ROOT=/work/smjifan1/bgam_eeg_project
export BRMS_MEEG=$PROJECT_ROOT/external_brms_meeg/code

export MC_ALPHA=0.05
export MC_DATA_SUBDIR=montecarlo_p5_gaussian_hardlabel
export MC_COMPARE_SUBDIR=montecarlo_p5_gaussian_hardlabel_sem_test
export MC_SUMMARY_SUBDIR=montecarlo_p5_gaussian_hardlabel_sem_test

cd "$PROJECT_ROOT"
Rscript code/09_run_montecarlo_sem.R "$SLURM_ARRAY_TASK_ID"
