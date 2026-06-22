#!/bin/bash
#SBATCH --job-name=bgam_signal_p5
#SBATCH --partition=ultralong
#SBATCH --array=1-2
#SBATCH --time=10-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH -o logs/bgam_signal_p5_%A_%a.log
#SBATCH -e logs/bgam_signal_p5_%A_%a.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=fjhmrqx@126.com

mkdir -p logs
module load R/4.4.2-gcc141-base

export PROJECT_ROOT=/work/smjifan1/bgam_eeg_project
export BRMS_MEEG=$PROJECT_ROOT/external_brms_meeg/code
export MC_N_PPT=5
export K_XY=20
export K_T=20

if [ "$SLURM_ARRAY_TASK_ID" = "1" ]; then
    export MC_SIGNAL_MODE=continuous_label
    export MC_DATA_SUBDIR=signal_sensitivity_p5_continuous_label
    export MC_COMPARE_SUBDIR=signal_sensitivity_p5_continuous_label
    export MC_SUMMARY_SUBDIR=signal_sensitivity_p5_continuous_label
else
    export MC_SIGNAL_MODE=hard_signal
    export MC_DATA_SUBDIR=signal_sensitivity_p5_hard_signal
    export MC_COMPARE_SUBDIR=signal_sensitivity_p5_hard_signal
    export MC_SUMMARY_SUBDIR=signal_sensitivity_p5_hard_signal
fi

cd "$PROJECT_ROOT"

# Generate one seed under the selected signal definition, then fit and compare.
Rscript code/00_generate_montecarlo.R 1
Rscript code/09_run_montecarlo_sem.R 1
