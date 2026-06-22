#!/bin/bash
#SBATCH -A hpc2n2025-235
#SBATCH --job-name=lap_rec_hpc2n
#SBATCH --array=1-4
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH -o logs/lap_rec_hpc2n_%A_%a.log
#SBATCH -e logs/lap_rec_hpc2n_%A_%a.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=fjhmrqx@126.com

mkdir -p logs
module load GCC/13.2.0
module load R/4.4.1

export PROJECT_ROOT=/proj/nobackup/hpc2n2025-235/jiahui/bgam_eeg_project
export R_LIBS_USER=/proj/nobackup/hpc2n2025-235/jiahui/Rlibs/4.4
export LAPLACE_DATA_SUBDIR=montecarlo_p5_gaussian_hardlabel
export LAPLACE_SEED=1
export LAPLACE_TIME_K=10
export LAPLACE_ITER=1000
export LAPLACE_WARMUP=500
export LAPLACE_ADAPT_TAG=0p999

K_VALUES=(10 15 20 30)
K="${K_VALUES[$((SLURM_ARRAY_TASK_ID - 1))]}"

export LAPLACE_K="$K"
export LAPLACE_FIT_SUBDIR="laplacian_basis_brms_ladder_hpc2n/k${K}"

cd "$PROJECT_ROOT"
Rscript code/13_laplacian_recovery_diagnostics.R
