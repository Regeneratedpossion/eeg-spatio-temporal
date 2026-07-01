#!/bin/bash
#SBATCH --job-name=bgam_ab_defctl_rescue_A
#SBATCH --partition=ultralong
#SBATCH --array=1,3
#SBATCH --time=10-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH -o logs/bgam_ab_defctl_rescue_A_%A_%a.log
#SBATCH -e logs/bgam_ab_defctl_rescue_A_%A_%a.log
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=fjhmrqx@126.com

mkdir -p logs /work/smjifan1/tmp
module load R/4.4.2-gcc141-base

export PROJECT_ROOT=/work/smjifan1/bgam_eeg_project
export TMPDIR=/work/smjifan1/tmp
export TMP=/work/smjifan1/tmp
export TEMP=/work/smjifan1/tmp
export SCREEN_SEED=1
export SCREEN_ITER=3000
export SCREEN_WARMUP=1000
export SCREEN_ADAPT_DELTA=0.80
export SCREEN_MAX_TREEDEPTH=10
export SCREEN_OUT_SUBDIR=basis_screening_full_ab_default_control_rescue_A

case "$SLURM_ARRAY_TASK_ID" in
  1)
    export SCREEN_DATA_SUBDIR=montecarlo_p5_gaussian_hardlabel
    export INNER_TASK_ID=1
    ;;
  3)
    export SCREEN_DATA_SUBDIR=montecarlo_p20_gaussian_hardlabel
    export INNER_TASK_ID=1
    ;;
  *)
    echo "Unsupported rescue array task: $SLURM_ARRAY_TASK_ID" >&2
    exit 2
    ;;
esac

cd "$PROJECT_ROOT"
echo "rescue data=$SCREEN_DATA_SUBDIR inner_task=$INNER_TASK_ID adapt=$SCREEN_ADAPT_DELTA td=$SCREEN_MAX_TREEDEPTH TMPDIR=$TMPDIR"
SLURM_ARRAY_TASK_ID="$INNER_TASK_ID" Rscript code/18_basis_screening_full_ab_control.R
