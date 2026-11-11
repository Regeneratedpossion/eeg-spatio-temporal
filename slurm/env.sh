# Shared R environment for all Slurm scripts on LIDO.
# Sourced by each script: source "$PROJECT_ROOT/slurm/env.sh"
# Rscript --vanilla ignores .Renviron, so the package library must be set in the shell.
# Export R_LIBS_USER before you run sbatch, for example in your ~/.bashrc.

module purge
module load R/4.4.2-gcc141-base

: "${R_LIBS_USER:?Set R_LIBS_USER before sbatch, e.g. export R_LIBS_USER=/path/to/your/R/library}"
export R_LIBS_USER
