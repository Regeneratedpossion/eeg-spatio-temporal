###########################################################
# 20_laplacian_sensitivity.R
#
# Entrypoint for the fixed K=20 + constant Laplacian workflow.
#
# Actions:
#   LAPLACE_ACTION=fit                 fit basis-score temporal GAMMs
#   LAPLACE_ACTION=diagnose_hard       threshold-based recovery summary
#   LAPLACE_ACTION=diagnose_continuous continuous Gaussian-field recovery
###########################################################

source(file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()), "code/utils.R"))
source(file.path(project_root(), "code/laplacian/config.R"))
source(file.path(project_root(), "code/laplacian/fit_basis_scores.R"))
source(file.path(project_root(), "code/laplacian/reconstruct_effect.R"))
source(file.path(project_root(), "code/laplacian/diagnose_hard.R"))
source(file.path(project_root(), "code/laplacian/diagnose_continuous.R"))

setup_cmdstan_cache()

laplace_action <- Sys.getenv("LAPLACE_ACTION", "fit")

switch(
    laplace_action,
    fit = run_fit(),
    diagnose_hard = run_hard_diagnostics(),
    diagnose_continuous = run_continuous_diagnostics(),
    stop("LAPLACE_ACTION must be one of: fit, diagnose_hard, diagnose_continuous.")
)
