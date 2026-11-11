# Fit one Bayesian GAM for one scenario and seed.

suppressPackageStartupMessages({
    library(brms)
    library(cmdstanr)
    library(dplyr)
    library(posterior)
})

# 1. Configuration -------------------------------------------

root <- normalizePath(
    Sys.getenv("PROJECT_ROOT", unset = here::here())
)

source(file.path(root, "code", "R", "utils.R"))
source(file.path(root, "code", "R", "fit_config.R"))

# The condition comes from the COND environment variable (1, 2, 3 or 4).
config <- ds_cr_seed_config(root)
seed <- validate_ds_cr_seed(
    Sys.getenv("SLURM_ARRAY_TASK_ID"), config
)

data_paths <- ds_cr_seed_data_paths(config, seed)
fit_paths <- ds_cr_seed_fit_paths(config, seed)


# 2. Check inputs and prepare model data --------------------

inputs <- check_fit_inputs(config, data_paths, fit_paths, seed)

model_data <- inputs$model_data
truth <- inputs$truth


# 3. CmdStan setup -------------------------------------------

setup_cmdstan(config)


# 4. Save fitting configuration ------------------------------

formula_ds_cr <- ds_cr_seed_formula()

saveRDS(list(
    seed = seed,
    data_subdir = data_paths$data_subdir,
    scenario = config$scenario,
    radius_rule = config$radius_rule,
    predictor_scaling = config$predictor_scaling,
    condition_contrasts = config$condition_contrasts,
    formula = formula_ds_cr,
    basis = c(config$spatial_basis, config$temporal_basis),
    k = config$k,
    chains = config$chains,
    cores = config$cores,
    threads_per_chain = config$threads_per_chain,
    iter = config$iter,
    warmup = config$warmup,
    adapt_delta = config$adapt_delta,
    max_treedepth = config$max_treedepth,
    metric = config$metric,
    prior_mode = config$prior_mode,
    hostname = Sys.info()[["nodename"]],
    input_md5 = tools::md5sum(c(
        data_paths$sim_file,
        data_paths$truth_file
    ))
), fit_paths$config_file)


# 5. Fit the Bayesian GAM ------------------------------------

started <- Sys.time()

fit <- brms::brm(
    formula = formula_ds_cr,
    data = model_data,
    family = gaussian(),
    chains = config$chains,
    cores = config$cores,
    iter = config$iter,
    warmup = config$warmup,
    seed = seed,
    backend = "cmdstanr",
    control = list(
        adapt_delta = config$adapt_delta,
        max_treedepth = config$max_treedepth
    ),
    metric = config$metric,
    save_warmup = TRUE,
    output_dir = fit_paths$sampler_dir,
    refresh = 100L
)

finished <- Sys.time()

saveRDS(fit, fit_paths$fit_file)


# 6. Diagnostics ---------------------------------------------

write_fit_diagnostics(
    fit, config, fit_paths, seed, started, finished
)

writeLines(
    as.character(finished),
    fit_paths$complete_file
)

cat("completed ds x cr fit for seed ", seed, "\n", sep = "")
