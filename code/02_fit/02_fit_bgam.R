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


# 2. Input and output paths ----------------------------------

required_inputs <- c(
    data_paths$sim_file,
    data_paths$truth_file
)

missing_inputs <- required_inputs[!file.exists(required_inputs)]

if (length(missing_inputs)) {
    stop(
        "Missing seed inputs: ",
        paste(missing_inputs, collapse = ", "),
        call. = FALSE
    )
}

if (
    dir.exists(fit_paths$seed_dir) &&
    length(list.files(
        fit_paths$seed_dir, all.files = TRUE, no.. = TRUE
    ))
) {
    stop(
        "Fit output is non-empty for seed ", seed,
        "; refusing to overwrite",
        call. = FALSE
    )
}

dir.create(
    fit_paths$sampler_dir,
    recursive = TRUE,
    showWarnings = FALSE
)


# 3. CmdStan setup -------------------------------------------

shared_cache_dir <- file.path(
    root, "results", "shared_model_cache"
)

dir.create(
    shared_cache_dir,
    recursive = TRUE,
    showWarnings = FALSE
)

# If CMDSTAN_PATH is empty, use the CmdStan that cmdstanr finds.
cmdstan_path <- Sys.getenv("CMDSTAN_PATH", unset = "")

if (nzchar(cmdstan_path)) {
    cmdstanr::set_cmdstan_path(cmdstan_path)
}

Sys.setenv(CMDSTAN_CACHE_DIR = shared_cache_dir)

options(
    cmdstanr_write_stan_file_dir = shared_cache_dir,
    mc.cores = config$cores
)


# 4. Prepare model data --------------------------------------

sim_data <- readRDS(data_paths$sim_file)
truth <- readRDS(data_paths$truth_file)

model_data <- aggregate_trials(sim_data) %>%
    mutate(participant = factor(participant))

contrasts(model_data$condition) <- matrix(
    config$condition_contrasts,
    ncol = 1L,
    dimnames = list(
        c("cond1", "cond2"),
        "condition_sum"
    )
)


# 5. Validate model data -------------------------------------

stopifnot(
    identical(truth$scenario, config$scenario),
    identical(truth$radius_rule, config$radius_rule),
    length(truth$active_electrodes) ==
        config$expected_active_electrodes,
    nlevels(model_data$participant) == config$n_ppt,
    n_distinct(model_data$channel) == config$expected_channels,
    n_distinct(model_data$time) == config$expected_time_points,
    all(is.finite(as.matrix(model_data[
        c("eeg_mean", "eeg_sem", "xproj", "yproj", "time")
    ]))),
    all(model_data$eeg_sem > 0),
    identical(
        unname(contrasts(model_data$condition)[, 1L]),
        config$condition_contrasts
    )
)


# 6. Save fitting configuration ------------------------------

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
    input_md5 = tools::md5sum(required_inputs)
), fit_paths$config_file)


# 7. Fit the Bayesian GAM ------------------------------------

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


# 8. Runtime diagnostics -------------------------------------

wall_seconds <- as.numeric(
    difftime(finished, started, units = "secs")
)

write.csv(data.frame(
    seed = seed,
    hostname = Sys.info()[["nodename"]],
    started = as.character(started),
    finished = as.character(finished),
    wall_seconds = wall_seconds,
    wall_hours = wall_seconds / 3600,
    allocated_cpus = config$cores,
    cpu_hours = wall_seconds / 3600 * config$cores,
    threads_per_chain = config$threads_per_chain
), fit_paths$runtime_file, row.names = FALSE)


# 9. NUTS diagnostics ----------------------------------------

nuts <- brms::nuts_params(fit)

write.csv(data.frame(
    seed = seed,
    divergences = sum(
        nuts$Value[nuts$Parameter == "divergent__"]
    ),
    treedepth_hits = sum(
        nuts$Value[nuts$Parameter == "treedepth__"] >=
            config$max_treedepth
    ),
    mean_leapfrog = mean(
        nuts$Value[nuts$Parameter == "n_leapfrog__"]
    ),
    mean_step_size = mean(
        nuts$Value[nuts$Parameter == "stepsize__"]
    )
), fit_paths$nuts_file, row.names = FALSE)


# 10. Convergence diagnostics --------------------------------

draw_summary <- posterior::summarise_draws(
    posterior::as_draws_df(fit)
)

write.csv(data.frame(
    seed = seed,
    max_rhat = max(draw_summary$rhat, na.rm = TRUE),
    min_ess_bulk = min(draw_summary$ess_bulk, na.rm = TRUE),
    min_ess_tail = min(draw_summary$ess_tail, na.rm = TRUE)
), fit_paths$convergence_file, row.names = FALSE)

writeLines(
    as.character(finished),
    fit_paths$complete_file
)

cat("completed ds x cr fit for seed ", seed, "\n", sep = "")
