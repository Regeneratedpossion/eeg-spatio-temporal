# Configuration for the Bayesian GAM fitting pipeline.

suppressPackageStartupMessages(library(brms))


# 1. Experiment and model settings --------------------------------

ds_cr_seed_config <- function(project_root, cond = Sys.getenv("COND")) {

    stopifnot(
        "COND must be one of 1, 2, 3, 4" =
            cond %in% c("1", "2", "3", "4")
    )

    project_root <- normalizePath(project_root)

    experiment_name <- sprintf(
        "single_node_cond%s_ds_cr_seeds001_100", cond
    )
    results_dir <- file.path(
        project_root, "results", experiment_name
    )

    list(
        project_root = project_root,
        cond = cond,
        experiment_name = experiment_name,
        experiment_dir = results_dir,
        results_dir = results_dir,
        seeds = 1:100,
        stage1_seeds = c(
            2L, 11L, 22L, 23L, 27L,
            45L, 50L, 75L, 96L, 97L
        ),
        scenario = paste0("cond_", cond),
        radius_rule = "half-max",
        n_ppt = 20L,
        expected_channels = 61L,
        expected_time_points = 51L,
        expected_active_electrodes = 9L,
        predictor_scaling = "none",
        condition_contrasts = c(-0.5, 0.5),
        chains = 4L,
        cores = 4L,
        threads_per_chain = 1L,
        iter = 3000L,
        warmup = 1000L,
        adapt_delta = 0.80,
        max_treedepth = 10L,
        metric = "diag_e",
        spatial_basis = "ds",
        temporal_basis = "cr",
        k = c(16L, 12L),
        prior_mode = "brms_default"
    )
}


# 2. Seed validation ---------------------------------------------

validate_ds_cr_seed <- function(seed, config) {

    value <- suppressWarnings(as.numeric(seed))

    if (
        length(value) != 1L ||
        !is.finite(value) ||
        value != floor(value) ||
        !(value %in% config$seeds)
    ) {
        stop("seed must be one integer in 1:100", call. = FALSE)
    }

    as.integer(value)
}


# 3. Input data paths --------------------------------------------

ds_cr_seed_data_paths <- function(config, seed) {

    seed <- validate_ds_cr_seed(seed, config)

    data_subdir <- sprintf(
        "single_node_cond%s_seeds001_100", config$cond
    )
    data_dir <- file.path(
        config$project_root, "data", data_subdir
    )

    list(
        data_subdir = data_subdir,
        sim_file = file.path(
            data_dir, sprintf("sim_%03d.rds", seed)
        ),
        truth_file = file.path(
            data_dir, sprintf("gt_%03d.rds", seed)
        )
    )
}


# 4. Model output paths ------------------------------------------

ds_cr_seed_fit_paths <- function(config, seed) {

    seed <- validate_ds_cr_seed(seed, config)

    seed_dir <- file.path(
        config$results_dir, sprintf("seed_%03d", seed)
    )

    list(
        seed_dir = seed_dir,
        fit_file = file.path(seed_dir, "fit.rds"),
        config_file = file.path(seed_dir, "config.rds"),
        runtime_file = file.path(seed_dir, "runtime.csv"),
        nuts_file = file.path(seed_dir, "nuts_summary.csv"),
        convergence_file = file.path(seed_dir, "convergence.csv"),
        sampler_dir = file.path(seed_dir, "sampler"),
        cache_dir = file.path(seed_dir, "cmdstan_cache"),
        complete_file = file.path(seed_dir, "complete.flag")
    )
}


# 5. Bayesian GAM formula ----------------------------------------

ds_cr_seed_formula <- function() {

    brms::bf(
        eeg_mean | se(eeg_sem, sigma = TRUE) ~
            condition +
            (1 + condition | participant) +
            t2(
                xproj, yproj, time,
                by = condition,
                bs = c("ds", "cr"),
                k = c(16, 12),
                d = c(2, 1),
                m = list(c(1, 0.5), c(2, 2)),
                full = FALSE
            )
    )
}


# 6. Input checks ------------------------------------------------

# Checks the input and output files, then returns the model data and
# the ground truth of one seed.
check_fit_inputs <- function(config, data_paths, fit_paths, seed) {

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

    sim_data <- readRDS(data_paths$sim_file)
    truth <- readRDS(data_paths$truth_file)

    model_data <- aggregate_trials(sim_data)
    model_data$participant <- factor(model_data$participant)

    contrasts(model_data$condition) <- matrix(
        config$condition_contrasts,
        ncol = 1L,
        dimnames = list(
            c("cond1", "cond2"),
            "condition_sum"
        )
    )

    stopifnot(
        identical(truth$scenario, config$scenario),
        identical(truth$radius_rule, config$radius_rule),
        length(truth$active_electrodes) ==
            config$expected_active_electrodes,
        nlevels(model_data$participant) == config$n_ppt,
        dplyr::n_distinct(model_data$channel) ==
            config$expected_channels,
        dplyr::n_distinct(model_data$time) ==
            config$expected_time_points,
        all(is.finite(as.matrix(model_data[
            c("eeg_mean", "eeg_sem", "xproj", "yproj", "time")
        ]))),
        all(model_data$eeg_sem > 0),
        identical(
            unname(contrasts(model_data$condition)[, 1L]),
            config$condition_contrasts
        )
    )

    list(model_data = model_data, truth = truth)
}


# 7. CmdStan setup -----------------------------------------------

# If CMDSTAN_PATH is empty, use the CmdStan that cmdstanr finds.
setup_cmdstan <- function(config) {

    shared_cache_dir <- file.path(
        config$project_root, "results", "shared_model_cache"
    )

    dir.create(
        shared_cache_dir,
        recursive = TRUE,
        showWarnings = FALSE
    )

    cmdstan_path <- Sys.getenv("CMDSTAN_PATH", unset = "")

    if (nzchar(cmdstan_path)) {
        cmdstanr::set_cmdstan_path(cmdstan_path)
    }

    Sys.setenv(CMDSTAN_CACHE_DIR = shared_cache_dir)

    options(
        cmdstanr_write_stan_file_dir = shared_cache_dir,
        mc.cores = config$cores
    )

    invisible(shared_cache_dir)
}


# 8. Fit diagnostics ---------------------------------------------

# Writes runtime, NUTS and convergence summaries of one fit.
write_fit_diagnostics <- function(fit, config, fit_paths, seed,
                                  started, finished) {

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

    draw_summary <- posterior::summarise_draws(
        posterior::as_draws_df(fit)
    )

    write.csv(data.frame(
        seed = seed,
        max_rhat = max(draw_summary$rhat, na.rm = TRUE),
        min_ess_bulk = min(draw_summary$ess_bulk, na.rm = TRUE),
        min_ess_tail = min(draw_summary$ess_tail, na.rm = TRUE)
    ), fit_paths$convergence_file, row.names = FALSE)

    invisible(NULL)
}
