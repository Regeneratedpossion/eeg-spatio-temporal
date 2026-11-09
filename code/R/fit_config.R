suppressPackageStartupMessages(library(brms))

ds_cr_seed_config <- function(project_root, cond = Sys.getenv("COND")) {
    stopifnot("COND must be one of 1, 2, 3, 4" = cond %in% c("4", "3", "2", "1"))
    project_root <- normalizePath(project_root)
    experiment_name <- sprintf("single_node_cond%s_ds_cr_seeds001_100", cond)
    list(
        project_root = project_root,
        cond = cond,
        experiment_name = experiment_name,
        experiment_dir = file.path(project_root, "results", experiment_name),
        results_dir = file.path(project_root, "results", experiment_name),
        seeds = 1:100,
        stage1_seeds = c(2L, 11L, 22L, 23L, 27L,
                         45L, 50L, 75L, 96L, 97L),
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

validate_ds_cr_seed <- function(seed, config) {
    value <- suppressWarnings(as.numeric(seed))
    if (length(value) != 1L || !is.finite(value) ||
        value != floor(value) || !(value %in% config$seeds)) {
        stop("seed must be one integer in 1:100", call. = FALSE)
    }
    as.integer(value)
}

ds_cr_seed_data_paths <- function(config, seed) {
    seed <- validate_ds_cr_seed(seed, config)
    data_subdir <- sprintf("single_node_cond%s_seeds001_100", config$cond)
    data_dir <- file.path(config$project_root, "data", data_subdir)
    list(
        data_subdir = data_subdir,
        sim_file = file.path(data_dir, sprintf("sim_%03d.rds", seed)),
        truth_file = file.path(data_dir, sprintf("gt_%03d.rds", seed))
    )
}

ds_cr_seed_fit_paths <- function(config, seed) {
    seed <- validate_ds_cr_seed(seed, config)
    seed_dir <- file.path(config$results_dir, sprintf("seed_%03d", seed))
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

ds_cr_seed_formula <- function() {
    brms::bf(
        eeg_mean | se(eeg_sem, sigma = TRUE) ~
            condition + (1 + condition | participant) +
            t2(
                xproj, yproj, time, by = condition,
                bs = c("ds", "cr"), k = c(16, 12), d = c(2, 1),
                m = list(c(1, 0.5), c(2, 2)), full = FALSE
            )
    )
}
