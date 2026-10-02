source(file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()), "code/utils.R"))
source(file.path(
    Sys.getenv("PROJECT_ROOT", unset = getwd()),
    "code/paper_config.R"
))
fit_dependency <- file.path(
    project_root(), "code", "paired_radius_fit_bgam_dependency.R"
)
if (!file.exists(fit_dependency)) {
    fit_dependency <- file.path(project_root(), "code/03_fit_bgam.R")
}
source(fit_dependency)

paired_radius_bgam_config <- function() {
    paper_analysis_config()$model
}

paired_radius_bgam_specs <- function() {
    data.frame(
        radius_rule = "half-max",
        data_subdir = "montecarlo_rhalfmax",
        jump_weight = 0.5094451,
        expected_active = 9L,
        stringsAsFactors = FALSE
    )
}

paired_radius_output_dir <- function(radius_rule) {
    file.path(project_root(), "results", "paired_radius_bgam_pilot",
              radius_rule)
}

run_paired_radius_bgam_fit <- function(radius_rule, seed = 1L) {
    seed <- as.integer(seed)
    if (length(seed) != 1L || is.na(seed) || seed < 1L) {
        stop("seed must be one positive integer", call. = FALSE)
    }
    config <- paired_radius_bgam_config()
    specs <- paired_radius_bgam_specs()
    spec <- specs[specs$radius_rule == radius_rule, , drop = FALSE]
    if (nrow(spec) != 1L) {
        stop("Unknown radius_rule: ", radius_rule, call. = FALSE)
    }

    sim_data <- readRDS(data_path(
        spec$data_subdir, sprintf("sim_%03d.rds", seed)
    ))
    ground_truth <- readRDS(data_path(
        spec$data_subdir, sprintf("gt_%03d.rds", seed)
    ))
    stopifnot(
        identical(ground_truth$scenario, "high"),
        identical(ground_truth$radius_rule, radius_rule),
        length(unique(sim_data$participant)) == 20L,
        length(ground_truth$active_electrodes) == spec$expected_active
    )

    out_dir <- paired_radius_output_dir(radius_rule)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    config_path <- file.path(out_dir, "model_config.rds")
    if (file.exists(config_path)) {
        stopifnot(identical(readRDS(config_path), config))
    } else {
        saveRDS(config, config_path)
    }

    started <- Sys.time()
    fit <- fit_bgam(
        sim_data,
        seed = config$seed + seed - 1L,
        k_xy = config$k_xy,
        k_t = config$k_t,
        bs_space = config$bs_space,
        bs_time = config$bs_time,
        full = config$full,
        adapt_delta = config$adapt_delta,
        max_treedepth = config$max_treedepth,
        chains = config$chains,
        cores = config$chains,
        iter = config$iter,
        warmup = config$warmup,
        threads = 4L
    )
    saveRDS(fit, file.path(out_dir, sprintf("fit_seed%03d.rds", seed)))
    write.csv(
        data.frame(
            radius_rule = radius_rule,
            started = as.character(started),
            finished = as.character(Sys.time()),
            elapsed_hours = as.numeric(Sys.time() - started, units = "hours")
        ),
        file.path(out_dir, sprintf("fit_runtime_seed%03d.csv", seed)),
        row.names = FALSE
    )
    invisible(fit)
}

if (!interactive() && sys.nframe() == 0) {
    args <- commandArgs(trailingOnly = TRUE)
    radius_rule <- if (length(args)) args[[1]] else "half-max"
    seed <- if (length(args) >= 2L) as.integer(args[[2]]) else 1L
    run_paired_radius_bgam_fit(radius_rule, seed)
}
