source(file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()), "code/utils.R"))
source(file.path(
    Sys.getenv("PROJECT_ROOT", unset = getwd()),
    "code/config/paper_config.R"
))
fit_dependency <- file.path(
    project_root(), "code", "paired_radius_fit_bgam_dependency.R"
)
if (!file.exists(fit_dependency)) {
    fit_dependency <- file.path(project_root(), "code/01_fit_bgam.R")
}
source(fit_dependency)

paired_radius_bgam_config <- function() {
    paper_analysis_config()$model
}

paired_radius_bgam_specs <- function() {
    data.frame(
        radius_rule = c("1sd", "half-max", "2sd"),
        data_subdir = c(
            "montecarlo_r1sd",
            "montecarlo_rhalfmax",
            "montecarlo_r2sd"
        ),
        jump_weight = c(0.7038613, 0.5094451, 0.1932698),
        expected_active = c(5L, 9L, 21L),
        stringsAsFactors = FALSE
    )
}

paired_radius_output_dir <- function(radius_rule) {
    file.path(project_root(), "results", "paired_radius_bgam_pilot",
              radius_rule)
}

run_paired_radius_bgam_fit <- function(radius_rule) {
    config <- paired_radius_bgam_config()
    specs <- paired_radius_bgam_specs()
    spec <- specs[specs$radius_rule == radius_rule, , drop = FALSE]
    if (nrow(spec) != 1L) {
        stop("Unknown radius_rule: ", radius_rule, call. = FALSE)
    }

    sim_data <- readRDS(data_path(spec$data_subdir, "sim_001.rds"))
    ground_truth <- readRDS(data_path(spec$data_subdir, "gt_001.rds"))
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
        seed = config$seed,
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
    saveRDS(fit, file.path(out_dir, "fit_seed001.rds"))
    write.csv(
        data.frame(
            radius_rule = radius_rule,
            started = as.character(started),
            finished = as.character(Sys.time()),
            elapsed_hours = as.numeric(Sys.time() - started, units = "hours")
        ),
        file.path(out_dir, "fit_runtime.csv"),
        row.names = FALSE
    )
    invisible(fit)
}

if (!interactive() && sys.nframe() == 0) {
    args <- commandArgs(trailingOnly = TRUE)
    radius_rule <- if (length(args)) args[[1]] else "1sd"
    run_paired_radius_bgam_fit(radius_rule)
}
