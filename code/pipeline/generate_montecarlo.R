###########################################################
# 00_generate_montecarlo.R
#
# Generate paired hard-radius Monte Carlo datasets.
#
# Run:
#   Rscript code/pipeline/generate_montecarlo.R
#
# Optional settings:
#   MC_N_PPT=20          participant count
#   MC_SCENARIO=high     legacy, low, medium, or high
#   MC_SIGNAL_MODE=hard_signal
###########################################################

source(file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()), "code/utils.R"))

# Load the paired data generator.
source(file.path(project_root(), "code/pipeline/simulation.R"))

paired_scenario_manifest <- function(params, spec, seeds) {
    list(
        scenario = params$scenario,
        n_ppt = params$n_ppt,
        n_trials = params$n_trials,
        signal_mode = "hard_signal",
        active_radius = spec$active_radius,
        radius_rule = spec$radius_rule,
        replicate_seeds = as.integer(seeds),
        paired_design_id = paste0(
            params$scenario, "-p", params$n_ppt,
            "-r", params$n_trials,
            "-seeds-", min(seeds), "-", max(seeds)
        ),
        onset_support_ms = params$ronset,
        bump_amplitude = params$bump_amplitude,
        omega_A = params$omega_A,
        omega_sigma = params$omega_sigma
    )
}

save_paired_replicate <- function(seed, datasets, output_dirs) {
    stopifnot(setequal(names(datasets), names(output_dirs)))
    for (radius_rule in names(output_dirs)) {
        dataset <- datasets[[radius_rule]]
        out_dir <- output_dirs[[radius_rule]]
        saveRDS(
            dataset$sim_data,
            file.path(out_dir, sprintf("sim_%03d.rds", seed))
        )
        saveRDS(
            dataset$ground_truth,
            file.path(out_dir, sprintf("gt_%03d.rds", seed))
        )
    }
    invisible(NULL)
}

write_or_validate_manifest <- function(out_dir, manifest) {
    manifest_path <- file.path(out_dir, "scenario_config.rds")
    if (file.exists(manifest_path)) {
        existing <- readRDS(manifest_path)
        if (!identical(existing, manifest)) {
            stop("Output directory contains a different scenario configuration: ",
                 out_dir)
        }
    } else {
        saveRDS(manifest, manifest_path)
    }
}

run_paired_radius_simulation <- function(seeds, n_ppt = 20L,
                                         scenario = "high") {
    params <- simulation_parameters(
        n_ppt = n_ppt,
        scenario = scenario,
        signal_mode = "hard_signal"
    )
    specs <- radius_specifications(params$sigma_space)
    output_dirs <- setNames(
        vapply(specs$data_subdir, data_path, character(1)),
        specs$radius_rule
    )
    for (i in seq_len(nrow(specs))) {
        out_dir <- output_dirs[[specs$radius_rule[i]]]
        dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
        write_or_validate_manifest(
            out_dir,
            paired_scenario_manifest(params, specs[i, ], seeds)
        )
    }

    for (seed in seeds) {
        message("Generating paired replicate ", seed, "...")
        datasets <- make_paired_datasets(
            seed = seed,
            n_ppt = n_ppt,
            scenario = scenario
        )
        save_paired_replicate(seed, datasets, output_dirs)
        rm(datasets)
        invisible(gc())
    }
    invisible(output_dirs)
}

main <- function() {
    n_sim <- as.integer(c(commandArgs(trailingOnly = TRUE), 30L)[1])
    n_ppt <- as.integer(Sys.getenv("MC_N_PPT", unset = "20"))
    scenario <- Sys.getenv("MC_SCENARIO", unset = "high")
    signal_mode <- Sys.getenv("MC_SIGNAL_MODE", unset = "hard_signal")
    if (!identical(signal_mode, "hard_signal")) {
        stop("MC_SIGNAL_MODE must be hard_signal", call. = FALSE)
    }
    run_paired_radius_simulation(seq_len(n_sim), n_ppt, scenario)
}

if (!interactive() && sys.nframe() == 0) {
    main()
}
