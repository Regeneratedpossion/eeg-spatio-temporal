###########################################################
# 00_generate_montecarlo.R
#
# Generate Monte Carlo datasets for later BGAM fitting and method comparison.
# Signal mode is controlled by MC_SIGNAL_MODE:
#   continuous_label = continuous Gaussian signal; threshold labels active electrodes
#   hard_signal      = thresholded signal; outside-threshold weights are 0
#
# Run:
#   Rscript code/00_generate_montecarlo.R
#
# Optional settings:
#   MC_N_PPT=20          participant count
#   MC_SCENARIO=legacy   legacy, low, medium, or high
#   MC_DATA_SUBDIR=...   output folder under data/
#   MC_SIGNAL_MODE=...   continuous_label or hard_signal
###########################################################

source(file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()), "code/utils.R"))

# Load make_dataset().
source(file.path(project_root(), "code/00_generate_data.R"))

# Generate and save one simulation and its ground truth.
generate_dataset_for_seed <- function(seed, out_dir, n_ppt, scenario,
                                      signal_mode) {
    sim_path <- file.path(out_dir, sprintf("sim_%03d.rds", seed))
    gt_path <- file.path(out_dir, sprintf("gt_%03d.rds", seed))

    dataset <- make_dataset(
        seed = seed,
        n_ppt = n_ppt,
        scenario = scenario,
        signal_mode = signal_mode
    )
    saveRDS(dataset$sim_data, sim_path)
    saveRDS(dataset$ground_truth, gt_path)
}

scenario_manifest <- function(n_ppt, scenario, signal_mode) {
    params <- simulation_parameters(
        n_ppt = n_ppt,
        scenario = scenario,
        signal_mode = signal_mode
    )
    list(
        scenario = scenario,
        n_ppt = n_ppt,
        n_trials = params$n_trials,
        signal_mode = signal_mode,
        onset_support_ms = params$ronset,
        bump_amplitude = params$bump_amplitude,
        omega_A = params$omega_A,
        omega_sigma = params$omega_sigma
    )
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

main <- function() {
    n_sim <- as.integer(c(commandArgs(trailingOnly = TRUE), 30L)[1])
    n_ppt <- as.integer(Sys.getenv("MC_N_PPT", unset = "20"))
    scenario <- Sys.getenv("MC_SCENARIO", unset = "legacy")
    signal_mode <- Sys.getenv("MC_SIGNAL_MODE", unset = "continuous_label")
    out_dir <- mc_paths()$data_dir
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    write_or_validate_manifest(
        out_dir,
        scenario_manifest(n_ppt, scenario, signal_mode)
    )

    invisible(lapply(seq_len(n_sim), generate_dataset_for_seed,
                     out_dir = out_dir, n_ppt = n_ppt,
                     scenario = scenario, signal_mode = signal_mode))
}

main()
