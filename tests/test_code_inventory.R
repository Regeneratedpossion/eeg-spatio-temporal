project_root <- Sys.getenv("PROJECT_ROOT", unset = getwd())

expected_r <- c(
    sprintf("code/%02d_%s.R", 0:22, c(
        "check_environment",
        "simulate_eeg_data",
        "generate_montecarlo",
        "fit_bgam",
        "compare_methods",
        "aggregate_results",
        "validate_paired_radius_pilot",
        "fit_paired_radius_bgam",
        "diagnose_paired_radius_bgam",
        "aggregate_paired_radius_bgam",
        "compute_delta_sesoi",
        "plot_paired_radius_truth",
        "plot_2sd_bgam_diagnostics",
        "plot_simulated_eeg",
        "plot_bgam_fit",
        "plot_posterior_probability",
        "plot_topomaps",
        "prior_predictive_check",
        "basis_screening",
        "basis_screening_full",
        "laplacian_sensitivity",
        "posterior_odds_sensitivity",
        "significance_alpha_sensitivity"
    )),
    "code/archive/09_run_montecarlo_sem.R",
    "code/archive/11_run_bgam_diagnostics.R",
    "code/archive/12_validate_heterogeneity_bgam.R",
    "code/archive/compute_trial_effective_sample_size.R",
    "code/archive/migrate_ground_truth.R",
    "code/ground_truth.R",
    "code/laplacian/config.R",
    "code/laplacian/diagnose_continuous.R",
    "code/laplacian/diagnose_hard.R",
    "code/laplacian/fit_basis_scores.R",
    "code/laplacian/reconstruct_effect.R",
    "code/paper_config.R",
    "code/utils.R"
)

actual_r <- list.files(
    file.path(project_root, "code"),
    pattern = "[.]R$",
    recursive = TRUE,
    full.names = FALSE
)
actual_r <- sort(file.path("code", actual_r))
stopifnot(identical(actual_r, sort(expected_r)))

numbered_r <- list.files(
    file.path(project_root, "code"),
    pattern = "^[0-9]{2}_.+[.]R$",
    full.names = FALSE
)
observed_numbers <- as.integer(sub("^([0-9]{2})_.*$", "\\1", numbered_r))
stopifnot(identical(sort(observed_numbers), 0:22))

removed_dirs <- file.path(
    project_root,
    "code",
    c("pipeline", "pilot", "diagnostics", "reporting",
      "sensitivity", "truth", "config")
)
stopifnot(!any(dir.exists(removed_dirs)))

expected_slurm <- sort(c(
    "slurm/formal/aggregate_paired_radius.sh",
    "slurm/formal/fit_paired_radius_dortmund.sh",
    "slurm/formal/fit_paired_radius_hpc2n.sh",
    "slurm/formal/generate_paired_radius.sh",
    "slurm/pilot/run_paired_radius_pilot.sh"
))
actual_slurm <- list.files(
    file.path(project_root, "slurm"),
    pattern = "[.]sh$",
    recursive = TRUE,
    full.names = FALSE
)
actual_slurm <- sort(file.path("slurm", actual_slurm))
stopifnot(identical(actual_slurm, expected_slurm))

cat("code inventory tests passed\n")
