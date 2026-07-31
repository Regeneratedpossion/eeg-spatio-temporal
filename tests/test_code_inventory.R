project_root <- Sys.getenv("PROJECT_ROOT", unset = getwd())

expected_r <- c(
    "code/archive/09_run_montecarlo_sem.R",
    "code/archive/compute_trial_effective_sample_size.R",
    "code/check_env.R",
    "code/config/paper_config.R",
    "code/diagnostics/02_diagnostics.R",
    "code/diagnostics/19_practical_null_fpr_rebound.R",
    "code/diagnostics/28_delta_sesoi_probe.R",
    "code/diagnostics/29_validate_heterogeneity_bgam.R",
    "code/diagnostics/35_diagnose_paired_radius_bgam_pilot.R",
    "code/pilot/00_generate_paired_radius_pilot.R",
    "code/pilot/00_validate_paired_radius_pilot.R",
    "code/pilot/36_aggregate_paired_radius_bgam_pilot.R",
    "code/pipeline/aggregate_results.R",
    "code/pipeline/compare_methods.R",
    "code/pipeline/fit_bgam.R",
    "code/pipeline/fit_paired_radius.R",
    "code/pipeline/generate_montecarlo.R",
    "code/pipeline/simulation.R",
    "code/reporting/05_slide_simulated_eeg_n20.R",
    "code/reporting/06_slide_bgam_fit_n20.R",
    "code/reporting/07_plot_paired_radius_truth_pilot.R",
    "code/reporting/07_slide_posterior_probability_n20.R",
    "code/reporting/08_slide_topomap_n20.R",
    "code/reporting/31_plot_all_electrode_sesoi_posterior.R",
    "code/reporting/37_plot_2sd_bgam_diagnostics.R",
    "code/sensitivity/10_prior_predictive_check.R",
    "code/sensitivity/11_basis_screening.R",
    "code/sensitivity/12_laplacian_basis_brms_fixed.R",
    "code/sensitivity/18_basis_screening_full_ab_sigma_true_control.R",
    "code/sensitivity/32_threshold_sensitivity.R",
    "code/sensitivity/33_alpha_sensitivity.R",
    "code/sensitivity/laplacian/00_config.R",
    "code/sensitivity/laplacian/01_fit_basis_scores.R",
    "code/sensitivity/laplacian/02_reconstruct_effect.R",
    "code/sensitivity/laplacian/03_diagnose_hard.R",
    "code/sensitivity/laplacian/04_diagnose_continuous.R",
    "code/truth/ground_truth.R",
    "code/truth/migrate_ground_truth.R",
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
