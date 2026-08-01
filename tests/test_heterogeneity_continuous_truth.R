script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[[1]])
project_root <- normalizePath(file.path(dirname(script_path), ".."))
Sys.setenv(PROJECT_ROOT = project_root)

source(file.path(project_root, "code", "01_simulate_eeg_data.R"))

stopifnot(exists("heterogeneity_scenarios", mode = "function"))
scenarios <- heterogeneity_scenarios()
stopifnot(identical(names(scenarios), c("legacy", "low", "medium", "high")))
stopifnot(scenarios$medium$omega_A == 0.45)
stopifnot(scenarios$medium$omega_sigma == 0.30)

params <- simulation_parameters(n_ppt = 2L, scenario = "medium")
sensors <- data.frame(
    channel = c("CZ", "PZ"),
    spatial_distance = c(0, 1),
    spatial_weight_raw = c(1, 0.5),
    spatial_weight = c(1, 0),
    active_radius = c(0.75, 0.75),
    radius_rule = c("1sd", "1sd")
)
participants <- data.frame(
    participant = c("participant_01", "participant_02"),
    onset_ms = c(140, 180),
    amplitude = c(1, 3),
    noise_scale = c(0.7, 1.4)
)

truth <- build_ground_truth(params, sensors, participants)
keep_idx <- seq(1, length(Xf), by = params$downsample_factor)
expected_seed_effect <- (
    participants$amplitude[1] *
        temporal_signal_at_onset(params, participants$onset_ms[1]) +
    participants$amplitude[2] *
        temporal_signal_at_onset(params, participants$onset_ms[2])
) / nrow(participants)
expected_seed_effect <- expected_seed_effect[keep_idx]

actual_seed_effect <- subset(
    truth$continuous_effect,
    channel == "CZ"
)$effect
stopifnot(isTRUE(all.equal(actual_seed_effect, expected_seed_effect,
                           tolerance = 0)))
stopifnot(isTRUE(all.equal(
    subset(truth$continuous_effect, channel == "PZ")$effect,
    rep(0, length(expected_seed_effect)),
    tolerance = 0
)))
stopifnot(identical(truth$active_electrodes, "CZ"))
stopifnot(identical(truth$active_radius, 0.75))
stopifnot(identical(truth$radius_rule, "1sd"))
stopifnot(identical(truth$signal_mode, "hard_signal"))
stopifnot(identical(truth$true_onset, 0.160))
stopifnot(identical(truth$true_offset, 0.342))
expected_active_times <- seq(0.16, 0.34, by = 0.01)
stopifnot(isTRUE(all.equal(
    truth$temporal_active_times,
    expected_active_times,
    tolerance = sqrt(.Machine$double.eps)
)))
stopifnot(is.null(truth$temporal_threshold))
expected_active_cells <- expand.grid(
    channel = "CZ",
    time = expected_active_times,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
)
stopifnot(identical(truth$active_cells$channel,
                    expected_active_cells$channel))
stopifnot(isTRUE(all.equal(
    truth$active_cells$time,
    expected_active_cells$time,
    tolerance = sqrt(.Machine$double.eps)
)))

cat("heterogeneity continuous-truth regression tests passed\n")
