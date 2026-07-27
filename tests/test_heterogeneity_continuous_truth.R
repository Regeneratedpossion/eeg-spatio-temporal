script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[[1]])
project_root <- normalizePath(file.path(dirname(script_path), ".."))
Sys.setenv(PROJECT_ROOT = project_root)

source(file.path(project_root, "code", "00_generate_data.R"))

stopifnot(exists("heterogeneity_scenarios", mode = "function"))
scenarios <- heterogeneity_scenarios()
stopifnot(identical(names(scenarios), c("legacy", "low", "medium", "high")))
stopifnot(scenarios$medium$omega_A == 0.45)
stopifnot(scenarios$medium$omega_sigma == 0.30)

params <- simulation_parameters(n_ppt = 2L, scenario = "medium")
sensors <- data.frame(
    channel = c("CZ", "PZ"),
    spatial_weight_raw = c(1, 0.5),
    spatial_weight = c(1, 0.5),
    active_threshold = c(0.75, 0.75)
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
    0.5 * expected_seed_effect,
    tolerance = 0
)))

cat("heterogeneity continuous-truth regression tests passed\n")
