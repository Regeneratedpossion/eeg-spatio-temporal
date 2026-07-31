script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[[1]])
project_root <- normalizePath(file.path(dirname(script_path), ".."))
Sys.setenv(PROJECT_ROOT = project_root)

source(file.path(project_root, "code", "pipeline", "simulation.R"))

params <- simulation_parameters(n_ppt = 2L, scenario = "high")
stopifnot(identical(params$signal_mode, "hard_signal"))

layout <- load_electrode_layout()
sensor_base <- compute_spatial_weight_base(layout, params)
specs <- radius_specifications(params$sigma_space)

stopifnot(identical(specs$radius_rule, c("1sd", "half-max", "2sd")))
stopifnot(isTRUE(all.equal(
    specs$active_radius,
    c(0.5, 0.5 * sqrt(2 * log(2)), 1),
    tolerance = 0
)))

expected_active <- c(`1sd` = 5L, `half-max` = 9L, `2sd` = 21L)
for (i in seq_len(nrow(specs))) {
    sensors <- apply_active_radius(
        sensor_base,
        specs$active_radius[i],
        specs$radius_rule[i]
    )
    active <- sensors$spatial_weight > 0
    stopifnot(sum(active) == expected_active[[specs$radius_rule[i]]])
    stopifnot(all(sensors$spatial_weight[!active] == 0))
}

stopifnot(inherits(
    try(apply_active_radius(sensor_base, 0, "1sd"), silent = TRUE),
    "try-error"
))
stopifnot(inherits(
    try(apply_active_radius(sensor_base, 0.5, "invalid"), silent = TRUE),
    "try-error"
))

# Deterministic materialisation must neither draw random numbers nor mutate the
# common base. Dyadic fixture values make exact noise recovery meaningful.
params$downsample_factor <- 1L
precision_fixture <- data.frame(
    time = c(0, 0.01),
    eeg_noise = c(0.123456, 0.654321),
    signal_unit = c(0.123456789, 0.987654321)
)
precision_result <- downsample_base_dataset(precision_fixture, params)
stopifnot(identical(
    precision_result$signal_unit,
    precision_fixture$signal_unit
))
stopifnot(identical(
    precision_result$eeg_noise,
    round(precision_fixture$eeg_noise, 4)
))

sensors <- data.frame(
    channel = c("CZ", "PZ"),
    xproj = c(0, 1),
    yproj = c(0, 0),
    spatial_distance = c(0, 1),
    spatial_weight_raw = c(1, 0.5)
)
base_data <- expand.grid(
    participant = c("participant_01", "participant_02"),
    channel = c("CZ", "PZ"),
    condition = c("cond1", "cond2"),
    trial = 1L,
    time = c(0, 10),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
)
base_data$xproj <- ifelse(base_data$channel == "CZ", 0, 1)
base_data$yproj <- 0
base_data$condition <- factor(
    base_data$condition,
    levels = c("cond1", "cond2")
)
base_data$eeg_noise <- as.numeric(seq_len(nrow(base_data)))
base_data$signal_unit <- ifelse(
    base_data$condition == "cond2",
    ifelse(base_data$time == 0, 2, 4),
    0
)
participants <- data.frame(
    participant = c("participant_01", "participant_02"),
    onset_ms = c(140, 180),
    amplitude = c(1, 2),
    noise_scale = c(1, 1),
    stringsAsFactors = FALSE
)
base <- list(
    base_data = base_data,
    participant_parameters = participants
)
base_copy <- base

set.seed(123)
rng_before <- .Random.seed
narrow <- materialize_dataset(base, params, sensors, 0.5, "1sd")
wide <- materialize_dataset(base, params, sensors, 1, "2sd")
stopifnot(identical(.Random.seed, rng_before))
stopifnot(identical(base, base_copy))

recover_noise <- function(result, radius_sensors) {
    weight <- radius_sensors$spatial_weight[
        match(result$sim_data$channel, radius_sensors$channel)
    ]
    result$sim_data$eeg - weight * base_data$signal_unit
}
stopifnot(identical(
    recover_noise(narrow, apply_active_radius(sensors, 0.5, "1sd")),
    base_data$eeg_noise
))
stopifnot(identical(
    recover_noise(wide, apply_active_radius(sensors, 1, "2sd")),
    base_data$eeg_noise
))

for (result in list(narrow, wide)) {
    gt <- result$ground_truth
    inactive <- !(gt$continuous_effect$channel %in% gt$active_electrodes)
    if (any(inactive)) {
        stopifnot(max(abs(gt$continuous_effect$effect[inactive])) == 0)
        stopifnot(sum(gt$continuous_effect$effect[inactive] != 0) == 0L)
    }
    stopifnot(identical(gt$signal_mode, "hard_signal"))
}

cat("paired hard-radius generator tests passed\n")
