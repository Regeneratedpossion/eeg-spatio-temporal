###########################################################
# 00_generate_data.R
#
# Generate one simulated two-condition EEG dataset.
#
# Main structure:
#   1. define_simulation_parameters()
#   2. load_electrode_layout()
#   3. compute_spatial_weights()
#   4. compute_noise_covariance()
#   5. generate_one_dataset()
#   6. downsample_dataset()
#   7. build_ground_truth()
###########################################################

source(file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()), "code/utils.R"))

library(dplyr)

# Archived helper code for noise and ERP template.
BRMS_MEEG <- brms_meeg_path()

# -------------------------------------------------------
# 1. Simulation parameters
# -------------------------------------------------------

heterogeneity_scenarios <- function() {
    list(
        legacy = list(
            omega_A = 0.00, omega_sigma = 0.00,
            ronset_ms = seq(from = 150, to = 170, by = 2)
        ),
        low = list(
            omega_A = 0.20, omega_sigma = 0.00,
            ronset_ms = seq(from = 140, to = 180, by = 2)
        ),
        medium = list(
            omega_A = 0.45, omega_sigma = 0.30,
            ronset_ms = seq(from = 140, to = 180, by = 2)
        ),
        high = list(
            omega_A = 0.70, omega_sigma = 0.60,
            ronset_ms = seq(from = 140, to = 180, by = 2)
        )
    )
}

define_simulation_parameters <- function(scenario = "legacy") {
    scenarios <- heterogeneity_scenarios()
    if (!scenario %in% names(scenarios)) {
        stop("Unknown heterogeneity scenario: ", scenario,
             ". Choose one of: ", paste(names(scenarios), collapse = ", "))
    }
    dispersion <- scenarios[[scenario]]

    list(
        n_trials = 20,
        n_ppt = 5,
        outvar = 0.5,
        srate = 500,
        ronset_ms = dispersion$ronset_ms,
        true_onset_sec = 0.160,
        true_offset_sec = 0.342,
        downsample_factor = 5,
        seed_electrode = "CZ",
        sigma_space = 0.5,
        bump_amplitude = 1.5,
        rho = 0.15,
        signal_mode = "continuous_label",
        scenario = scenario,
        omega_A = dispersion$omega_A,
        omega_sigma = dispersion$omega_sigma
    )
}

default_params <- define_simulation_parameters()

simulation_parameters <- function(n_ppt = default_params$n_ppt,
                                  signal_mode = Sys.getenv(
                                      "MC_SIGNAL_MODE",
                                      default_params$signal_mode
                                  ),
                                  scenario = Sys.getenv(
                                      "MC_SCENARIO",
                                      "legacy"
                                  )) {
    params <- define_simulation_parameters(scenario)
    params$n_ppt <- n_ppt
    params$signal_mode <- signal_mode
    params$ronset <- params$ronset_ms
    params$true_onset <- params$true_onset_sec
    params$true_offset <- params$true_offset_sec
    params$ronset_ms <- NULL
    params$true_onset_sec <- NULL
    params$true_offset_sec <- NULL
    params
}

draw_participant_parameters <- function(params) {
    data.frame(
        participant = paste0(
            "participant_", sprintf("%02d", seq_len(params$n_ppt))
        ),
        onset_ms = sample(params$ronset, params$n_ppt, replace = TRUE),
        amplitude = rlnorm(
            params$n_ppt,
            meanlog = log(params$bump_amplitude) - 0.5 * params$omega_A^2,
            sdlog = params$omega_A
        ),
        noise_scale = rlnorm(
            params$n_ppt,
            meanlog = -params$omega_sigma^2,
            sdlog = params$omega_sigma
        ),
        stringsAsFactors = FALSE
    )
}

# -------------------------------------------------------
# 2. Electrode layout and simulation assets
# -------------------------------------------------------

# Reuse the electrode coordinates from the pilot data.
load_electrode_layout <- function() {
    sim_pilot <- readRDS(data_path("sim_data_pilot.rds"))
    sim_pilot %>%
        distinct(channel, xproj, yproj)
}

# Load the external noise generator, ERP template, and frequency-power profile.
load_noise_and_erp_assets <- function() {
    source(file.path(BRMS_MEEG, "eeg_noise.R"))
    source(file.path(BRMS_MEEG, "erp_template.R"))
    meanpower <<- unlist(read.table(file.path(BRMS_MEEG, "meanpower.txt")))
    invisible(meanpower)
}

# -------------------------------------------------------
# 3. Spatial weights and noise covariance
# -------------------------------------------------------

# Part A: spatial weights for the true signal.
# `continuous_label`: continuous Gaussian signal; threshold only labels active electrodes.
# `hard_signal`: thresholded Gaussian signal; outside-threshold weights are set to 0.

compute_spatial_weights <- function(sensors, params) {
    stopifnot(params$signal_mode %in% c("continuous_label", "hard_signal"))

    seed_coord <- sensors %>%
        filter(channel == params$seed_electrode)

    sensors <- sensors %>%
        mutate(
            spatial_distance = sqrt((xproj - seed_coord$xproj)^2 +
                                    (yproj - seed_coord$yproj)^2),
            spatial_weight_raw = exp(-(spatial_distance^2) /
                                     (2 * params$sigma_space^2))
        )

    active_threshold <- mean(sensors$spatial_weight_raw) +
        sd(sensors$spatial_weight_raw)

    sensors %>%
        mutate(
            active_threshold = active_threshold,
            signal_mode = params$signal_mode,
            spatial_weight = if_else(
                params$signal_mode == "hard_signal" &
                    spatial_weight_raw < active_threshold,
                0,
                spatial_weight_raw
            )
        )
}

# Part B: spatial covariance for the noise.
# rho controls how similar the noise is between nearby electrodes.
# L is used later to generate spatially correlated noise.
compute_noise_covariance <- function(sensors, params) {
    distance_matrix <- as.matrix(dist(sensors %>% select(xproj, yproj)))
    sigma <- exp(-distance_matrix^2 / (2 * params$rho^2))
    rownames(sigma) <- colnames(sigma) <- sensors$channel

    list(
        Sigma = sigma,
        L = t(chol(sigma))
    )
}

# -------------------------------------------------------
# 4. Generate one dataset
# -------------------------------------------------------

temporal_signal_at_onset <- function(params, onset_ms) {
    onset_index <- which(Xf == onset_ms)

    c(
        rep(0, onset_index - 2),
        erp,
        rep(0, Nf - onset_index - length(erp) + 2)
    )
}

build_temporal_signal <- function(params, onset_ms = NULL) {
    participant_onset <- if (is.null(onset_ms)) {
        sample(x = params$ronset, size = 1)
    } else {
        onset_ms
    }
    signal <- temporal_signal_at_onset(params, participant_onset)
    attr(signal, "onset_ms") <- participant_onset
    signal
}

generate_condition_noise <- function(n_elec, params, L) {
    temporal_noise <- sapply(seq_len(n_elec), function(e)
        eeg_noise(
            frames = Nf,
            srate = params$srate,
            outvar = params$outvar,
            meanpower = meanpower
        )
    )

    spatial_noise <- matrix(
        rnorm(Nf * n_elec, sd = sqrt(params$outvar)),
        nrow = Nf
    ) %*% t(L)

    temporal_noise + spatial_noise
}

build_trial_rows <- function(participant_id, trial, sensors,
                             temporal_signal, participant_amplitude,
                             noise_c1, noise_c2, params) {
    rows <- lapply(seq_len(nrow(sensors)), function(e) {
        signal_c1 <- noise_c1[, e]
        signal_c2 <- participant_amplitude *
            sensors$spatial_weight[e] * temporal_signal +
            noise_c2[, e]

        data.frame(
            participant = participant_id,
            channel = sensors$channel[e],
            xproj = sensors$xproj[e],
            yproj = sensors$yproj[e],
            condition = rep(c("cond1", "cond2"), each = Nf),
            trial = trial,
            time = rep(Xf, 2),
            eeg = c(signal_c1, signal_c2)
        )
    })

    bind_rows(rows)
}

generate_one_dataset <- function(params, sensors, L) {
    n_elec <- nrow(sensors)
    trial_data <- list()
    participant_onsets_ms <- numeric(params$n_ppt)
    participant_parameters <- if (params$scenario == "legacy") {
        data.frame(
            participant = paste0(
                "participant_", sprintf("%02d", seq_len(params$n_ppt))
            ),
            onset_ms = NA_real_,
            amplitude = params$bump_amplitude,
            noise_scale = 1,
            stringsAsFactors = FALSE
        )
    } else {
        draw_participant_parameters(params)
    }

    for (P in seq_len(params$n_ppt)) {
        participant <- participant_parameters[P, ]
        participant_id <- participant$participant
        # The legacy branch samples here, preserving the historical RNG order.
        temporal_signal <- if (params$scenario == "legacy") {
            build_temporal_signal(params)
        } else {
            build_temporal_signal(params, participant$onset_ms)
        }
        participant_onsets_ms[P] <- attr(temporal_signal, "onset_ms")
        participant_parameters$onset_ms[P] <- participant_onsets_ms[P]

        for (trial in seq_len(params$n_trials)) {
            if (params$scenario == "legacy") {
                noise_c1 <- generate_condition_noise(n_elec, params, L)
                noise_c2 <- generate_condition_noise(n_elec, params, L)
            } else {
                noise_c1 <- participant$noise_scale *
                    generate_condition_noise(n_elec, params, L)
                noise_c2 <- participant$noise_scale *
                    generate_condition_noise(n_elec, params, L)
            }

            trial_data[[length(trial_data) + 1]] <- build_trial_rows(
                participant_id = participant_id,
                trial = trial,
                sensors = sensors,
                temporal_signal = temporal_signal,
                participant_amplitude = participant$amplitude,
                noise_c1 = noise_c1,
                noise_c2 = noise_c2,
                params = params
            )
        }
    }

    sim_data <- bind_rows(trial_data) %>%
        mutate(
            time = time / 1000,
            condition = factor(condition, levels = c("cond1", "cond2"))
        )

    attr(sim_data, "participant_onsets_ms") <- participant_onsets_ms
    attr(sim_data, "participant_parameters") <- participant_parameters
    sim_data
}

# -------------------------------------------------------
# 5. Ground truth
# -------------------------------------------------------

# Downsample from 500 Hz to 100 Hz by keeping every fifth time point.
downsample_dataset <- function(sim_data, params) {
    all_times <- sort(unique(sim_data$time))
    keep_times <- all_times[seq(1, length(all_times),
                                by = params$downsample_factor)]

    sim_data %>%
        filter(time %in% keep_times) %>%
        mutate(across(where(is.numeric), ~round(.x, 4)))
}

# Save ground-truth values for evaluation. The continuous field is the exact
# average of the participant waveforms injected into this generated dataset.
build_ground_truth <- function(params, sensors, participant_parameters) {
    active_threshold <- unique(sensors$active_threshold)
    generation_waveforms <- vapply(
        seq_len(nrow(participant_parameters)),
        function(i) {
            participant_parameters$amplitude[i] *
                temporal_signal_at_onset(
                    params,
                    participant_parameters$onset_ms[i]
                )
        },
        numeric(length(Xf))
    )
    mean_effect_waveform <- rowMeans(generation_waveforms)

    reference_waveforms <- vapply(
        params$ronset,
        function(onset_ms) temporal_signal_at_onset(params, onset_ms),
        numeric(length(Xf))
    )
    reference_waveform <- rowMeans(reference_waveforms)

    keep_idx <- seq(1, length(Xf), by = params$downsample_factor)
    analysis_time <- Xf[keep_idx] / 1000
    analysis_effect_waveform <- mean_effect_waveform[keep_idx]
    analysis_reference_waveform <- reference_waveform[keep_idx]
    temporal_threshold <- mean(analysis_reference_waveform) +
        sd(analysis_reference_waveform)
    temporal_active <- analysis_time[
        analysis_reference_waveform >= temporal_threshold
    ]
    continuous_ground_truth <- expand.grid(
        channel = sensors$channel,
        time = analysis_time,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    ) %>%
        left_join(
            sensors %>% select(channel, spatial_weight),
            by = "channel"
        ) %>%
        mutate(
            temporal_weight = analysis_effect_waveform[
                match(time, analysis_time)
            ] / params$bump_amplitude,
            effect = spatial_weight * analysis_effect_waveform[
                match(time, analysis_time)
            ]
        )

    participant_onsets_ms <- participant_parameters$onset_ms
    list(
        true_onset = if (params$scenario == "legacy") {
            params$true_onset
        } else {
            min(temporal_active)
        },
        true_offset = if (params$scenario == "legacy") {
            params$true_offset
        } else {
            max(temporal_active)
        },
        seed_electrode = params$seed_electrode,
        sigma_space = params$sigma_space,
        bump_amplitude = params$bump_amplitude,
        rho = params$rho,
        outvar = params$outvar,
        downsample_factor = params$downsample_factor,
        signal_mode = params$signal_mode,
        scenario = params$scenario,
        onset_support_ms = params$ronset,
        omega_A = params$omega_A,
        omega_sigma = params$omega_sigma,
        active_threshold = active_threshold,
        participant_onsets_ms = participant_onsets_ms,
        temporal_threshold = temporal_threshold,
        temporal_active_times = temporal_active,
        population_mean_waveform = data.frame(
            time = analysis_time,
            temporal_weight = analysis_effect_waveform /
                params$bump_amplitude
        ),
        continuous_effect = continuous_ground_truth,
        active_cells = expand.grid(
            channel = sensors %>%
                filter(spatial_weight_raw >= active_threshold) %>%
                pull(channel),
            time = temporal_active,
            KEEP.OUT.ATTRS = FALSE,
            stringsAsFactors = FALSE
        ),
        participant_parameters = participant_parameters,
        all_channels = sensors$channel,
        active_electrodes = sensors %>%
            filter(spatial_weight_raw >= active_threshold) %>%
            pull(channel)
    )
}

make_dataset <- function(seed = 42, n_ppt = default_params$n_ppt,
                         scenario = Sys.getenv("MC_SCENARIO", "legacy"),
                         signal_mode = Sys.getenv("MC_SIGNAL_MODE",
                                                  "continuous_label")) {
    set.seed(seed)

    params <- simulation_parameters(
        n_ppt = n_ppt,
        scenario = scenario,
        signal_mode = signal_mode
    )
    sensors <- compute_spatial_weights(load_electrode_layout(), params)
    noise_cov <- compute_noise_covariance(sensors, params)

    sim_data <- generate_one_dataset(params, sensors, noise_cov$L)
    participant_parameters <- attr(sim_data, "participant_parameters")
    sim_data <- downsample_dataset(sim_data, params)
    ground_truth <- build_ground_truth(
        params,
        sensors,
        participant_parameters = participant_parameters
    )

    list(
        sim_data = sim_data,
        ground_truth = ground_truth
    )
}

# Load noise and ERP template before generating data.
load_noise_and_erp_assets()

# Run only when this file is executed directly.
if (!interactive() && sys.nframe() == 0) {
    pilot <- make_dataset(seed = 42)
    saveRDS(pilot$sim_data, data_path("sim_data_pilot.rds"))
    saveRDS(pilot$ground_truth, data_path("ground_truth_pilot.rds"))
}
