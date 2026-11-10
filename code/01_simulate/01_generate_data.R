###########################################################
# 01_generate_data.R
#
# Functions that simulate one two-condition EEG dataset of the thesis (scenario cond_4, cond_3, cond_2 or cond_1).
# Run code/01_simulate/01_generate_dataset.R to write the datasets.
#
# Main structure:
#   1. heterogeneity_scenarios() and simulation_parameters()
#   2. load_electrode_layout()
#   3. compute_spatial_weight_base() and compute_noise_covariance()
#   4. generate_dataset_base() and materialize_dataset()
#   5. build_ground_truth()
#   6. make_dataset()
###########################################################

source(file.path(Sys.getenv("PROJECT_ROOT", unset = here::here()), "code/R/utils.R"))
source(file.path(project_root(), "code/R/superpopulation_truth.R"))

library(dplyr)

# Archived helper code for noise and ERP template.
BRMS_MEEG <- brms_meeg_path()

# -------------------------------------------------------
# 1. Simulation parameters
# -------------------------------------------------------

heterogeneity_scenarios <- function() {
    onsets <- seq(from = 140, to = 180, by = 2)
    list(
        # Scenario 1: weak effect, little variation between participants
        cond_1 = list(omega_A = 0.260, omega_sigma = 0.00, bump_amplitude = 0.5, ronset_ms = onsets),
        # Scenario 2: weak effect, large variation
        cond_2 = list(omega_A = 0.700, omega_sigma = 0.00, bump_amplitude = 0.5, ronset_ms = onsets),
        # Scenario 3: strong effect, little variation
        cond_3 = list(omega_A = 0.260, omega_sigma = 0.00, bump_amplitude = 1.5, ronset_ms = onsets),
        # Scenario 4: strong effect, large variation
        cond_4 = list(omega_A = 0.700, omega_sigma = 0.00, bump_amplitude = 1.5, ronset_ms = onsets)
    )
}

define_simulation_parameters <- function(scenario) {
    scenarios <- heterogeneity_scenarios()
    if (!scenario %in% names(scenarios)) {
        stop("Unknown heterogeneity scenario: ", scenario,
             ". Choose one of: ", paste(names(scenarios), collapse = ", "))
    }
    dispersion <- scenarios[[scenario]]

    list(
        n_trials = 20,
        n_ppt = 20,
        outvar = 0.5,
        srate = 500,
        ronset_ms = dispersion$ronset_ms,
        true_onset_sec = 0.160,
        true_offset_sec = 0.340,
        downsample_factor = 5,
        seed_electrode = "CZ",
        sigma_space = 0.5,
        bump_amplitude = dispersion$bump_amplitude,
        rho = 0.15,
        signal_mode = "hard_signal",
        scenario = scenario,
        omega_A = dispersion$omega_A,
        omega_sigma = dispersion$omega_sigma,
        sigma_meanlog_mult = 0.5
    )
}

simulation_parameters <- function(scenario, n_ppt = 20L, signal_mode = "hard_signal") {
    params <- define_simulation_parameters(scenario)
    if (!identical(signal_mode, "hard_signal")) {
        stop("signal_mode is fixed to hard_signal", call. = FALSE)
    }
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
            meanlog = -params$sigma_meanlog_mult * params$omega_sigma^2,
            sdlog = params$omega_sigma
        ),
        stringsAsFactors = FALSE
    )
}

# -------------------------------------------------------
# 2. Electrode layout and simulation assets
# -------------------------------------------------------

# Electrode coordinates, in the order that the noise matrix uses.
load_electrode_layout <- function() {
    read.csv(file.path(project_root(), "resources", "electrodes", "electrode_layout_61ch.csv"),
             stringsAsFactors = FALSE)[, c("channel", "xproj", "yproj")]
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

# Part A: raw Gaussian weights and deterministic hard-radius truncation at the half-maximum radius.
half_max_radius <- function(sigma_space) sigma_space * sqrt(2 * log(2))

compute_spatial_weight_base <- function(sensors, params) {
    seed_coord <- sensors %>%
        filter(channel == params$seed_electrode)
    if (nrow(seed_coord) != 1L) {
        stop("Expected exactly one seed electrode", call. = FALSE)
    }

    sensors %>%
        mutate(
            spatial_distance = sqrt((xproj - seed_coord$xproj)^2 +
                                    (yproj - seed_coord$yproj)^2),
            spatial_weight_raw = exp(-(spatial_distance^2) /
                                     (2 * params$sigma_space^2))
        )
}

apply_active_radius <- function(sensors, r_active, radius_rule) {
    valid_rules <- c("1sd", "half-max", "2sd")
    if (length(r_active) != 1L || !is.finite(r_active) || r_active <= 0) {
        stop("r_active must be one finite positive number", call. = FALSE)
    }
    if (length(radius_rule) != 1L || !radius_rule %in% valid_rules) {
        stop("Unknown radius_rule: ", radius_rule, call. = FALSE)
    }

    sensors %>%
        mutate(
            active_radius = r_active,
            radius_rule = radius_rule,
            spatial_weight = if_else(
                spatial_distance <= r_active,
                spatial_weight_raw,
                0
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

build_trial_base_rows <- function(participant_id, trial, sensors,
                                  temporal_signal, participant_amplitude,
                                  noise_c1, noise_c2) {
    rows <- lapply(seq_len(nrow(sensors)), function(e) {
        data.frame(
            participant = participant_id,
            channel = sensors$channel[e],
            xproj = sensors$xproj[e],
            yproj = sensors$yproj[e],
            condition = rep(c("cond1", "cond2"), each = Nf),
            trial = trial,
            time = rep(Xf, 2),
            eeg_noise = c(noise_c1[, e], noise_c2[, e]),
            signal_unit = c(
                rep(0, Nf),
                participant_amplitude * temporal_signal
            )
        )
    })

    bind_rows(rows)
}

generate_dataset_base <- function(params, sensors, L) {
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

            trial_data[[length(trial_data) + 1]] <- build_trial_base_rows(
                participant_id = participant_id,
                trial = trial,
                sensors = sensors,
                temporal_signal = temporal_signal,
                participant_amplitude = participant$amplitude,
                noise_c1 = noise_c1,
                noise_c2 = noise_c2
            )
        }
    }

    base_data <- bind_rows(trial_data) %>%
        mutate(
            time = time / 1000,
            condition = factor(condition, levels = c("cond1", "cond2"))
        )
    base_data <- downsample_base_dataset(base_data, params)

    list(
        base_data = base_data,
        participant_parameters = participant_parameters
    )
}

# -------------------------------------------------------
# 5. Ground truth
# -------------------------------------------------------

# Preserve the exact temporal signal used to construct continuous ground truth.
# Noise and public coordinates retain the historical four-decimal storage.
downsample_base_dataset <- function(base_data, params) {
    all_times <- sort(unique(base_data$time))
    keep_times <- all_times[seq(1, length(all_times),
                                by = params$downsample_factor)]
    numeric_columns <- names(base_data)[vapply(base_data, is.numeric, logical(1))]
    rounded_columns <- setdiff(numeric_columns, "signal_unit")

    base_data %>%
        filter(time %in% keep_times) %>%
        mutate(across(all_of(rounded_columns), ~round(.x, 4)))
}

# Save ground-truth values for evaluation.
# The continuous effect field is the theoretical population truth, not the seed's sample mean.
build_ground_truth <- function(params, sensors, participant_parameters) {
    active_radius <- unique(sensors$active_radius)
    radius_rule <- unique(sensors$radius_rule)
    stopifnot(length(active_radius) == 1L, length(radius_rule) == 1L)
    population_waveform <- population_effect_waveform(
        params$bump_amplitude,
        params$ronset,
        params$downsample_factor,
        file.path(BRMS_MEEG, "erp_template.R")
    )
    analysis_time <- population_waveform$time
    analysis_effect_waveform <- population_waveform$effect_waveform
    temporal_active <- analysis_time[
        analysis_time >= params$true_onset &
            analysis_time <= params$true_offset
    ]
    continuous_ground_truth <- expand.grid(
        channel = sensors$channel,
        time = analysis_time,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    ) %>%
        left_join(
            sensors %>% select(
                channel, spatial_distance, spatial_weight_raw, spatial_weight
            ),
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
        true_onset = params$true_onset,
        true_offset = params$true_offset,
        seed_electrode = params$seed_electrode,
        sigma_space = params$sigma_space,
        bump_amplitude = params$bump_amplitude,
        rho = params$rho,
        outvar = params$outvar,
        downsample_factor = params$downsample_factor,
        signal_mode = "hard_signal",
        active_radius = active_radius,
        radius_rule = radius_rule,
        scenario = params$scenario,
        onset_support_ms = params$ronset,
        omega_A = params$omega_A,
        omega_sigma = params$omega_sigma,
        participant_onsets_ms = participant_onsets_ms,
        temporal_active_times = temporal_active,
        population_mean_waveform = data.frame(
            time = analysis_time,
            temporal_weight = analysis_effect_waveform /
                params$bump_amplitude
        ),
        continuous_effect = continuous_ground_truth,
        active_cells = expand.grid(
            channel = sensors %>%
                filter(spatial_weight > 0) %>%
                pull(channel),
            time = temporal_active,
            KEEP.OUT.ATTRS = FALSE,
            stringsAsFactors = FALSE
        ),
        participant_parameters = participant_parameters,
        all_channels = sensors$channel,
        active_electrodes = sensors %>%
            filter(spatial_weight > 0) %>%
            pull(channel)
    )
}

materialize_dataset <- function(base, params, sensors, r_active, radius_rule) {
    radius_sensors <- apply_active_radius(sensors, r_active, radius_rule)
    weights <- radius_sensors$spatial_weight[
        match(base$base_data$channel, radius_sensors$channel)
    ]
    sim_data <- base$base_data %>%
        mutate(eeg = eeg_noise + weights * signal_unit) %>%
        select(
            participant, channel, xproj, yproj,
            condition, trial, time, eeg
        )
    ground_truth <- build_ground_truth(
        params,
        radius_sensors,
        participant_parameters = base$participant_parameters
    )

    list(sim_data = sim_data, ground_truth = ground_truth)
}

make_dataset <- function(seed, scenario, n_ppt = 20L) {
    set.seed(seed)
    params <- simulation_parameters(n_ppt = n_ppt, scenario = scenario)
    sensors <- compute_spatial_weight_base(load_electrode_layout(), params)
    noise_cov <- compute_noise_covariance(sensors, params)
    base <- generate_dataset_base(params, sensors, noise_cov$L)
    materialize_dataset(base, params, sensors, half_max_radius(params$sigma_space), "half-max")
}

# Load noise and ERP template before generating data.
load_noise_and_erp_assets()
