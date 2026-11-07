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
            bump_amplitude = 1.5,
            ronset_ms = seq(from = 150, to = 170, by = 2)
        ),
        low = list(
            omega_A = 0.20, omega_sigma = 0.00,
            bump_amplitude = 1.5,
            ronset_ms = seq(from = 140, to = 180, by = 2)
        ),
        medium = list(
            omega_A = 0.45, omega_sigma = 0.30,
            bump_amplitude = 1.5,
            ronset_ms = seq(from = 140, to = 180, by = 2)
        ),
        # Condition D (already generated/fitted as "high"): high absolute
        # heterogeneity (s ~= 1.193 uV), large amplitude.
        high = list(
            omega_A = 0.70, omega_sigma = 0.60,
            bump_amplitude = 1.5,
            ronset_ms = seq(from = 140, to = 180, by = 2)
        ),
        # Same as "high" but the noise-scale centring gives E(sigma_p^2) = 1.
        high_fixed = list(
            omega_A = 0.70, omega_sigma = 0.60,
            bump_amplitude = 1.5,
            ronset_ms = seq(from = 140, to = 180, by = 2),
            sigma_meanlog_mult = 1.0
        ),
        # 2x2 heterogeneity x amplitude design (Ladislas). omega_A is chosen
        # per condition so that the IMPLIED ABSOLUTE participant-amplitude SD
        # (s = bump_amplitude * sqrt(exp(omega_A^2) - 1)) lands on one of two
        # fixed tiers: s_low ~= 0.397 uV (cond_A, cond_B) or s_high ~= 1.193
        # uV (cond_C, matching "high"/cond_D). Cohen's d = bump_amplitude / s
        # is therefore independently tunable from the absolute heterogeneity
        # tier. cond_A and "high" share the same d by design (paired
        # contrast: same standardised effect, 3x the absolute heterogeneity
        # and 3x the amplitude).
        cond_A = list(
            omega_A = 0.700, omega_sigma = 0.60,
            bump_amplitude = 0.5,
            ronset_ms = seq(from = 140, to = 180, by = 2)
        ),
        cond_B = list(
            omega_A = 0.260, omega_sigma = 0.60,
            bump_amplitude = 1.5,
            ronset_ms = seq(from = 140, to = 180, by = 2)
        ),
        cond_C = list(
            omega_A = 1.379, omega_sigma = 0.60,
            bump_amplitude = 0.5,
            ronset_ms = seq(from = 140, to = 180, by = 2)
        ),
        # Noise-heterogeneity contrast: identical to "high"/cond_D (mu=1.5,
        # omega_A=0.700, s~=1.193) except omega_sigma=0, to isolate the
        # effect of noise-scale heterogeneity alone.
        cond_E = list(
            omega_A = 0.700, omega_sigma = 0.00,
            bump_amplitude = 1.5,
            ronset_ms = seq(from = 140, to = 180, by = 2)
        ),
        # Same noise-heterogeneity contrast (omega_sigma=0), paired with
        # cond_B, cond_C and cond_A respectively.
        cond_F = list(
            omega_A = 0.260, omega_sigma = 0.00,
            bump_amplitude = 1.5,
            ronset_ms = seq(from = 140, to = 180, by = 2)
        ),
        cond_G = list(
            omega_A = 0.26, omega_sigma = 0.00,
            bump_amplitude = 0.5,
            ronset_ms = seq(from = 140, to = 180, by = 2)
        ),
        cond_H = list(
            omega_A = 0.700, omega_sigma = 0.00,
            bump_amplitude = 0.5,
            ronset_ms = seq(from = 140, to = 180, by = 2)
        ),
        # New 2x2 design entirely at omega_sigma=0: cond_I (mu=0.5,
        # omega_A=0.26) pairs with cond_F (mu=1.5, omega_A=0.26), cond_H
        # (mu=0.5, omega_A=0.70) and cond_E (mu=1.5, omega_A=0.70).
        cond_I = list(
            omega_A = 0.260, omega_sigma = 0.00,
            bump_amplitude = 0.5,
            ronset_ms = seq(from = 140, to = 180, by = 2)
        ),
        # Same as cond_A/cond_B/cond_C but with the variance-centred
        # noise-scale convention (sigma_meanlog_mult = 1.0, matching
        # "high_fixed"): E(sigma_p^2) = 1 instead of E(sigma_p) = 1.
        cond_A_fixed = list(
            omega_A = 0.700, omega_sigma = 0.60,
            bump_amplitude = 0.5,
            ronset_ms = seq(from = 140, to = 180, by = 2),
            sigma_meanlog_mult = 1.0
        ),
        cond_B_fixed = list(
            omega_A = 0.260, omega_sigma = 0.60,
            bump_amplitude = 1.5,
            ronset_ms = seq(from = 140, to = 180, by = 2),
            sigma_meanlog_mult = 1.0
        ),
        cond_C_fixed = list(
            omega_A = 0.26, omega_sigma = 0.60,
            bump_amplitude = 0.5,
            ronset_ms = seq(from = 140, to = 180, by = 2),
            sigma_meanlog_mult = 1.0
        )
    )
}

define_simulation_parameters <- function(scenario = "legacy",
                                         bump_amplitude = NULL) {
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
        true_offset_sec = 0.340,
        downsample_factor = 5,
        seed_electrode = "CZ",
        sigma_space = 0.5,
        bump_amplitude = if (is.null(bump_amplitude)) {
            dispersion$bump_amplitude
        } else {
            bump_amplitude
        },
        rho = 0.15,
        signal_mode = "hard_signal",
        scenario = scenario,
        omega_A = dispersion$omega_A,
        omega_sigma = dispersion$omega_sigma,
        sigma_meanlog_mult = if (is.null(dispersion$sigma_meanlog_mult)) {
            0.5
        } else {
            dispersion$sigma_meanlog_mult
        }
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
                                  ),
                                  bump_amplitude = {
                                      env_val <- Sys.getenv(
                                          "MC_BUMP_AMPLITUDE", unset = NA
                                      )
                                      if (is.na(env_val)) NULL else as.numeric(env_val)
                                  }) {
    params <- define_simulation_parameters(scenario, bump_amplitude)
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

# Part A: raw Gaussian weights and deterministic hard-radius truncation.
radius_specifications <- function(sigma_space) {
    data.frame(
        radius_rule = c("1sd", "half-max", "2sd"),
        active_radius = c(
            sigma_space,
            sigma_space * sqrt(2 * log(2)),
            sigma_space * 2
        ),
        data_subdir = c(
            "montecarlo_r1sd",
            "montecarlo_rhalfmax",
            "montecarlo_r2sd"
        ),
        expected_active = c(5L, 9L, 21L),
        stringsAsFactors = FALSE
    )
}

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

# Downsample from 500 Hz to 100 Hz by keeping every fifth time point.
downsample_dataset <- function(sim_data, params) {
    all_times <- sort(unique(sim_data$time))
    keep_times <- all_times[seq(1, length(all_times),
                                by = params$downsample_factor)]

    sim_data %>%
        filter(time %in% keep_times) %>%
        mutate(across(where(is.numeric), ~round(.x, 4)))
}

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

make_paired_datasets <- function(seed = 1L, n_ppt = 20L,
                                 scenario = "high",
                                 bump_amplitude = NULL) {
    set.seed(seed)
    params <- simulation_parameters(
        n_ppt = n_ppt, scenario = scenario, signal_mode = "hard_signal",
        bump_amplitude = bump_amplitude
    )
    sensors <- compute_spatial_weight_base(load_electrode_layout(), params)
    noise_cov <- compute_noise_covariance(sensors, params)
    base <- generate_dataset_base(params, sensors, noise_cov$L)
    specs <- radius_specifications(params$sigma_space)

    results <- lapply(seq_len(nrow(specs)), function(i) {
        materialize_dataset(
            base, params, sensors,
            specs$active_radius[i], specs$radius_rule[i]
        )
    })
    names(results) <- specs$radius_rule
    results
}

make_dataset <- function(seed = 42, n_ppt = default_params$n_ppt,
                         scenario = Sys.getenv("MC_SCENARIO", "legacy"),
                         signal_mode = Sys.getenv("MC_SIGNAL_MODE",
                                                  "hard_signal"),
                         radius_rule = "1sd") {
    if (!identical(signal_mode, "hard_signal")) {
        stop("signal_mode is fixed to hard_signal", call. = FALSE)
    }
    paired <- make_paired_datasets(seed, n_ppt, scenario)
    if (!radius_rule %in% names(paired)) {
        stop("Unknown radius_rule: ", radius_rule, call. = FALSE)
    }
    paired[[radius_rule]]
}

# Load noise and ERP template before generating data.
load_noise_and_erp_assets()

# Run only when this file is executed directly.
if (!interactive() && sys.nframe() == 0) {
    pilot <- make_dataset(seed = 42)
    saveRDS(pilot$sim_data, data_path("sim_data_pilot.rds"))
    saveRDS(pilot$ground_truth, data_path("ground_truth_pilot.rds"))
}
