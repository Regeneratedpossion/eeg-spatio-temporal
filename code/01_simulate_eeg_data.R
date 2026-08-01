###########################################################
# 01_simulate_eeg_data.R
#
# Generate one simulated two-condition EEG dataset.
#
# Main structure:
#   1. simulation_parameters()
#   2. load_electrode_layout()
#   3. compute_spatial_weight_base() and apply_active_radius()
#   4. compute_noise_cholesky()
#   5. generate_dataset_base()
#   6. build_ground_truth() and materialize_dataset()
###########################################################

project_root_path <- Sys.getenv("PROJECT_ROOT", unset = getwd())
source(file.path(project_root_path, "code/utils.R"))
source(file.path(project_root_path, "code/paper_config.R"))

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

simulation_parameters <- function(n_ppt = 5L,
                                  scenario = Sys.getenv(
                                      "MC_SCENARIO",
                                      "legacy"
                                  )) {
    scenarios <- heterogeneity_scenarios()
    if (!scenario %in% names(scenarios)) {
        stop("Unknown heterogeneity scenario: ", scenario,
             ". Choose one of: ", paste(names(scenarios), collapse = ", "))
    }
    dispersion <- scenarios[[scenario]]

    list(
        n_trials = 20,
        n_ppt = n_ppt,
        outvar = 0.5,
        srate = 500,
        ronset = dispersion$ronset_ms,
        true_onset = 0.160,
        true_offset = 0.342,
        downsample_factor = 5,
        seed_electrode = "CZ",
        sigma_space = 0.5,
        bump_amplitude = 1.5,
        rho = 0.15,
        scenario = scenario,
        omega_A = dispersion$omega_A,
        omega_sigma = dispersion$omega_sigma
    )
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

# Read the fixed 61-channel layout without depending on generated data.
load_electrode_layout <- function() {
    read.csv(
        resource_path("electrodes", "electrode_layout_61ch.csv"),
        stringsAsFactors = FALSE
    )
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
    config <- paper_analysis_config()$radii
    data.frame(
        radius_rule = config$radius_rule,
        active_radius = sigma_space * config$active_radius_multiplier,
        data_subdir = config$data_subdir,
        expected_active = config$expected_active,
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
    valid_rules <- paper_analysis_config()$radii$radius_rule
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
compute_noise_cholesky <- function(sensors, params) {
    distance_matrix <- as.matrix(dist(sensors %>% select(xproj, yproj)))
    sigma <- exp(-distance_matrix^2 / (2 * params$rho^2))
    t(chol(sigma))
}

# -------------------------------------------------------
# 4. Generate one dataset
# -------------------------------------------------------

temporal_signal_at_onset <- function(onset_ms) {
    onset_index <- which(Xf == onset_ms)

    c(
        rep(0, onset_index - 2),
        erp,
        rep(0, Nf - onset_index - length(erp) + 2)
    )
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
    n_elec <- nrow(sensors)
    electrode_index <- rep(seq_len(n_elec), each = 2 * Nf)
    condition_time <- rep(c("cond1", "cond2"), each = Nf)
    signal_by_condition <- c(
        rep(0, Nf),
        participant_amplitude * temporal_signal
    )

    data.frame(
        participant = participant_id,
        channel = sensors$channel[electrode_index],
        xproj = sensors$xproj[electrode_index],
        yproj = sensors$yproj[electrode_index],
        condition = rep(condition_time, times = n_elec),
        trial = trial,
        time = rep(rep(Xf, 2), times = n_elec),
        eeg_noise = as.vector(rbind(noise_c1, noise_c2)),
        signal_unit = rep(signal_by_condition, times = n_elec)
    )
}

generate_dataset_base <- function(params, sensors, L) {
    n_elec <- nrow(sensors)
    trial_data <- list()
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
        onset_ms <- if (params$scenario == "legacy") {
            sample(x = params$ronset, size = 1)
        } else {
            participant$onset_ms
        }
        participant_parameters$onset_ms[P] <- onset_ms
        temporal_signal <- temporal_signal_at_onset(onset_ms)

        for (trial in seq_len(params$n_trials)) {
            noise_c1 <- participant$noise_scale *
                generate_condition_noise(n_elec, params, L)
            noise_c2 <- participant$noise_scale *
                generate_condition_noise(n_elec, params, L)

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

# Save ground-truth values for evaluation. The continuous field is the exact
# average of the participant waveforms injected into this generated dataset.
build_ground_truth <- function(params, sensors, participant_parameters) {
    active_radius <- sensors$active_radius[1]
    radius_rule <- sensors$radius_rule[1]
    generation_waveforms <- vapply(
        seq_len(nrow(participant_parameters)),
        function(i) {
            participant_parameters$amplitude[i] *
                temporal_signal_at_onset(participant_parameters$onset_ms[i])
        },
        numeric(length(Xf))
    )
    mean_effect_waveform <- rowMeans(generation_waveforms)

    keep_idx <- seq(1, length(Xf), by = params$downsample_factor)
    analysis_time <- Xf[keep_idx] / 1000
    analysis_effect_waveform <- mean_effect_waveform[keep_idx]
    # This fixed reference window is used for evaluation. It is not the exact
    # non-zero temporal support of every realised population field.
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
            effect_waveform = analysis_effect_waveform[
                match(time, analysis_time)
            ],
            temporal_weight = effect_waveform / params$bump_amplitude,
            effect = spatial_weight * effect_waveform
        ) %>%
        select(-effect_waveform)

    active_electrodes <- sensors %>%
        filter(spatial_weight > 0) %>%
        pull(channel)
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
        participant_onsets_ms = participant_parameters$onset_ms,
        temporal_active_times = temporal_active,
        population_mean_waveform = data.frame(
            time = analysis_time,
            temporal_weight = analysis_effect_waveform /
                params$bump_amplitude
        ),
        continuous_effect = continuous_ground_truth,
        active_cells = expand.grid(
            channel = active_electrodes,
            time = temporal_active,
            KEEP.OUT.ATTRS = FALSE,
            stringsAsFactors = FALSE
        ),
        participant_parameters = participant_parameters,
        all_channels = sensors$channel,
        active_electrodes = active_electrodes
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
                                 scenario = "high") {
    set.seed(seed)
    params <- simulation_parameters(n_ppt = n_ppt, scenario = scenario)
    sensors <- compute_spatial_weight_base(load_electrode_layout(), params)
    noise_cholesky <- compute_noise_cholesky(sensors, params)
    base <- generate_dataset_base(params, sensors, noise_cholesky)
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

make_dataset <- function(seed = 42, n_ppt = 5L,
                         scenario = Sys.getenv("MC_SCENARIO", "legacy"),
                         radius_rule = "1sd") {
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
