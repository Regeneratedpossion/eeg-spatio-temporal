# Formal configuration for the reproducible paired-radius paper analysis.
# This module intentionally uses base R only and has no source-time effects.
paper_analysis_config <- function() {
    list(
        simulation = list(
            scenario = "high",
            n_ppt = 20L,
            n_trials = 20L,
            n_reps = 20L,
            signal_mode = "hard_signal",
            true_onset = 0.160,
            true_offset = 0.340
        ),
        radii = list(
            radius_rule = "half-max",
            active_radius_multiplier = sqrt(2 * log(2)),
            data_subdir = "montecarlo_rhalfmax",
            expected_active = 9L
        ),
        model = list(
            k_xy = 16L,
            k_t = 12L,
            bs_space = "ps",
            bs_time = "ps",
            full = FALSE,
            chains = 4L,
            iter = 3000L,
            warmup = 1000L,
            adapt_delta = 0.80,
            max_treedepth = 10L,
            se_variable = "eeg_sem",
            sigma_true = TRUE,
            prior_mode = "brms_default",
            seed = 1L
        ),
        evaluation = list(
            reference_onset = 0.160,
            reference_offset = 0.340,
            timing_grid = seq(0.16, 0.34, by = 0.01)
        )
    )
}
