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
            true_offset = 0.342
        ),
        radii = list(
            radius_rule = c("1sd", "half-max", "2sd"),
            active_radius_multiplier = c(1, sqrt(2 * log(2)), 2),
            data_subdir = c(
                "montecarlo_r1sd",
                "montecarlo_rhalfmax",
                "montecarlo_r2sd"
            ),
            expected_active = c(5L, 9L, 21L)
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
            reference_offset = 0.342,
            timing_grid = seq(0.16, 0.34, by = 0.01)
        )
    )
}
