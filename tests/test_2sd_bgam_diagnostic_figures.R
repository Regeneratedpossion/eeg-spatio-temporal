test_project_root <- Sys.getenv("PROJECT_ROOT", unset = getwd())

source(file.path(
    test_project_root, "code", "17_plot_2sd_bgam_diagnostics.R"
))

available_channels <- c("AF1", "C2", "CPZ", "CZ", "FCZ", "FP1", "PZ")
stopifnot(identical(
    diagnostic_channels(available_channels),
    c("CZ", "FCZ", "CPZ", "C2", "PZ", "FP1")
))

ground_truth_fixture <- list(true_onset = 0.20, true_offset = 0.30)
available_times <- seq(0, 0.50, by = 0.01)
stopifnot(identical(
    diagnostic_landmarks(ground_truth_fixture, available_times),
    c(0.20, 0.25, 0.30)
))
stopifnot(identical(
    topomap_diagnostic_times(available_times),
    c(0.16, 0.25, 0.34)
))
stopifnot(topomap_colour_limit(data.frame(value = c(-2, 0.5, 1.5))) == 2)
contour_fixture <- data.frame(
    time_label = rep(c("narrow", "variable"), each = 2),
    value = c(0, 0.1, 0, 1)
)
contour_filtered <- variable_contour_data(
    contour_fixture,
    breaks = c(-1, -0.5, 0, 0.5, 1)
)
stopifnot(identical(unique(contour_filtered$time_label), "variable"))

grid_fixture <- data.frame(
    channel = c("CZ", "PZ"),
    time = c(0.25, 0.25),
    stringsAsFactors = FALSE
)
draw_fixture <- cbind(
    c(0.20, 0.15, 0.05, 0.11),
    c(0.20, 0.20, 0.20, 0.20)
)
summary_fixture <- summarise_difference_draws(
    draw_fixture,
    grid_fixture,
    sesoi = 0.10,
    odds_threshold = 10
)
stopifnot(
    isTRUE(all.equal(summary_fixture$p_sig, c(0.75, 1), tolerance = 0)),
    identical(summary_fixture$above_threshold, c(FALSE, TRUE)),
    isTRUE(all.equal(
        summary_fixture$mean,
        colMeans(draw_fixture),
        tolerance = 0
    )),
    isTRUE(all.equal(
        summary_fixture$lo,
        apply(draw_fixture, 2, quantile, probs = 0.025),
        tolerance = 0
    )),
    isTRUE(all.equal(
        summary_fixture$hi,
        apply(draw_fixture, 2, quantile, probs = 0.975),
        tolerance = 0
    )),
    isTRUE(all.equal(
        attr(summary_fixture, "probability_threshold"),
        10 / 11,
        tolerance = 0
    ))
)

plot_source <- paste(
    readLines(file.path(
        test_project_root, "code", "17_plot_2sd_bgam_diagnostics.R"
    )),
    collapse = "\n"
)
for (required_text in c(
    "posterior_epred(",
    "re_formula = NA",
    "ROPE_UV",
    "paired_radius_superpopulation_truth(",
    '"bgam_fit_superpopulation.png"',
    '"posterior_probability_sesoi.png"',
    '"topomap_diagnostics_superpopulation.png"',
    '"diagnostic_plot_data.rds"'
)) {
    stopifnot(grepl(required_text, plot_source, fixed = TRUE))
}
for (required_topomap_text in c(
    "topomap_times <- topomap_diagnostic_times(times)",
    "shared_limit <- topomap_colour_limit(interpolated)",
    '"1. Ground truth (superpopulation)"',
    '"2. Empirical difference"',
    '"3. Noise residual (Empirical - Ground truth)"',
    '"4. BGAM estimate"',
    '"5. Empirical - BGAM"',
    "topomap_times = topomap_times",
    "topomap_colour_limit = shared_limit"
)) {
    stopifnot(grepl(required_topomap_text, plot_source, fixed = TRUE))
}

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
})
sim_model_fixture <- expand.grid(
    channel = c("A", "B"),
    condition = factor(c("cond1", "cond2"), levels = c("cond1", "cond2")),
    time = 0.25,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
) %>%
    mutate(
        xproj = if_else(channel == "A", 0, 1),
        yproj = 0,
        eeg_mean = if_else(condition == "cond2", 1, 0)
    )
posterior_map_fixture <- data.frame(
    channel = c("A", "B"), time = 0.25,
    xproj = c(0, 1), yproj = c(0, 0),
    mean_diff = c(0.9, 0.1), p_sig = c(1, 0), odds = c(8000, 0.000125)
)
population_truth_fixture <- data.frame(
    channel = c("A", "B"), time = 0.25, effect = c(1, 0)
)
topomap_fixture <- prepare_topomap_data(
    sim_model_fixture,
    posterior_map_fixture,
    population_truth_fixture,
    list(active_electrodes = "A"),
    landmarks = 0.25
)
stopifnot(
    nrow(topomap_fixture) == 10L,
    all(c("xproj", "yproj") %in% names(topomap_fixture)),
    !any(grepl("\\.[xy]$", names(topomap_fixture)))
)

cat("2sd BGAM diagnostic figure tests passed\n")
