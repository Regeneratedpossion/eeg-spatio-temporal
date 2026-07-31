test_project_root <- Sys.getenv("PROJECT_ROOT", unset = getwd())

source(file.path(
    test_project_root, "code", "37_plot_2sd_bgam_diagnostics.R"
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

cat("2sd BGAM diagnostic figure tests passed\n")
