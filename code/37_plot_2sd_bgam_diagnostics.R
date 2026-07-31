diagnostic_channels <- function(available_channels) {
    requested <- c("CZ", "FCZ", "CPZ", "C2", "PZ", "FP1")
    missing <- setdiff(requested, available_channels)
    if (length(missing)) {
        stop(
            "Missing diagnostic channels: ",
            paste(missing, collapse = ", "),
            call. = FALSE
        )
    }
    requested
}

diagnostic_landmarks <- function(ground_truth, available_times) {
    if (is.null(ground_truth$true_onset) ||
        is.null(ground_truth$true_offset)) {
        stop("ground_truth must contain true_onset and true_offset",
             call. = FALSE)
    }
    targets <- c(
        ground_truth$true_onset,
        0.25,
        ground_truth$true_offset
    )
    unname(vapply(
        targets,
        function(target) {
            available_times[which.min(abs(available_times - target))]
        },
        numeric(1)
    ))
}

summarise_difference_draws <- function(
        diff_draws, grid, sesoi = 0.10, odds_threshold = 10) {
    diff_draws <- as.matrix(diff_draws)
    if (ncol(diff_draws) != nrow(grid)) {
        stop("Each grid row must correspond to one draw-matrix column",
             call. = FALSE)
    }
    if (length(sesoi) != 1L || !is.finite(sesoi)) {
        stop("sesoi must be one finite number", call. = FALSE)
    }
    if (length(odds_threshold) != 1L || !is.finite(odds_threshold) ||
        odds_threshold <= 0) {
        stop("odds_threshold must be one positive finite number",
             call. = FALSE)
    }
    probability_threshold <- odds_threshold / (1 + odds_threshold)
    output <- data.frame(
        grid,
        mean = colMeans(diff_draws),
        lo = unname(apply(diff_draws, 2, quantile, probs = 0.025)),
        hi = unname(apply(diff_draws, 2, quantile, probs = 0.975)),
        p_sig = colMeans(diff_draws > sesoi),
        stringsAsFactors = FALSE
    )
    output$above_threshold <- output$p_sig > probability_threshold
    attr(output, "sesoi") <- sesoi
    attr(output, "odds_threshold") <- odds_threshold
    attr(output, "probability_threshold") <- probability_threshold
    output
}
