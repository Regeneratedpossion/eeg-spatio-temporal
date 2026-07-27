###########################################################
# Hard-threshold diagnostics.
#
# This keeps the original active/inactive evaluation:
# active channel x time cells are compared against posterior odds.
###########################################################

temporal_recovery_pass <- function(times, seed_odds, ground_truth) {
    clusters <- find_clusters(times, seed_odds, ODDS_THRESHOLD)
    estimated <- c(clusters$cluster_onset[1], clusters$cluster_offset[1])
    truth <- c(ground_truth$true_onset, ground_truth$true_offset)

    error_ms <- abs(estimated - truth) * 1000
    error_ms[is.na(error_ms)] <- Inf
    max(error_ms) <= BIAS_TOLERANCE_MS
}

posterior_cell_summaries <- function(effect_array) {
    prob_above_rope <- apply(effect_array > ROPE_UV, c(2, 3), mean)
    min_prob <- 1 / dim(effect_array)[1]

    list(
        mean = apply(effect_array, c(2, 3), mean),
        lower = apply(effect_array, c(2, 3), quantile, probs = 0.025),
        odds = pmax(prob_above_rope, min_prob) /
            pmax(1 - prob_above_rope, min_prob)
    )
}

mask_f1 <- function(predicted, truth) {
    TP <- sum(predicted & truth)
    FP <- sum(predicted & !truth)
    FN <- sum(!predicted & truth)

    2 * TP / (2 * TP + FP + FN)
}

convergence_pass <- function(summary) {
    max(summary$n_divergent, na.rm = TRUE) == 0 &&
        max(summary$n_treedepth, na.rm = TRUE) == 0 &&
        max(summary$max_rhat, na.rm = TRUE) < 1.01 &&
        min(summary$min_ess_bulk, na.rm = TRUE) > 400
}

recovery_summary <- function(draws, truth) {
    ci <- quantile(draws, c(0.025, 0.975))

    list(
        ci_covers_truth = ci[1] <= truth && truth <= ci[2],
        recovery_pct = 100 * mean(draws) / truth
    )
}

run_hard_diagnostics <- function() {
    # 1. Load reconstructed posterior effect map and ground truth.
    rec <- load_laplace_reconstruction()
    out_dir <- file.path(project_root(), "results", "diagnostics",
                         "laplacian_basis_fixed")
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

    effect_array <- rec$effect_array
    times <- rec$times
    channels <- rec$channels
    gt <- rec$gt
    seed_idx <- match(gt$seed_electrode, channels)
    active_channels <- channels %in% gt$active_electrodes

    # 2. Summarise the posterior effect map at every time x channel cell.
    # Posterior probability: mean(draws > ROPE_UV).
    posterior <- posterior_cell_summaries(effect_array)

    # 3. Temporal recovery:
    # does the detected seed-channel cluster recover the true onset/offset?
    temporal_pass <- temporal_recovery_pass(times, posterior$odds[, seed_idx], gt)

    # 4. Spatial recovery:
    # at the true peak time, compare predicted active channels with true channels.
    seed_truth <- continuous_truth_cells(
        gt,
        channels = gt$seed_electrode,
        times = times
    )
    peak_time <- seed_truth$time[which.max(seed_truth$effect)]
    peak_idx <- match(peak_time, times)
    pred_active_space <- posterior$lower[peak_idx, ] > 0
    spatial_F1 <- mask_f1(pred_active_space, active_channels)

    # 5. Spatiotemporal recovery:
    # compare predicted and true active masks over all time x channel cells.
    pred_active <- posterior$odds > ODDS_THRESHOLD
    truth_active <- outer(times >= gt$true_onset & times <= gt$true_offset,
                          active_channels, `&`)
    st_F1 <- mask_f1(pred_active, truth_active)

    # 6. Amplitude recovery:
    # compare posterior peak/window amplitude with the known simulated truth.
    seed_draws <- effect_array[, , seed_idx]
    seed_peak_draws <- seed_draws[, peak_idx]
    window_idx <- which(times >= gt$true_onset & times <= gt$true_offset)
    window_draws <- rowMeans(seed_draws[, window_idx, drop = FALSE])
    truth <- continuous_truth_summary(
        gt,
        channel = gt$seed_electrode,
        peak_time = peak_time,
        window_times = times[window_idx]
    )

    # 7. Convergence, top-channel, and amplitude checks.
    spatial_top_pass <- channels[which.max(posterior$mean[peak_idx, ])] %in%
        cz_neighbours
    peak_recovery <- recovery_summary(seed_peak_draws, truth$peak)
    mean_recovery <- recovery_summary(window_draws, truth$window_mean)

    # 8. Save one-row hard-threshold diagnostic summary.
    hard_summary <- data.frame(
        fit_file = rec$tag,
        dataset = data_subdir,
        convergence_pass = convergence_pass(rec$conv_summary),
        temporal_pass = temporal_pass,
        spatial_top_pass = spatial_top_pass,
        peak_ci_covers_truth = peak_recovery$ci_covers_truth,
        mean_ci_covers_truth = mean_recovery$ci_covers_truth,
        spatial_F1 = spatial_F1,
        spatiotemporal_F1 = st_F1,
        peak_recovery_pct = peak_recovery$recovery_pct,
        mean_recovery_pct = mean_recovery$recovery_pct
    )

    write.csv(hard_summary, file.path(out_dir, "summary_hard.csv"), row.names = FALSE)
    print(hard_summary)
}
