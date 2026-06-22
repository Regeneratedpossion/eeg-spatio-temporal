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
    list(
        mean = apply(effect_array, c(2, 3), mean),
        lower = apply(effect_array, c(2, 3), quantile, probs = 0.025),
        prob_above_rope = apply(effect_array > ROPE_UV, c(2, 3), mean)
    )
}

probability_to_odds <- function(probability, n_draws) {
    min_prob <- 1 / n_draws
    pmax(probability, min_prob) / pmax(1 - probability, min_prob)
}

run_hard_diagnostics <- function() {
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

    f1_score <- function(TP, FP, FN) 2 * TP / (2 * TP + FP + FN)

    # Posterior probability: mean(draws > ROPE_UV).
    posterior <- posterior_cell_summaries(effect_array)
    odds_mat <- probability_to_odds(posterior$prob_above_rope,
                                    n_draws = dim(effect_array)[1])

    # Temporal recovery at the seed electrode.
    temporal_pass <- temporal_recovery_pass(times, odds_mat[, seed_idx], gt)

    # Spatial recovery at the true peak time.
    peak_time <- pick_time((gt$true_onset + gt$true_offset) / 2, times)
    peak_idx <- match(peak_time, times)
    pred_active_space <- posterior$lower[peak_idx, ] > 0
    spatial_F1 <- f1_score(
        TP = sum(pred_active_space & active_channels),
        FP = sum(pred_active_space & !active_channels),
        FN = sum(!pred_active_space & active_channels)
    )

    # Spatiotemporal recovery over all channel x time cells.
    pred_active <- odds_mat > ODDS_THRESHOLD
    truth_active <- outer(times >= gt$true_onset & times <= gt$true_offset,
                          active_channels, `&`)
    st_F1 <- f1_score(
        TP = sum(pred_active & truth_active),
        FP = sum(pred_active & !truth_active),
        FN = sum(!pred_active & truth_active)
    )

    # Amplitude recovery at the seed electrode.
    seed_effect <- effect_array[, , seed_idx]
    seed_peak_draws <- effect_array[, peak_idx, seed_idx]
    window_idx <- which(times >= gt$true_onset & times <= gt$true_offset)
    window_draws <- rowMeans(seed_effect[, window_idx, drop = FALSE])
    truth_time <- seq(gt$true_onset, gt$true_offset, length.out = 5000)
    truth_mean <- mean(truth_temporal_bump(truth_time, gt))
    cz_neighbours <- c("CZ", "FCZ", "CPZ", "C1", "C2", "FC1", "FC2",
                       "CP1", "CP2", "Cz", "FCz", "CPz")

    convergence_pass <- max(rec$conv_summary$n_divergent, na.rm = TRUE) == 0 &&
        max(rec$conv_summary$n_treedepth, na.rm = TRUE) == 0 &&
        max(rec$conv_summary$max_rhat, na.rm = TRUE) < 1.01 &&
        min(rec$conv_summary$min_ess_bulk, na.rm = TRUE) > 400
    spatial_top_pass <- channels[which.max(posterior$mean[peak_idx, ])] %in%
        cz_neighbours
    peak_ci <- quantile(seed_peak_draws, c(0.025, 0.975))
    mean_ci <- quantile(window_draws, c(0.025, 0.975))
    peak_recovery_pct <- 100 * mean(seed_peak_draws) / gt$bump_amplitude
    mean_recovery_pct <- 100 * mean(window_draws) / truth_mean

    summary <- data.frame(
        fit_file = rec$tag,
        dataset = data_subdir,
        convergence_pass = convergence_pass,
        temporal_pass = temporal_pass,
        spatial_top_pass = spatial_top_pass,
        peak_ci_covers_truth = peak_ci[1] <= gt$bump_amplitude &&
            gt$bump_amplitude <= peak_ci[2],
        mean_ci_covers_truth = mean_ci[1] <= truth_mean && truth_mean <= mean_ci[2],
        spatial_F1 = spatial_F1,
        spatiotemporal_F1 = st_F1,
        peak_recovery_pct = peak_recovery_pct,
        mean_recovery_pct = mean_recovery_pct
    )

    summary_path <- file.path(out_dir, "summary_hard.csv")
    write.csv(summary, summary_path, row.names = FALSE)
    cat("Saved:", summary_path, "\n")
    print(summary)
}
