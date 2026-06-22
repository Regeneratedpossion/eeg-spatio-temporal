###########################################################
# Hard-threshold diagnostics.
#
# This keeps the original active/inactive evaluation:
# active channel x time cells are compared against posterior odds.
###########################################################

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

    # Posterior summaries for each time x channel cell.
    posterior_mean <- matrix(NA_real_, length(times), length(channels))
    posterior_lower <- matrix(NA_real_, length(times), length(channels))
    prob_above_rope <- matrix(NA_real_, length(times), length(channels))

    for (t in seq_along(times)) {
        for (ch in seq_along(channels)) {
            draws <- effect_array[, t, ch]
            posterior_mean[t, ch] <- mean(draws)
            posterior_lower[t, ch] <- quantile(draws, 0.025)
            prob_above_rope[t, ch] <- mean(draws > ROPE_UV)
        }
    }

    # Convert posterior probability to odds, avoiding exact 0 or Inf.
    min_prob <- 1 / dim(effect_array)[1]
    odds_mat <- pmax(prob_above_rope, min_prob) /
        pmax(1 - prob_above_rope, min_prob)

    # Temporal recovery at the seed electrode.
    clusters <- find_clusters(times, odds_mat[, seed_idx], ODDS_THRESHOLD)
    onset_est <- clusters$cluster_onset[1]
    offset_est <- clusters$cluster_offset[1]
    temporal_error_ms <- abs(c(
        onset_est - gt$true_onset,
        offset_est - gt$true_offset
    )) * 1000
    temporal_error_ms[is.na(temporal_error_ms)] <- Inf
    temporal_pass <- max(temporal_error_ms) <= BIAS_TOLERANCE_MS

    # Spatial recovery at the true peak time.
    peak_time <- pick_time((gt$true_onset + gt$true_offset) / 2, times)
    peak_idx <- match(peak_time, times)
    pred_active_space <- posterior_lower[peak_idx, ] > 0
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
    spatial_top_pass <- channels[which.max(posterior_mean[peak_idx, ])] %in%
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
