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
    conv_summary <- rec$conv_summary
    seed_ch <- gt$seed_electrode
    active_channels <- gt$active_electrodes
    true_on <- gt$true_onset
    true_off <- gt$true_offset

    mean_mat <- apply(effect_array, c(2, 3), mean)
    lo_mat <- apply(effect_array, c(2, 3), quantile, probs = 0.025)
    hi_mat <- apply(effect_array, c(2, 3), quantile, probs = 0.975)
    p_sig_mat <- apply(effect_array > ROPE_UV, c(2, 3), mean)
    p_floor <- 1 / dim(effect_array)[1]
    odds_mat <- pmax(p_sig_mat, p_floor) / pmax(1 - p_sig_mat, p_floor)

    full <- expand.grid(time = times, channel = channels, KEEP.OUT.ATTRS = FALSE) %>%
        mutate(
            mean_diff = as.vector(mean_mat),
            lo = as.vector(lo_mat),
            hi = as.vector(hi_mat),
            p_sig = as.vector(p_sig_mat),
            odds = as.vector(odds_mat),
            pred_active = odds > ODDS_THRESHOLD,
            truth_active = channel %in% active_channels &
                time >= true_on & time <= true_off
        )

    seed_curve <- full %>% filter(channel == seed_ch) %>% arrange(time)
    clusters <- find_clusters(seed_curve$time, seed_curve$odds, ODDS_THRESHOLD)
    onset_est <- clusters$cluster_onset[1]
    offset_est <- clusters$cluster_offset[1]
    temporal_pass <- !is.na(onset_est) &&
        abs((onset_est - true_on) * 1000) <= BIAS_TOLERANCE_MS &&
        abs((offset_est - true_off) * 1000) <= BIAS_TOLERANCE_MS

    peak_time <- pick_time((true_on + true_off) / 2, times)
    peak_space <- full %>% filter(time == peak_time) %>% arrange(desc(mean_diff))
    pred_active_space <- peak_space$channel[peak_space$lo > 0]

    spatial_TP <- length(intersect(pred_active_space, active_channels))
    spatial_FP <- length(setdiff(pred_active_space, active_channels))
    spatial_FN <- length(setdiff(active_channels, pred_active_space))
    spatial_precision <- spatial_TP / (spatial_TP + spatial_FP)
    spatial_recall <- spatial_TP / (spatial_TP + spatial_FN)
    spatial_F1 <- 2 * spatial_precision * spatial_recall /
        (spatial_precision + spatial_recall)

    st_TP <- sum(full$pred_active & full$truth_active, na.rm = TRUE)
    st_FP <- sum(full$pred_active & !full$truth_active, na.rm = TRUE)
    st_FN <- sum(!full$pred_active & full$truth_active, na.rm = TRUE)
    st_precision <- st_TP / (st_TP + st_FP)
    st_recall <- st_TP / (st_TP + st_FN)
    st_F1 <- 2 * st_precision * st_recall / (st_precision + st_recall)

    seed_idx <- match(seed_ch, channels)
    peak_idx <- match(peak_time, times)
    seed_peak_draws <- effect_array[, peak_idx, seed_idx]
    window_idx <- which(times >= true_on & times <= true_off)
    window_draws <- rowMeans(effect_array[, window_idx, seed_idx, drop = FALSE][, , 1])
    truth_time <- seq(true_on, true_off, length.out = 5000)
    truth_mean <- mean(truth_temporal_bump(truth_time, gt))
    cz_neighbours <- c("CZ", "FCZ", "CPZ", "C1", "C2", "FC1", "FC2",
                       "CP1", "CP2", "Cz", "FCz", "CPz")

    convergence_pass <- max(conv_summary$n_divergent, na.rm = TRUE) == 0 &&
        max(conv_summary$n_treedepth, na.rm = TRUE) == 0 &&
        max(conv_summary$max_rhat, na.rm = TRUE) < 1.01 &&
        min(conv_summary$min_ess_bulk, na.rm = TRUE) > 400
    spatial_top_pass <- peak_space$channel[1] %in% cz_neighbours
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
