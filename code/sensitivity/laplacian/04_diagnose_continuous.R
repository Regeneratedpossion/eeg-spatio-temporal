###########################################################
# Continuous Gaussian-field diagnostics.
#
# The simulated truth is a smooth Gaussian effect over sensors, not a binary
# active/inactive mask. These diagnostics compare the full continuous field.
###########################################################

# Correlation is undefined if one vector has no variation.
correlation_or_na <- function(x, y) {
    if (sd(x) * sd(y) == 0) return(NA_real_)
    cor(x, y)
}

continuous_truth_map <- function(ground_truth, times, channels) {
    continuous_truth_matrix(ground_truth, times, channels)
}

# Cosine similarity compares the direction/shape of two full maps.
cosine_similarity <- function(x, y) {
    sum(x * y) / (sqrt(sum(x^2)) * sqrt(sum(y^2)))
}

# Soft Dice compares overlapping continuous signal mass, not binary labels.
soft_dice_score <- function(truth, estimate) {
    truth_scaled <- truth / max(truth)
    estimate_scaled <- pmax(estimate, 0) / max(pmax(estimate, 0))

    2 * sum(pmin(truth_scaled, estimate_scaled)) /
        (sum(truth_scaled) + sum(estimate_scaled))
}

# Recovery weighted by the true signal strength.
weighted_recovery <- function(estimate, truth, index) {
    weighted_estimate <- sum(estimate[index] * truth[index]) / sum(truth[index])
    weighted_truth <- sum(truth[index]^2) / sum(truth[index])

    100 * weighted_estimate / weighted_truth
}

run_continuous_diagnostics <- function() {
    # 1. Load reconstructed posterior effect map and ground truth.
    rec <- load_laplace_reconstruction()
    out_dir <- file.path(project_root(), "results", "diagnostics",
                         "laplacian_continuous")
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

    # 2. Posterior mean estimate and continuous Gaussian truth.
    estimate_mat <- apply(rec$effect_array, c(2, 3), mean)
    truth_mat <- continuous_truth_map(rec$gt, rec$times, rec$channels)

    truth_vec <- as.vector(truth_mat)
    estimate_vec <- as.vector(estimate_mat)
    active_idx <- truth_vec > 0
    truth_top25_idx <- truth_vec >= 0.25 * max(truth_vec)

    # 3. Peak-time spatial comparison.
    peak_idx <- which.max(apply(truth_mat, 1, max))

    # 4. Continuous-map recovery metrics.
    truth_weighted_recovery <- weighted_recovery(estimate_vec, truth_vec, active_idx)
    top25_recovery <- 100 * mean(estimate_vec[truth_top25_idx]) /
        mean(truth_vec[truth_top25_idx])

    # 5. Save one-row continuous diagnostic summary.
    continuous_summary <- data.frame(
        fit_file = rec$tag,
        dataset = data_subdir,
        spatiotemporal_cor_all = correlation_or_na(truth_vec, estimate_vec),
        spatiotemporal_cor_effect_window = correlation_or_na(
            truth_vec[active_idx], estimate_vec[active_idx]
        ),
        spatiotemporal_rmse = sqrt(mean((estimate_vec - truth_vec)^2)),
        spatiotemporal_mae = mean(abs(estimate_vec - truth_vec)),
        spatiotemporal_cosine = cosine_similarity(truth_vec, estimate_vec),
        spatiotemporal_soft_dice = soft_dice_score(truth_vec, estimate_vec),
        peak_spatial_cor = correlation_or_na(truth_mat[peak_idx, ],
                                             estimate_mat[peak_idx, ]),
        truth_peak_uV = max(truth_vec),
        estimate_peak_uV = max(estimate_vec),
        peak_recovery_pct = 100 * max(estimate_vec) / max(truth_vec),
        truth_weighted_mean_recovery_pct = truth_weighted_recovery,
        truth_top25_mean_recovery_pct = top25_recovery
    )

    write.csv(continuous_summary, file.path(out_dir, "summary_continuous.csv"),
              row.names = FALSE)
    print(continuous_summary)
}
