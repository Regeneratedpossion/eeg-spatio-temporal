###########################################################
# Continuous Gaussian-field diagnostics.
#
# This is better matched to the simulated signal than hard active/inactive
# labels because the spatial Gaussian effect is continuous.
###########################################################

safe_cor <- function(x, y) {
    if (sd(x) == 0 || sd(y) == 0) return(NA_real_)
    cor(x, y)
}

run_continuous_diagnostics <- function() {
    rec <- load_laplace_reconstruction()
    out_dir <- file.path(project_root(), "results", "diagnostics",
                         "laplacian_continuous")
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

    estimate_mat <- apply(rec$effect_array, c(2, 3), mean)

    sensors <- rec$sim_model %>%
        distinct(channel, xproj, yproj) %>%
        filter(channel %in% rec$channels) %>%
        arrange(match(channel, rec$channels))

    seed_coord <- sensors %>% filter(channel == rec$gt$seed_electrode)
    spatial_distance <- sqrt((sensors$xproj - seed_coord$xproj)^2 +
                             (sensors$yproj - seed_coord$yproj)^2)
    spatial_weight <- exp(-(spatial_distance^2) / (2 * rec$gt$sigma_space^2))
    temporal_weight <- truth_temporal_bump(rec$times, rec$gt)
    truth_mat <- outer(temporal_weight, spatial_weight)
    dimnames(truth_mat) <- list(as.character(rec$times), rec$channels)

    truth_vec <- as.vector(truth_mat)
    estimate_vec <- as.vector(estimate_mat)
    active_idx <- truth_vec > 0
    truth_top25_idx <- truth_vec >= 0.25 * max(truth_vec)
    peak_time <- pick_time((rec$gt$true_onset + rec$gt$true_offset) / 2,
                           rec$times)
    peak_idx <- match(peak_time, rec$times)

    cosine <- sum(truth_vec * estimate_vec) /
        (sqrt(sum(truth_vec^2)) * sqrt(sum(estimate_vec^2)))
    truth_scaled <- truth_vec / max(truth_vec)
    estimate_scaled <- pmax(estimate_vec, 0) / max(pmax(estimate_vec, 0))
    soft_dice <- 2 * sum(pmin(truth_scaled, estimate_scaled)) /
        (sum(truth_scaled) + sum(estimate_scaled))
    truth_weighted_estimate <- sum(estimate_vec[active_idx] * truth_vec[active_idx]) /
        sum(truth_vec[active_idx])
    truth_weighted_truth <- sum(truth_vec[active_idx]^2) / sum(truth_vec[active_idx])
    truth_weighted_recovery <- 100 * truth_weighted_estimate / truth_weighted_truth
    top25_recovery <- 100 * mean(estimate_vec[truth_top25_idx]) /
        mean(truth_vec[truth_top25_idx])

    summary <- data.frame(
        fit_file = rec$tag,
        dataset = data_subdir,
        spatiotemporal_cor_all = safe_cor(truth_vec, estimate_vec),
        spatiotemporal_cor_effect_window = safe_cor(truth_vec[active_idx],
                                                    estimate_vec[active_idx]),
        spatiotemporal_rmse = sqrt(mean((estimate_vec - truth_vec)^2)),
        spatiotemporal_mae = mean(abs(estimate_vec - truth_vec)),
        spatiotemporal_cosine = cosine,
        spatiotemporal_soft_dice = soft_dice,
        peak_spatial_cor = safe_cor(truth_mat[peak_idx, ],
                                    estimate_mat[peak_idx, ]),
        truth_peak_uV = max(truth_vec),
        estimate_peak_uV = max(estimate_vec),
        peak_recovery_pct = 100 * max(estimate_vec) / max(truth_vec),
        truth_weighted_mean_recovery_pct = truth_weighted_recovery,
        truth_top25_mean_recovery_pct = top25_recovery
    )

    summary_path <- file.path(out_dir, "summary_continuous.csv")
    write.csv(summary, summary_path, row.names = FALSE)
    cat("Saved:", summary_path, "\n")
    print(summary)
}
