###########################################################
# 13_laplacian_recovery_diagnostics.R
#
# Full recovery diagnostics for the Laplacian low-rank basis model.
# Reconstructs channel x time posterior effects from basis-score
# posteriors, then computes convergence, temporal, spatial,
# spatiotemporal, peak, and window-mean recovery diagnostics.
###########################################################

suppressPackageStartupMessages({
    library(brms)
    library(dplyr)
    library(tidyr)
    library(posterior)
})

source(file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()), "code/utils.R"))

fit_subdir <- Sys.getenv("LAPLACE_FIT_SUBDIR", "laplacian_basis_brms_fixed")
fit_dir <- file.path(project_root(), "results", fit_subdir)
out_dir <- file.path(project_root(), "results", "diagnostics",
                     "laplacian_basis_fixed")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

data_subdir <- Sys.getenv("LAPLACE_DATA_SUBDIR",
                          "montecarlo_p5_gaussian_hardlabel")
seed <- as.integer(Sys.getenv("LAPLACE_SEED", "1"))
adapt_tag <- Sys.getenv("LAPLACE_ADAPT_TAG", "0p999")
k_basis <- as.integer(Sys.getenv("LAPLACE_K", "10"))
time_k <- as.integer(Sys.getenv("LAPLACE_TIME_K", "10"))
iter <- as.integer(Sys.getenv("LAPLACE_ITER", "3000"))
warmup <- as.integer(Sys.getenv("LAPLACE_WARMUP", "1000"))
tag <- sprintf("laplacian_k%d_seed%03d_iter%d_warmup%d_adapt%s",
               k_basis, seed, iter, warmup, adapt_tag)
tag_suffix <- Sys.getenv("LAPLACE_TAG_SUFFIX", "")
if (nzchar(tag_suffix)) {
    tag <- paste0(tag, "_", tag_suffix)
}

cond_levels <- c("cond1", "cond2")
gt <- readRDS(data_path(data_subdir, sprintf("gt_%03d.rds", seed)))
sim <- readRDS(data_path(data_subdir, sprintf("sim_%03d.rds", seed)))
sim_model <- aggregate_trials(sim)
times <- sort(unique(sim_model$time))
channels <- sort(unique(sim_model$channel))

summary_file <- file.path(
    fit_dir,
    sprintf("summary_all_iter%d_warmup%d_adapt%s.csv", iter, warmup, adapt_tag)
)
conv_summary <- read.csv(summary_file)

fit_files <- conv_summary$fit_file[order(conv_summary$basis_task)]
if (!all(file.exists(fit_files))) {
    stop("Some Laplacian fit files are missing.")
}

write_diag <- function(x, name) {
    path <- file.path(out_dir, sprintf("%s_%s.csv", name, tag))
    write.csv(x, path, row.names = FALSE)
    cat("Saved:", path, "\n")
}

f1_metrics <- function(TP, FP, FN) {
    precision <- if ((TP + FP) > 0) TP / (TP + FP) else NA_real_
    recall <- if ((TP + FN) > 0) TP / (TP + FN) else NA_real_
    F1 <- if (is.na(precision) || is.na(recall) || precision + recall == 0) {
        NA_real_
    } else {
        2 * precision * recall / (precision + recall)
    }
    c(TP = TP, FP = FP, FN = FN,
      precision = precision, recall = recall, F1 = F1)
}

truth_window_mean <- function(ground_truth, n = 5000) {
    true_on <- ground_truth$true_onset
    true_off <- ground_truth$true_offset
    peak_t <- 0.250
    sigma_t <- 0.050

    z <- (seq(true_on, true_off, length.out = n) - peak_t) / sigma_t
    val <- exp(-z^2 / 2)
    val_trunc <- exp(-1.5^2 / 2)
    val[abs(z) > 1.5] <- val_trunc
    val_norm <- pmax((val - val_trunc) / (1 - val_trunc), 0)

    ground_truth$bump_amplitude * mean(val_norm)
}

odds_from_p <- function(p, n_draws) {
    floor <- 1 / n_draws
    pmax(p, floor) / pmax(1 - p, floor)
}

condition_score_draws <- function(obj, times) {
    fit <- obj$fit
    base <- data.frame(
        time = times,
        score = 0,
        score_sem = 1,
        participant = obj$score_data$participant[1]
    )
    nd1 <- base
    nd2 <- base
    nd1$condition <- factor("cond1", levels = cond_levels)
    nd2$condition <- factor("cond2", levels = cond_levels)

    ep1 <- posterior_epred(fit, newdata = nd1, re_formula = NA)
    ep2 <- posterior_epred(fit, newdata = nd2, re_formula = NA)
    ep2 - ep1
}

first_obj <- readRDS(fit_files[1])
channels <- first_obj$channels
n_draws <- nrow(posterior_epred(
    first_obj$fit,
    newdata = data.frame(
        time = times[1],
        score = 0,
        score_sem = 1,
        participant = first_obj$score_data$participant[1],
        condition = factor("cond1", levels = cond_levels)
    ),
    re_formula = NA
))

effect_array <- array(
    0,
    dim = c(n_draws, length(times), length(channels)),
    dimnames = list(NULL, as.character(times), channels)
)

for (file in fit_files) {
    obj <- readRDS(file)
    score_diff <- condition_score_draws(obj, times)
    basis_vector <- obj$basis_vector[channels]
    for (j in seq_along(channels)) {
        effect_array[, , j] <- effect_array[, , j] +
            score_diff * basis_vector[j]
    }
}

mean_mat <- apply(effect_array, c(2, 3), mean)
lo_mat <- apply(effect_array, c(2, 3), quantile, probs = 0.025)
hi_mat <- apply(effect_array, c(2, 3), quantile, probs = 0.975)
p_sig_mat <- apply(effect_array > ROPE_UV, c(2, 3), mean)
odds_mat <- odds_from_p(p_sig_mat, n_draws)

full <- expand.grid(
    time = times,
    channel = channels,
    KEEP.OUT.ATTRS = FALSE
) %>%
    mutate(
        mean_diff = as.vector(mean_mat),
        lo = as.vector(lo_mat),
        hi = as.vector(hi_mat),
        p_sig = as.vector(p_sig_mat),
        odds = as.vector(odds_mat),
        pred_active = odds > ODDS_THRESHOLD,
        truth_active = channel %in% gt$active_electrodes &
            time >= gt$true_onset & time <= gt$true_offset
    )

target1 <- data.frame(
        fit_file = tag,
        dataset = data_subdir,
        metric = c("basis_models", "max_divergent", "max_treedepth_hits",
                   "max_rhat", "max_n_rhat_gt_1.01", "min_ess_bulk",
                   "min_ess_tail", "total_elapsed_sec"),
        value = c(
            nrow(conv_summary),
            max(conv_summary$n_divergent, na.rm = TRUE),
            max(conv_summary$n_treedepth, na.rm = TRUE),
            max(conv_summary$max_rhat, na.rm = TRUE),
            max(conv_summary$n_rhat_gt_1_01, na.rm = TRUE),
            min(conv_summary$min_ess_bulk, na.rm = TRUE),
            min(conv_summary$min_ess_tail, na.rm = TRUE),
            sum(conv_summary$elapsed_sec, na.rm = TRUE)
        )
    )
write_diag(target1, "target1_convergence")

seed_ch <- gt$seed_electrode
seed_curve <- full %>%
    filter(channel == seed_ch) %>%
    arrange(time)
clusters <- find_clusters(seed_curve$time, seed_curve$odds, ODDS_THRESHOLD)
onset_est <- clusters$cluster_onset[1]
offset_est <- clusters$cluster_offset[1]
target2 <- data.frame(
    fit_file = tag,
    dataset = data_subdir,
    seed_electrode = seed_ch,
    true_onset_ms = gt$true_onset * 1000,
    true_offset_ms = gt$true_offset * 1000,
    est_onset_ms = onset_est * 1000,
    est_offset_ms = offset_est * 1000,
    bias_onset_ms = (onset_est - gt$true_onset) * 1000,
    bias_offset_ms = (offset_est - gt$true_offset) * 1000,
    pass = !is.na(onset_est) &&
        abs((onset_est - gt$true_onset) * 1000) <= BIAS_TOLERANCE_MS &&
        abs((offset_est - gt$true_offset) * 1000) <= BIAS_TOLERANCE_MS
)
write_diag(target2, "target2_temporal")
write_diag(seed_curve, "target2_temporal_full")

peak_time <- pick_time((gt$true_onset + gt$true_offset) / 2, times)
peak_space <- full %>%
    filter(time == peak_time) %>%
    arrange(desc(mean_diff))
pred_active_space <- peak_space$channel[peak_space$lo > 0]
spatial <- f1_metrics(
    TP = length(intersect(pred_active_space, gt$active_electrodes)),
    FP = length(setdiff(pred_active_space, gt$active_electrodes)),
    FN = length(setdiff(gt$active_electrodes, pred_active_space))
)
cz_neighbours <- c("CZ", "FCZ", "CPZ", "C1", "C2", "FC1", "FC2",
                   "CP1", "CP2", "Cz", "FCz", "CPz")
target3 <- data.frame(
    fit_file = tag,
    dataset = data_subdir,
    peak_time_ms = peak_time * 1000,
    seed_electrode = seed_ch,
    top_channel = peak_space$channel[1],
    top_channel_match = peak_space$channel[1] %in% cz_neighbours,
    n_pred_active = length(pred_active_space),
    n_gt_active = length(gt$active_electrodes),
    TP = spatial["TP"],
    FP = spatial["FP"],
    FN = spatial["FN"],
    precision = spatial["precision"],
    recall = spatial["recall"],
    F1 = spatial["F1"]
)
write_diag(target3, "target3_spatial")
write_diag(peak_space, "target3_spatial_full")

st <- f1_metrics(
    TP = sum(full$pred_active & full$truth_active, na.rm = TRUE),
    FP = sum(full$pred_active & !full$truth_active, na.rm = TRUE),
    FN = sum(!full$pred_active & full$truth_active, na.rm = TRUE)
)
target3b <- data.frame(
    fit_file = tag,
    dataset = data_subdir,
    n_pred_active = sum(full$pred_active, na.rm = TRUE),
    n_gt_active = sum(full$truth_active, na.rm = TRUE),
    TP = st["TP"],
    FP = st["FP"],
    FN = st["FN"],
    precision = st["precision"],
    recall = st["recall"],
    F1 = st["F1"]
)
write_diag(target3b, "target3b_spatiotemporal")
write_diag(full, "target3b_spatiotemporal_full")

seed_idx <- match(seed_ch, channels)
peak_idx <- match(peak_time, times)
seed_peak_draws <- effect_array[, peak_idx, seed_idx]
target4 <- data.frame(
    fit_file = tag,
    dataset = data_subdir,
    seed_electrode = seed_ch,
    truth_peak_uV = gt$bump_amplitude,
    est_peak_uV = mean(seed_peak_draws),
    est_lo95_uV = unname(quantile(seed_peak_draws, 0.025)),
    est_hi95_uV = unname(quantile(seed_peak_draws, 0.975)),
    recovery_pct = 100 * mean(seed_peak_draws) / gt$bump_amplitude,
    ci_covers_truth = unname(quantile(seed_peak_draws, 0.025)) <= gt$bump_amplitude &&
        gt$bump_amplitude <= unname(quantile(seed_peak_draws, 0.975))
)
write_diag(target4, "target4_peak_amplitude")

window_idx <- which(times >= gt$true_onset & times <= gt$true_offset)
window_draws <- rowMeans(effect_array[, window_idx, seed_idx, drop = FALSE][, , 1])
truth_mean <- truth_window_mean(gt)
target4b <- data.frame(
    fit_file = tag,
    dataset = data_subdir,
    window_onset_ms = gt$true_onset * 1000,
    window_offset_ms = gt$true_offset * 1000,
    truth_peak_uV = gt$bump_amplitude,
    truth_mean_uV = truth_mean,
    est_mean_uV = mean(window_draws),
    lo95_uV = unname(quantile(window_draws, 0.025)),
    hi95_uV = unname(quantile(window_draws, 0.975)),
    recovery_pct = 100 * mean(window_draws) / truth_mean,
    ci_covers_truth = unname(quantile(window_draws, 0.025)) <= truth_mean &&
        truth_mean <= unname(quantile(window_draws, 0.975)),
    n_window_points = length(window_idx)
)
write_diag(target4b, "target4b_mean_amplitude")

summary <- data.frame(
    fit_file = tag,
    dataset = data_subdir,
    convergence_pass = max(conv_summary$n_divergent, na.rm = TRUE) == 0 &&
        max(conv_summary$n_treedepth, na.rm = TRUE) == 0 &&
        max(conv_summary$max_rhat, na.rm = TRUE) < 1.01 &&
        min(conv_summary$min_ess_bulk, na.rm = TRUE) > 400,
    temporal_pass = target2$pass,
    spatial_top_pass = target3$top_channel_match,
    peak_ci_covers_truth = target4$ci_covers_truth,
    mean_ci_covers_truth = target4b$ci_covers_truth,
    spatial_F1 = target3$F1,
    spatiotemporal_F1 = target3b$F1,
    peak_recovery_pct = target4$recovery_pct,
    mean_recovery_pct = target4b$recovery_pct
)
write_diag(summary, "summary")
print(summary)
