###########################################################
# 02_diagnostics.R
#
# Full diagnostics for one fitted BGAM model.
#
# Outputs:
#   Target 1  convergence
#   Target 2  onset / offset at seed electrode
#   Target 3  spatial localisation at peak time
#   Target 3B spatiotemporal localisation over channel x time
#   Target 4  peak amplitude at seed electrode
#   Target 4B mean amplitude over the true signal window
###########################################################

library(brms)
library(dplyr)
library(tidyr)
library(tidybayes)
library(posterior)

source(file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()), "code/utils.R"))

# Input paths for the fitted model and hard-threshold dataset.
fit_file <- Sys.getenv("FIT_FILE", "bgam_fit_v5_k20_sem.rds")
data_subdir <- Sys.getenv("MC_DATA_SUBDIR", "montecarlo_n20_hard")
seed <- as.integer(Sys.getenv("SEED", "1"))
dataset <- data_subdir
sim_file <- file.path(data_subdir, sprintf("sim_%03d.rds", seed))
gt_file <- file.path(data_subdir, sprintf("gt_%03d.rds", seed))

# Read k from the fit name so diagnostics can work with different basis sizes.
resolve_k <- function(file) {
    hit <- regmatches(file, regexpr("k\\d{2}(\\d{2})?", file, perl = TRUE))
    if (!nzchar(hit)) return(list(k_xy = 20L, k_t = 20L))

    digits <- sub("^k", "", hit)
    k_xy <- as.integer(substr(digits, 1, 2))
    k_t <- if (nchar(digits) == 4) as.integer(substr(digits, 3, 4)) else k_xy
    list(k_xy = k_xy, k_t = k_t)
}

# Diagnostic settings.
k <- resolve_k(fit_file)
k_xy <<- k$k_xy
k_t <<- k$k_t
full <<- grepl("_full", fit_file, fixed = TRUE)
tag <- tools::file_path_sans_ext(fit_file)
out_dir <- diag_path("june_full")
cond_levels <- c("cond1", "cond2")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Compute posterior condition difference; return summary or raw draws.
posterior_condition_diff <- function(fit, newdata, group_vars, summarise = TRUE) {
    draws <- newdata %>%
        mutate(eeg_sem = 1) %>%
        add_epred_draws(fit, re_formula = NA) %>%
        ungroup() %>%
        select(all_of(group_vars), condition, .draw, .epred) %>%
        pivot_wider(names_from = condition, values_from = .epred) %>%
        mutate(diff = cond2 - cond1)

    if (!summarise) return(draws)

    n_draws <- length(unique(draws$.draw))

    draws %>%
        group_by(across(all_of(group_vars))) %>%
        summarise(
            mean_diff = mean(diff),
            lo = quantile(diff, 0.025),
            hi = quantile(diff, 0.975),
            p_sig = mean(diff > ROPE_UV),
            odds = pmax(p_sig, 1 / n_draws) /
                pmax(1 - p_sig, 1 / n_draws),
            .groups = "drop"
        )
}

# Compute precision, recall, and F1 from spatial or spatiotemporal counts.
f1_metrics <- function(TP, FP, FN) {
    precision <- TP / (TP + FP)
    recall <- TP / (TP + FN)
    F1 <- 2 * precision * recall / (precision + recall)

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

recovery_summary <- function(draws, truth) {
    est <- mean(draws)
    lo <- unname(quantile(draws, 0.025))
    hi <- unname(quantile(draws, 0.975))

    list(
        est = est,
        lo = lo,
        hi = hi,
        recovery_pct = 100 * est / truth,
        ci_covers_truth = lo <= truth && truth <= hi
    )
}

# Check MCMC convergence and build the Target 1 table.
convergence_diagnostics <- function(fit, fit_file, dataset) {
    np <- nuts_params(fit)
    n_div <- sum(np$Value[np$Parameter == "divergent__"])
    max_td <- fit$fit@stan_args[[1]]$control$max_treedepth
    n_tree <- sum(np$Value[np$Parameter == "treedepth__"] >= max_td)
    summ <- posterior::summarise_draws(as_draws_df(fit))

    target1 <- data.frame(
        fit_file = fit_file,
        dataset = dataset,
        metric = c("divergent_transitions", "max_treedepth_hits",
                   "max_rhat", "n_rhat_gt_1.01", "min_ess_bulk",
                   "mean_ess_bulk", "min_ess_tail"),
        value = c(n_div, n_tree,
                  max(summ$rhat, na.rm = TRUE),
                  sum(summ$rhat > 1.01, na.rm = TRUE),
                  min(summ$ess_bulk, na.rm = TRUE),
                  mean(summ$ess_bulk, na.rm = TRUE),
                  min(summ$ess_tail, na.rm = TRUE))
    )

    list(target1 = target1, n_div = n_div, summ = summ)
}

# Estimate onset and offset at the seed electrode and build Target 2.
temporal_diagnostics <- function(fit, sim_model, times, seed_ch,
                                 true_on, true_off, fit_file, dataset) {
    seed_row <- sim_model %>%
        filter(channel == seed_ch) %>%
        select(xproj, yproj) %>%
        slice(1)

    newd_cz <- expand_grid(
        time = times,
        condition = factor(cond_levels, levels = cond_levels)
    ) %>%
        mutate(xproj = seed_row$xproj, yproj = seed_row$yproj)

    diff_cz <- posterior_condition_diff(fit, newd_cz, "time") %>% arrange(time)
    clusters <- find_clusters(diff_cz$time, diff_cz$odds, ODDS_THRESHOLD)

    onset_est <- clusters$cluster_onset[1]
    offset_est <- clusters$cluster_offset[1]
    bias_onset_ms <- (onset_est - true_on) * 1000
    bias_offset_ms <- (offset_est - true_off) * 1000

    target2 <- data.frame(
        fit_file = fit_file,
        dataset = dataset,
        seed_electrode = seed_ch,
        true_onset_ms = true_on * 1000,
        true_offset_ms = true_off * 1000,
        est_onset_ms = onset_est * 1000,
        est_offset_ms = offset_est * 1000,
        bias_onset_ms = bias_onset_ms,
        bias_offset_ms = bias_offset_ms,
        pass = !is.na(onset_est) &&
            abs(bias_onset_ms) <= BIAS_TOLERANCE_MS &&
            abs(bias_offset_ms) <= BIAS_TOLERANCE_MS
    )

    list(target2 = target2, diff_cz = diff_cz, seed_row = seed_row)
}

# Evaluate spatial and spatiotemporal localisation and build Target 3/3B.
spatial_diagnostics <- function(fit, sim_model, times, true_on, true_off,
                                seed_ch, ground_truth, fit_file, dataset) {
    peak_time <- pick_time((true_on + true_off) / 2, times)
    elec <- sim_model %>% select(channel, xproj, yproj) %>% distinct()

    newd_space <- expand_grid(
        channel = elec$channel,
        condition = factor(cond_levels, levels = cond_levels)
    ) %>%
        left_join(elec, by = "channel") %>%
        mutate(time = peak_time)

    diff_space <- posterior_condition_diff(fit, newd_space, "channel") %>%
        arrange(desc(mean_diff))

    cz_neighbours <- c("CZ", "FCZ", "CPZ", "C1", "C2", "FC1", "FC2",
                       "CP1", "CP2", "Cz", "FCz", "CPz")
    top_channel <- diff_space$channel[1]
    pred_active <- diff_space$channel[diff_space$lo > 0]
    spatial <- f1_metrics(
        TP = length(intersect(pred_active, ground_truth$active_electrodes)),
        FP = length(setdiff(pred_active, ground_truth$active_electrodes)),
        FN = length(setdiff(ground_truth$active_electrodes, pred_active))
    )

    target3 <- data.frame(
        fit_file = fit_file,
        dataset = dataset,
        peak_time_ms = peak_time * 1000,
        seed_electrode = seed_ch,
        top_channel = top_channel,
        top_channel_match = top_channel %in% cz_neighbours,
        n_pred_active = length(pred_active),
        n_gt_active = length(ground_truth$active_electrodes),
        TP = spatial["TP"],
        FP = spatial["FP"],
        FN = spatial["FN"],
        precision = spatial["precision"],
        recall = spatial["recall"],
        F1 = spatial["F1"]
    )

    newd_st <- expand_grid(
        channel = elec$channel,
        time = times,
        condition = factor(cond_levels, levels = cond_levels)
    ) %>%
        left_join(elec, by = "channel")

    diff_st <- posterior_condition_diff(fit, newd_st, c("channel", "time")) %>%
        mutate(pred_active = odds > ODDS_THRESHOLD,
               truth_active = channel %in% ground_truth$active_electrodes &
                   time >= true_on & time <= true_off)

    st <- f1_metrics(
        TP = sum(diff_st$pred_active & diff_st$truth_active, na.rm = TRUE),
        FP = sum(diff_st$pred_active & !diff_st$truth_active, na.rm = TRUE),
        FN = sum(!diff_st$pred_active & diff_st$truth_active, na.rm = TRUE)
    )

    target3b <- data.frame(
        fit_file = fit_file,
        dataset = dataset,
        n_pred_active = sum(diff_st$pred_active, na.rm = TRUE),
        n_gt_active = sum(diff_st$truth_active, na.rm = TRUE),
        TP = st["TP"],
        FP = st["FP"],
        FN = st["FN"],
        precision = st["precision"],
        recall = st["recall"],
        F1 = st["F1"]
    )

    list(target3 = target3, target3b = target3b,
         diff_space = diff_space, diff_st = diff_st,
         spatial = spatial, st = st)
}

# Evaluate peak and window-mean amplitude and build Target 4/4B.
amplitude_diagnostics <- function(fit, diff_space, seed_row, times,
                                  true_on, true_off, seed_ch, truth_peak,
                                  ground_truth, fit_file, dataset) {
    cz_space <- diff_space %>% filter(channel == seed_ch) %>% slice(1)
    peak <- recovery_summary(cz_space$mean_diff, truth_peak)

    target4 <- data.frame(
        fit_file = fit_file,
        dataset = dataset,
        seed_electrode = seed_ch,
        truth_peak_uV = truth_peak,
        est_peak_uV = peak$est,
        est_lo95_uV = cz_space$lo,
        est_hi95_uV = cz_space$hi,
        recovery_pct = peak$recovery_pct,
        ci_covers_truth = cz_space$lo <= truth_peak && truth_peak <= cz_space$hi
    )

    window_times <- times[times >= true_on & times <= true_off]
    newd_window <- expand_grid(
        time = window_times,
        condition = factor(cond_levels, levels = cond_levels)
    ) %>%
        mutate(xproj = seed_row$xproj, yproj = seed_row$yproj)

    window_draws <- posterior_condition_diff(fit, newd_window, "time",
                                             summarise = FALSE) %>%
        group_by(.draw) %>%
        summarise(mean_in_window = mean(diff), .groups = "drop")

    truth_mean <- truth_window_mean(ground_truth)
    
    window <- recovery_summary(window_draws$mean_in_window, truth_mean)
    target4b <- data.frame(
        fit_file = fit_file,
        dataset = dataset,
        window_onset_ms = true_on * 1000,
        window_offset_ms = true_off * 1000,
        truth_peak_uV = truth_peak,
        truth_mean_uV = truth_mean,
        est_mean_uV = window$est,
        lo95_uV = window$lo,
        hi95_uV = window$hi,
        recovery_pct = window$recovery_pct,
        ci_covers_truth = window$ci_covers_truth,
        n_window_points = length(window_times)
    )

    list(target4 = target4, target4b = target4b)
}

write_diag <- function(x, name, out_dir, tag) {
    path <- file.path(out_dir, sprintf("%s_%s.csv", name, tag))
    write.csv(x, path, row.names = FALSE)
    cat("Saved:", path, "\n")
}

# Run all diagnostic targets and save their output tables.
run_all_diagnostics <- function(fit, sim_model, times, seed_ch, true_on,
                                true_off, truth_peak, ground_truth,
                                fit_file, dataset, out_dir, tag) {
    conv <- convergence_diagnostics(fit, fit_file, dataset)
    print(conv$target1)
    write_diag(conv$target1, "target1_convergence", out_dir, tag)

    temp <- temporal_diagnostics(fit, sim_model, times, seed_ch,
                                 true_on, true_off, fit_file, dataset)
    print(temp$target2)
    write_diag(temp$target2, "target2_temporal", out_dir, tag)
    write_diag(temp$diff_cz, "target2_temporal_full", out_dir, tag)

    spat <- spatial_diagnostics(fit, sim_model, times, true_on, true_off,
                                seed_ch, ground_truth, fit_file, dataset)
    print(spat$target3)
    write_diag(spat$target3, "target3_spatial", out_dir, tag)
    write_diag(spat$diff_space, "target3_spatial_full", out_dir, tag)

    print(spat$target3b)
    write_diag(spat$target3b, "target3b_spatiotemporal", out_dir, tag)
    write_diag(spat$diff_st, "target3b_spatiotemporal_full", out_dir, tag)

    amp <- amplitude_diagnostics(fit, spat$diff_space, temp$seed_row, times,
                                 true_on, true_off, seed_ch, truth_peak,
                                 ground_truth, fit_file, dataset)
    print(amp$target4)
    write_diag(amp$target4, "target4_peak_amplitude", out_dir, tag)

    print(amp$target4b)
    write_diag(amp$target4b, "target4b_mean_amplitude", out_dir, tag)

    summary <- data.frame(
        fit_file = fit_file,
        dataset = dataset,
        convergence_pass = conv$n_div == 0 &&
            max(conv$summ$rhat, na.rm = TRUE) < 1.01 &&
            min(conv$summ$ess_bulk, na.rm = TRUE) > 400,
        temporal_pass = temp$target2$pass,
        spatial_top_pass = spat$target3$top_channel_match,
        peak_ci_covers_truth = amp$target4$ci_covers_truth,
        mean_ci_covers_truth = amp$target4b$ci_covers_truth,
        spatial_F1 = spat$spatial["F1"],
        spatiotemporal_F1 = spat$st["F1"]
    )
    print(summary)
    write_diag(summary, "summary", out_dir, tag)

    invisible(summary)
}

# Load data.
fit <- readRDS(fit_path(fit_file))
sim_data <- readRDS(data_path(sim_file))
ground_truth <- readRDS(data_path(gt_file))
sim_model <- aggregate_trials(sim_data)
times <- sort(unique(sim_model$time))

seed_ch <- ground_truth$seed_electrode
true_on <- ground_truth$true_onset
true_off <- ground_truth$true_offset
truth_peak <- ground_truth$bump_amplitude

run_all_diagnostics(fit, sim_model, times, seed_ch, true_on,
                    true_off, truth_peak, ground_truth,
                    fit_file, dataset, out_dir, tag)
