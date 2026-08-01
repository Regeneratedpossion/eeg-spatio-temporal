###########################################################
# Practical-null FPR and post-offset rebound diagnostics
#
# Computes posterior evidence for cond2 - cond1 > 0 on all
# electrode x time cells, then summarizes:
#   1. practical-null false-positive rate
#   2. post-offset rebound on active/effect electrodes
###########################################################

library(brms)
library(dplyr)
library(tidyr)
library(readr)

source(file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()), "code/utils.R"))

inputs <- slide_inputs()

out_subdir <- Sys.getenv("DIAG_OUT_SUBDIR", "diagnostic_figures_39962984_4")
out_dir <- file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()),
                     "results", "diagnostics", out_subdir)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

tag <- Sys.getenv("DIAG_TAG", Sys.getenv("SLIDE_TAG", "bgam_diagnostic"))
k_space <- as.numeric(Sys.getenv("K_SPACE", "2"))
odds_threshold <- as.numeric(Sys.getenv("ODDS_THRESHOLD", "10"))

sim_data <- readRDS(inputs$sim_file)
gt <- readRDS(inputs$gt_file)
fit <- readRDS(inputs$fit_file)

sim_model <- aggregate_trials(sim_data)

sensors <- sim_model %>%
    select(channel, xproj, yproj) %>%
    distinct()

seed_xy <- sensors %>%
    filter(channel == gt$seed_electrode) %>%
    slice(1)

sensor_mask <- sensors %>%
    mutate(
        distance_to_seed = sqrt((xproj - seed_xy$xproj)^2 +
                                (yproj - seed_xy$yproj)^2),
        practical_null = distance_to_seed > k_space * gt$sigma_space,
        active_label = channel %in% gt$active_electrodes
    )

times <- sort(unique(sim_model$time))
cond_levels <- c("cond1", "cond2")

pred_grid <- expand_grid(
    channel = sort(unique(sim_model$channel)),
    time = times,
    condition = factor(cond_levels, levels = cond_levels)
) %>%
    left_join(sensors, by = "channel") %>%
    mutate(eeg_sem = 1)

message("Computing posterior_epred for ", nrow(pred_grid), " prediction rows...")
ep <- posterior_epred(fit, newdata = pred_grid, re_formula = NA)

idx <- pred_grid %>%
    mutate(row_id = row_number()) %>%
    select(row_id, channel, time, condition) %>%
    pivot_wider(names_from = condition, values_from = row_id)

diff_draws <- ep[, idx$cond2, drop = FALSE] - ep[, idx$cond1, drop = FALSE]

posterior_cells <- idx %>%
    transmute(
        channel,
        time,
        p_sig = colMeans(diff_draws > 0),
        odds = if_else(p_sig >= 1, Inf, p_sig / (1 - p_sig)),
        flagged = odds > odds_threshold
    ) %>%
    left_join(sensor_mask, by = "channel") %>%
    mutate(
        in_true_window = time >= gt$true_onset & time <= gt$true_offset,
        post_offset = time > gt$true_offset,
        late_rebound_window = time >= 0.45 & time <= 0.50
    )

fpr_overall <- posterior_cells %>%
    filter(practical_null) %>%
    summarise(
        k_space = k_space,
        odds_threshold = odds_threshold,
        n_null_cells = n(),
        n_flagged_null_cells = sum(flagged),
        practical_null_fpr = mean(flagged),
        .groups = "drop"
    )

fpr_by_electrode <- posterior_cells %>%
    filter(practical_null) %>%
    group_by(channel, distance_to_seed) %>%
    summarise(
        n_time = n(),
        n_flagged = sum(flagged),
        fpr = mean(flagged),
        max_p_sig = max(p_sig),
        max_odds = max(odds),
        first_flagged_time = if_else(any(flagged), min(time[flagged]), NA_real_),
        last_flagged_time = if_else(any(flagged), max(time[flagged]), NA_real_),
        .groups = "drop"
    ) %>%
    arrange(desc(fpr), distance_to_seed)

raw_diff <- sim_model %>%
    group_by(channel, time, condition) %>%
    summarise(eeg_mean = mean(eeg_mean), .groups = "drop") %>%
    pivot_wider(names_from = condition, values_from = eeg_mean) %>%
    mutate(raw_diff = cond2 - cond1)

post_offset_rebound <- posterior_cells %>%
    left_join(raw_diff, by = c("channel", "time")) %>%
    filter(active_label, post_offset) %>%
    group_by(channel, distance_to_seed) %>%
    summarise(
        n_post_offset_time = n(),
        n_post_offset_flagged = sum(flagged),
        post_offset_flagged_rate = mean(flagged),
        max_post_offset_p_sig = max(p_sig),
        max_post_offset_odds = max(odds),
        first_post_offset_flagged_time = if_else(any(flagged), min(time[flagged]), NA_real_),
        last_post_offset_flagged_time = if_else(any(flagged), max(time[flagged]), NA_real_),
        mean_raw_diff_post_offset = mean(raw_diff),
        max_raw_diff_post_offset = max(raw_diff),
        mean_raw_diff_late_045_050 = mean(raw_diff[late_rebound_window], na.rm = TRUE),
        max_raw_diff_late_045_050 = max(raw_diff[late_rebound_window], na.rm = TRUE),
        .groups = "drop"
    ) %>%
    arrange(desc(post_offset_flagged_rate), channel)

posterior_out <- file.path(out_dir, paste0("posterior_cell_flags_", tag, ".csv"))
fpr_overall_out <- file.path(out_dir, paste0("practical_null_fpr_overall_", tag, ".csv"))
fpr_by_electrode_out <- file.path(out_dir, paste0("practical_null_fpr_by_electrode_", tag, ".csv"))
rebound_out <- file.path(out_dir, paste0("post_offset_rebound_", tag, ".csv"))

write_csv(posterior_cells, posterior_out)
write_csv(fpr_overall, fpr_overall_out)
write_csv(fpr_by_electrode, fpr_by_electrode_out)
write_csv(post_offset_rebound, rebound_out)

cat("\n=== Practical-null FPR ===\n")
print(fpr_overall)

cat("\n=== Top null electrodes by FPR ===\n")
print(head(fpr_by_electrode, 12))

cat("\n=== Post-offset rebound on active electrodes ===\n")
print(post_offset_rebound)

cat("\nWrote:\n")
cat(posterior_out, "\n")
cat(fpr_overall_out, "\n")
cat(fpr_by_electrode_out, "\n")
cat(rebound_out, "\n")
