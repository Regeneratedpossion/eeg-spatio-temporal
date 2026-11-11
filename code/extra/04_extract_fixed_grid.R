###########################################################
# 04_extract_fixed_grid.R
#
# Per-seed extraction for the aggregated 3-plot diagnostic
# (fit curve / posterior probability / topomap) across the
# scenario<COND> batch.
#
# For ONE seed, computes:
#   - six_channel.csv : posterior mean_diff + p_sig for the
#     6 diagnostic channels (CZ, FCZ, CPZ, C2, PZ, FP1) at
#     every time point, plus the matching empirical
#     (participant-level) cond2-cond1 difference.
#   - topomap_3t.csv   : posterior mean_diff + empirical
#     cond2-cond1 for ALL channels at the 3 fixed topomap
#     times (~160/250/340 ms).
#   - sesoi.csv        : the baseline-calibrated global SESOI
#     for this seed (for QC only; the aggregate plot uses the
#     across-seed average).
#
# Usage:
#   Rscript code/extra/04_extract_fixed_grid.R <seed_int>
#
# Heavy step is posterior_epred() on the fit; this script does
# NOT keep the full posterior_map / draws in the output, only
# the fixed-grid summaries needed downstream, so per-seed
# output is a few KB instead of hundreds of MB.
###########################################################

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(brms)
})

args <- commandArgs(trailingOnly = TRUE)
seed_task <- if (length(args) >= 1) {
    as.integer(args[[1]])
} else {
    as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA))
}
if (is.na(seed_task)) {
    stop("Seed not provided: pass as arg[1] or set SLURM_ARRAY_TASK_ID",
         call. = FALSE)
}

source(file.path(Sys.getenv("PROJECT_ROOT", unset = here::here()), "code/R/utils.R"))
source(file.path(project_root(), "code/R/settings.R"))
source(file.path(project_root(), "code/R/plot_diagnostics.R"))

seed_str <- sprintf("%03d", seed_task)
result_dir <- file.path(
    project_root(), "results",
    fixed_grid_batch(),
    paste0("seed_", seed_str)
)
out_dir <- file.path(result_dir, "fixed_grid")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("=== seed", seed_task, "===\n")

data_subdir <- fixed_grid_batch()   # the data folder has the same name as the results folder

sim_data <- readRDS(data_path(data_subdir, paste0("sim_", seed_str, ".rds")))
ground_truth <- readRDS(data_path(data_subdir, paste0("gt_", seed_str, ".rds")))
fit <- readRDS(file.path(result_dir, "fit.rds"))
sim_model <- aggregate_trials(sim_data)

# -----------------------------------------------------------
# Baseline calibration (needed to reproduce the same global
# SESOI logic as the original single-fit script). Restrict the
# posterior_epred call to time <= 0.10 across all channels.
# -----------------------------------------------------------
baseline_end <- 0.10
sensors <- sim_model %>% distinct(channel, xproj, yproj) %>% arrange(channel)
all_times <- sort(unique(sim_model$time))
baseline_times <- all_times[all_times >= 0 & all_times <= baseline_end]

baseline_grid <- expand_grid(channel = sensors$channel, time = baseline_times) %>%
    left_join(sensors, by = "channel") %>%
    mutate(eeg_sem = 1)
b1 <- baseline_grid %>% mutate(condition = factor("cond1", levels = c("cond1", "cond2")))
b2 <- baseline_grid %>% mutate(condition = factor("cond2", levels = c("cond1", "cond2")))
baseline_draws <- posterior_epred(fit, newdata = b2, re_formula = NA) -
    posterior_epred(fit, newdata = b1, re_formula = NA)
baseline_map <- baseline_grid %>% mutate(mean_diff = colMeans(baseline_draws))
calibration_cell <- baseline_map %>% slice_max(mean_diff, n = 1, with_ties = FALSE)
global_sesoi <- max(0, calibration_cell$mean_diff)
rm(baseline_draws, baseline_grid, b1, b2)
gc(FALSE)

write.csv(
    data.frame(
        seed = seed_task,
        channel = calibration_cell$channel,
        time = calibration_cell$time,
        global_sesoi = global_sesoi
    ),
    file.path(out_dir, "sesoi.csv"),
    row.names = FALSE
)

# -----------------------------------------------------------
# Six diagnostic channels x all times (fit curve + posterior
# probability plot inputs).
# -----------------------------------------------------------
channels6 <- diagnostic_channels(unique(as.character(sim_model$channel)))
draw_bundle <- prepare_six_channel_draws(sim_model, fit, channels6)
posterior_summary6 <- summarise_difference_draws(
    draw_bundle$difference_draws,
    draw_bundle$grid %>% select(channel, time),
    sesoi = global_sesoi,
    odds_threshold = MAIN_ODDS
)
empirical6 <- participant_difference_data(sim_model, channels6)

six_out <- posterior_summary6 %>%
    select(channel, time, model_mean_diff = mean, model_p_sig = p_sig) %>%
    left_join(
        empirical6$grand %>% rename(empirical_diff = diff),
        by = c("channel", "time")
    ) %>%
    mutate(seed = seed_task, .before = 1)
write.csv(six_out, file.path(out_dir, "six_channel.csv"), row.names = FALSE)
rm(draw_bundle, posterior_summary6, empirical6)
gc(FALSE)

# -----------------------------------------------------------
# All channels x 3 fixed topomap times (~160/250/340 ms).
# -----------------------------------------------------------
topomap_times <- topomap_diagnostic_times(all_times)
topo_grid <- expand_grid(channel = sensors$channel, time = topomap_times) %>%
    left_join(sensors, by = "channel") %>%
    mutate(eeg_sem = 1)
t1 <- topo_grid %>% mutate(condition = factor("cond1", levels = c("cond1", "cond2")))
t2 <- topo_grid %>% mutate(condition = factor("cond2", levels = c("cond1", "cond2")))
topo_draws <- posterior_epred(fit, newdata = t2, re_formula = NA) -
    posterior_epred(fit, newdata = t1, re_formula = NA)
model_topo <- topo_grid %>%
    transmute(
        channel, time, xproj, yproj,
        model_mean_diff = colMeans(topo_draws),
        model_p_sig = colMeans(topo_draws > global_sesoi)
    )
rm(topo_draws, t1, t2)
gc(FALSE)

empirical_topo <- sim_model %>%
    filter(time %in% topomap_times) %>%
    group_by(channel, time, xproj, yproj, condition) %>%
    summarise(value = mean(eeg_mean), .groups = "drop") %>%
    pivot_wider(names_from = condition, values_from = value) %>%
    transmute(channel, time, xproj, yproj, empirical_diff = cond2 - cond1)

topo_out <- model_topo %>%
    left_join(
        empirical_topo %>% select(channel, time, empirical_diff),
        by = c("channel", "time")
    ) %>%
    mutate(seed = seed_task, .before = 1)
write.csv(topo_out, file.path(out_dir, "topomap_3t.csv"), row.names = FALSE)

cat("done seed", seed_task, "-> ", out_dir, "\n")
