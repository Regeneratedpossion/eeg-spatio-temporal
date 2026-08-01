###########################################################
# Delta-based SESOI and posterior diagnostic analysis
#
# Checks whether a baseline-derived threshold on the condition
# difference reduces false positives in the full BGAM posterior map.
# Also writes the practical-null diagnostics and all-electrode posterior
# plot that previously required separate scripts and a second prediction pass.
###########################################################

library(brms)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

source(file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()), "code/utils.R"))

inputs <- slide_inputs()

out_subdir <- Sys.getenv("DIAG_OUT_SUBDIR", "delta_sesoi_probe")
out_dir <- file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()),
                     "results", "diagnostics", out_subdir)
fig_dir <- file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()),
                     "results", "figures", out_subdir)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

tag <- Sys.getenv("DIAG_TAG", Sys.getenv("SLIDE_TAG", "bgam_delta_sesoi"))
odds_threshold <- as.numeric(Sys.getenv("ODDS_THRESHOLD", "10"))
k_space <- as.numeric(Sys.getenv("K_SPACE", "2"))
baseline_ends <- as.numeric(strsplit(Sys.getenv("SESOI_BASELINE_ENDS",
                                                "0.10,0.15"), ",")[[1]])

sim_data <- readRDS(inputs$sim_file)
gt <- readRDS(inputs$gt_file)
fit <- readRDS(inputs$fit_file)

sim_model <- aggregate_trials(sim_data)
times <- sort(unique(sim_model$time))
cond_levels <- c("cond1", "cond2")

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
colnames(diff_draws) <- paste(idx$channel, idx$time, sep = "___")

diff_mean <- idx %>%
    transmute(
        cell_id = row_number(),
        channel,
        time,
        diff_mean = colMeans(diff_draws)
    )

make_sesoi <- function(baseline_end) {
    baseline_cells <- diff_mean %>%
        filter(time >= 0, time <= baseline_end)

    global_sesoi <- max(0, max(baseline_cells$diff_mean, na.rm = TRUE))

    electrode_sesoi <- baseline_cells %>%
        group_by(channel) %>%
        summarise(sesoi = max(0, max(diff_mean, na.rm = TRUE)),
                  .groups = "drop")

    bind_rows(
        tibble(
            baseline_end = baseline_end,
            threshold_type = "zero",
            channel = sort(unique(diff_mean$channel)),
            sesoi = 0
        ),
        tibble(
            baseline_end = baseline_end,
            threshold_type = "global_delta_sesoi",
            channel = sort(unique(diff_mean$channel)),
            sesoi = global_sesoi
        ),
        electrode_sesoi %>%
            mutate(baseline_end = baseline_end,
                   threshold_type = "electrode_delta_sesoi") %>%
            select(baseline_end, threshold_type, channel, sesoi)
    )
}

sesoi_table <- bind_rows(lapply(baseline_ends, make_sesoi))

compute_cells <- function(thresholds) {
    thresholds <- thresholds %>%
        arrange(channel)

    stopifnot(identical(thresholds$channel, idx$channel[match(thresholds$channel,
                                                             idx$channel)]))

    cell_threshold <- idx %>%
        select(channel, time) %>%
        left_join(thresholds %>% select(channel, sesoi), by = "channel") %>%
        pull(sesoi)

    p_sig <- colMeans(sweep(diff_draws, 2, cell_threshold, ">"))
    odds <- ifelse(p_sig >= 1, Inf, p_sig / (1 - p_sig))

    idx %>%
        transmute(channel, time, p_sig = p_sig, odds = odds,
                  flagged = odds > odds_threshold) %>%
        mutate(
            baseline_end = unique(thresholds$baseline_end),
            threshold_type = unique(thresholds$threshold_type)
        )
}

posterior_cells <- sesoi_table %>%
    group_by(baseline_end, threshold_type) %>%
    group_split() %>%
    lapply(compute_cells) %>%
    bind_rows() %>%
    left_join(sensor_mask, by = "channel") %>%
    mutate(
        in_true_window = time >= gt$true_onset & time <= gt$true_offset,
        post_offset = time > gt$true_offset,
        late_rebound_window = time >= 0.45 & time <= 0.50
    )

overall_summary <- posterior_cells %>%
    group_by(baseline_end, threshold_type) %>%
    summarise(
        odds_threshold = odds_threshold,
        n_cells = n(),
        n_flagged = sum(flagged),
        flagged_rate = mean(flagged),
        null_cells = sum(practical_null),
        null_flagged = sum(flagged & practical_null),
        practical_null_fpr = mean(flagged[practical_null]),
        active_post_offset_rate = mean(flagged[active_label & post_offset]),
        active_true_window_rate = mean(flagged[active_label & in_true_window]),
        .groups = "drop"
    ) %>%
    arrange(baseline_end, threshold_type)

by_electrode <- posterior_cells %>%
    group_by(baseline_end, threshold_type, channel, distance_to_seed,
             practical_null, active_label) %>%
    summarise(
        fpr_or_flag_rate = mean(flagged),
        max_p_sig = max(p_sig),
        max_odds = max(odds),
        first_flagged_time = if (any(flagged)) min(time[flagged]) else NA_real_,
        last_flagged_time = if (any(flagged)) max(time[flagged]) else NA_real_,
        .groups = "drop"
    ) %>%
    arrange(baseline_end, threshold_type, desc(fpr_or_flag_rate))

# The zero-threshold result is identical for every baseline choice. Keep one
# copy for the practical-null and post-offset diagnostics.
zero_cells <- posterior_cells %>%
    filter(threshold_type == "zero") %>%
    arrange(baseline_end) %>%
    distinct(channel, time, .keep_all = TRUE)

fpr_overall <- zero_cells %>%
    filter(practical_null) %>%
    summarise(
        k_space = k_space,
        odds_threshold = odds_threshold,
        n_null_cells = n(),
        n_flagged_null_cells = sum(flagged),
        practical_null_fpr = mean(flagged),
        .groups = "drop"
    )

fpr_by_electrode <- zero_cells %>%
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

post_offset_rebound <- zero_cells %>%
    left_join(raw_diff, by = c("channel", "time")) %>%
    filter(active_label, post_offset) %>%
    group_by(channel, distance_to_seed) %>%
    summarise(
        n_post_offset_time = n(),
        n_post_offset_flagged = sum(flagged),
        post_offset_flagged_rate = mean(flagged),
        max_post_offset_p_sig = max(p_sig),
        max_post_offset_odds = max(odds),
        first_post_offset_flagged_time = if_else(
            any(flagged), min(time[flagged]), NA_real_
        ),
        last_post_offset_flagged_time = if_else(
            any(flagged), max(time[flagged]), NA_real_
        ),
        mean_raw_diff_post_offset = mean(raw_diff),
        max_raw_diff_post_offset = max(raw_diff),
        mean_raw_diff_late_045_050 = mean(
            raw_diff[late_rebound_window], na.rm = TRUE
        ),
        max_raw_diff_late_045_050 = max(
            raw_diff[late_rebound_window], na.rm = TRUE
        ),
        .groups = "drop"
    ) %>%
    arrange(desc(post_offset_flagged_rate), channel)

selected_channels <- intersect(
    c("CZ", "FCZ", "CPZ", "C2", "PZ", "FP1", "FC6"),
    unique(as.character(posterior_cells$channel))
)

plot_cells <- posterior_cells %>%
    filter(channel %in% selected_channels) %>%
    mutate(
        channel = factor(channel, levels = selected_channels),
        baseline = paste0("baseline 0-", baseline_end, " s"),
        threshold_type = factor(
            threshold_type,
            levels = c("zero", "global_delta_sesoi", "electrode_delta_sesoi"),
            labels = c("Delta > 0", "global SESOI", "electrode SESOI")
        )
    )

p_threshold <- odds_threshold / (odds_threshold + 1)
landmarks <- data.frame(t = c(gt$true_onset,
                              (gt$true_onset + gt$true_offset) / 2,
                              gt$true_offset))

p <- ggplot(plot_cells, aes(x = time, y = p_sig, color = threshold_type)) +
    geom_vline(data = landmarks, aes(xintercept = t),
               inherit.aes = FALSE, linetype = "dashed",
               color = "grey60", linewidth = 0.35) +
    geom_hline(yintercept = p_threshold, linetype = "dashed",
               color = "darkgreen", linewidth = 0.35) +
    geom_line(linewidth = 0.65) +
    facet_grid(baseline ~ channel) +
    coord_cartesian(ylim = c(0, 1)) +
    scale_color_manual(values = c("Delta > 0" = "black",
                                  "global SESOI" = "#2166ac",
                                  "electrode SESOI" = "#b2182b")) +
    labs(
        title = "Delta-based SESOI probe",
        subtitle = paste0("Dashed horizontal line: posterior odds > ",
                          odds_threshold),
        x = "Time (s)",
        y = expression(Pr(Delta > threshold ~ "|" ~ data)),
        color = NULL
    ) +
    theme_minimal(base_size = 9) +
    theme(strip.text = element_text(face = "bold"),
          legend.position = "bottom")

# All-electrode posterior plot for one requested SESOI definition.
plot_baseline_end <- as.numeric(Sys.getenv("SESOI_BASELINE_END", "0.10"))
plot_threshold_type <- Sys.getenv(
    "SESOI_THRESHOLD_TYPE", "electrode_delta_sesoi"
)
stopifnot(plot_threshold_type %in%
              c("global_delta_sesoi", "electrode_delta_sesoi"))

all_cells <- posterior_cells %>%
    filter(abs(baseline_end - plot_baseline_end) < 1e-10,
           threshold_type == plot_threshold_type)
plot_thresholds <- sesoi_table %>%
    filter(abs(baseline_end - plot_baseline_end) < 1e-10,
           threshold_type == plot_threshold_type) %>%
    select(channel, sesoi)

truth_cells <- continuous_truth_cells(gt) %>%
    rename(true_effect = effect) %>%
    left_join(plot_thresholds, by = "channel") %>%
    mutate(truth_active = true_effect > sesoi)

make_bands <- function(data, active_col) {
    data %>%
        arrange(channel, time) %>%
        group_by(channel) %>%
        mutate(active = .data[[active_col]],
               run = cumsum(c(TRUE, diff(active) != 0))) %>%
        filter(active) %>%
        group_by(channel, run) %>%
        summarise(xmin = min(time), xmax = max(time), .groups = "drop")
}

detected_bands <- make_bands(all_cells, "flagged")
truth_bands <- make_bands(truth_cells, "truth_active")
channel_order <- all_cells %>%
    distinct(channel, distance_to_seed) %>%
    arrange(distance_to_seed, channel) %>%
    pull(channel)

all_cells <- all_cells %>%
    mutate(channel = factor(channel, levels = channel_order))
detected_bands <- detected_bands %>%
    mutate(channel = factor(channel, levels = channel_order))
truth_bands <- truth_bands %>%
    mutate(channel = factor(channel, levels = channel_order))

threshold_summary <- plot_thresholds %>%
    summarise(minimum = min(sesoi), median = median(sesoi), maximum = max(sesoi))
threshold_label <- if (plot_threshold_type == "global_delta_sesoi") {
    "global SESOI"
} else {
    "electrode-specific SESOI"
}
threshold_file_label <- if (plot_threshold_type == "global_delta_sesoi") {
    "global_sesoi"
} else {
    "electrode_sesoi"
}
probability_axis <- if (plot_threshold_type == "global_delta_sesoi") {
    expression(Pr(Delta > delta[global] ~ "|" ~ data))
} else {
    expression(Pr(Delta > delta[c] ~ "|" ~ data))
}

all_electrode_plot <- ggplot(all_cells, aes(time, p_sig)) +
    geom_rect(
        data = truth_bands,
        aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
        fill = "#377EB8", alpha = 0.10, inherit.aes = FALSE
    ) +
    geom_rect(
        data = detected_bands,
        aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
        fill = "#1EC27D", alpha = 0.20, inherit.aes = FALSE
    ) +
    geom_hline(yintercept = p_threshold, colour = "darkgreen",
               linetype = "dashed", linewidth = 0.28) +
    geom_line(colour = "black", linewidth = 0.38) +
    facet_wrap(~ channel, ncol = 7) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
        title = paste("All-electrode posterior probability with", threshold_label),
        subtitle = sprintf(
            paste0("High heterogeneity, seed 1; %s from 0-%.2f s baseline; ",
                   "SESOI min/median/max = %.3f / %.3f / %.3f uV"),
            threshold_label, plot_baseline_end, threshold_summary$minimum,
            threshold_summary$median, threshold_summary$maximum
        ),
        caption = paste0(
            "Blue background: true effect > electrode SESOI; green background: ",
            "posterior odds > ", odds_threshold, "; dashed horizontal: P = ",
            round(p_threshold, 3), ". Electrodes are ordered by distance from CZ."
        ),
        x = "Time (s)",
        y = probability_axis
    ) +
    theme_minimal(base_size = 8) +
    theme(
        strip.text = element_text(face = "bold", size = 7),
        axis.text = element_text(size = 6),
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 8.5),
        plot.caption = element_text(size = 7, hjust = 0)
    )

write_csv(sesoi_table,
          file.path(out_dir, paste0("delta_sesoi_thresholds_", tag, ".csv")))
write_csv(posterior_cells,
          file.path(out_dir, paste0("delta_sesoi_posterior_cells_", tag, ".csv")))
write_csv(overall_summary,
          file.path(out_dir, paste0("delta_sesoi_overall_summary_", tag, ".csv")))
write_csv(by_electrode,
          file.path(out_dir, paste0("delta_sesoi_by_electrode_", tag, ".csv")))
write_csv(zero_cells,
          file.path(out_dir, paste0("posterior_cell_flags_", tag, ".csv")))
write_csv(fpr_overall,
          file.path(out_dir, paste0("practical_null_fpr_overall_", tag, ".csv")))
write_csv(fpr_by_electrode,
          file.path(out_dir,
                    paste0("practical_null_fpr_by_electrode_", tag, ".csv")))
write_csv(post_offset_rebound,
          file.path(out_dir, paste0("post_offset_rebound_", tag, ".csv")))

ggsave(file.path(fig_dir, paste0("delta_sesoi_facets_", tag, ".pdf")),
       p, width = 12, height = 6.8, device = cairo_pdf)
ggsave(file.path(fig_dir, paste0("delta_sesoi_facets_", tag, ".png")),
       p, width = 12, height = 6.8, dpi = 220,
       device = function(...) grDevices::png(..., type = "cairo"))
ggsave(
    file.path(fig_dir, paste0("posterior_p_sig_facets_", tag,
                             "_all61_", threshold_file_label, ".pdf")),
    all_electrode_plot, width = 11.5, height = 14.5, device = cairo_pdf
)

cat("\n=== Delta SESOI overall summary ===\n")
print(overall_summary)

cat("\n=== Threshold summary ===\n")
print(sesoi_table %>%
          group_by(baseline_end, threshold_type) %>%
          summarise(min_sesoi = min(sesoi), median_sesoi = median(sesoi),
                    max_sesoi = max(sesoi), .groups = "drop"))

cat("\n=== Practical-null FPR ===\n")
print(fpr_overall)

cat("\n=== Top null electrodes by FPR ===\n")
print(head(fpr_by_electrode, 12))

cat("\n=== Post-offset rebound on active electrodes ===\n")
print(post_offset_rebound)

cat("\nWrote outputs to:\n")
cat(out_dir, "\n")
cat(fig_dir, "\n")
