###########################################################
# Delta-based SESOI probe
#
# Checks whether a baseline-derived threshold on the condition
# difference reduces false positives in the full BGAM posterior map.
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

write_csv(sesoi_table,
          file.path(out_dir, paste0("delta_sesoi_thresholds_", tag, ".csv")))
write_csv(posterior_cells,
          file.path(out_dir, paste0("delta_sesoi_posterior_cells_", tag, ".csv")))
write_csv(overall_summary,
          file.path(out_dir, paste0("delta_sesoi_overall_summary_", tag, ".csv")))
write_csv(by_electrode,
          file.path(out_dir, paste0("delta_sesoi_by_electrode_", tag, ".csv")))

ggsave(file.path(fig_dir, paste0("delta_sesoi_facets_", tag, ".pdf")),
       p, width = 12, height = 6.8, device = cairo_pdf)
ggsave(file.path(fig_dir, paste0("delta_sesoi_facets_", tag, ".png")),
       p, width = 12, height = 6.8, dpi = 220,
       device = function(...) grDevices::png(..., type = "cairo"))

cat("\n=== Delta SESOI overall summary ===\n")
print(overall_summary)

cat("\n=== Threshold summary ===\n")
print(sesoi_table %>%
          group_by(baseline_end, threshold_type) %>%
          summarise(min_sesoi = min(sesoi), median_sesoi = median(sesoi),
                    max_sesoi = max(sesoi), .groups = "drop"))

cat("\nWrote outputs to:\n")
cat(out_dir, "\n")
cat(fig_dir, "\n")
