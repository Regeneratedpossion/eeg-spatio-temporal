###########################################################
# 14_plot_bgam_fit.R
#
# Slide figure:
#   BGAM fit to the data, P = 20, k = (20,20)
###########################################################

library(brms)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidybayes)
library(patchwork)

source(file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()), "code/utils.R"))

inputs <- slide_inputs()
dir.create(inputs$out_dir, showWarnings = FALSE, recursive = TRUE)
slide_tag <- Sys.getenv("SLIDE_TAG", "bgam_fit_v5_k20_n20_hard_sem")

basis <- basis_settings()
k_xy <<- basis$k_xy
k_t <<- basis$k_t
full <<- basis$full

sim_data <- readRDS(inputs$sim_file)
fit <- readRDS(inputs$fit_file)

# Shared trial-aggregated data used for both raw curves and prediction grid.
sim_model <- aggregate_trials(sim_data)

# Plot settings
channels <- intersect(
    c("CZ", "FCZ", "CPZ", "C2", "PZ", "FP1"),
    unique(as.character(sim_data$channel))
)
cond_levels <- c("cond1", "cond2")
landmarks <- data.frame(t = c(0.160, 0.250, 0.340))
cond_cols <- c(cond1 = "#3498db", cond2 = "#e74c3c")

as_channel_factor <- function(data) {
    data %>% mutate(channel = factor(channel, levels = channels))
}

# Raw participant curves: one mean per participant x channel x condition x time.
participant_mean <- sim_model %>%
    filter(channel %in% channels) %>%
    transmute(participant, channel, condition, time, eeg = eeg_mean) %>%
    as_channel_factor()

# Build prediction grid for selected channels, then draw posterior fitted values.
posterior_predictions <- function(sim_model, fit) {
    elec <- sim_model %>%
        filter(channel %in% channels) %>%
        select(channel, xproj, yproj) %>%
        distinct()

    pred_grid <- expand_grid(
        channel = channels,
        time = sort(unique(sim_model$time)),
        condition = factor(cond_levels, levels = cond_levels)
    ) %>%
        left_join(elec, by = "channel") %>%
        mutate(eeg_sem = 1)

    add_epred_draws(pred_grid, fit, re_formula = NA)
}

summarise_prediction <- function(epred) {
    # Posterior mean and 95% credible interval for each condition curve.
    epred %>%
        group_by(channel, condition, time) %>%
        summarise(mean = mean(.epred),
                  lo = quantile(.epred, 0.025),
                  hi = quantile(.epred, 0.975),
                  .groups = "drop") %>%
        as_channel_factor()
}

# Difference is always cond2 - cond1, matching the simulated signal.
condition_difference <- function(data, value_col) {
    data %>%
        pivot_wider(names_from = condition, values_from = {{ value_col }}) %>%
        mutate(diff = cond2 - cond1) %>%
        as_channel_factor()
}

summarise_prediction_difference <- function(epred) {
    # Posterior mean and 95% credible interval for cond2 - cond1.
    epred %>%
        ungroup() %>%
        select(channel, time, .draw, condition, .epred) %>%
        condition_difference(.epred) %>%
        group_by(channel, time) %>%
        summarise(mean = mean(diff),
                  lo = quantile(diff, 0.025),
                  hi = quantile(diff, 0.975),
                  .groups = "drop")
}

plot_fit <- function(raw_pp, epred_summ, epred_diff) {
    # Dotted line: grand mean across participants.
    raw_grand <- raw_pp %>%
        group_by(channel, condition, time) %>%
        summarise(eeg = mean(eeg), .groups = "drop")

    raw_diff_pp <- condition_difference(raw_pp, eeg)
    raw_diff_grand <- condition_difference(raw_grand, eeg)

    # Panel 1: raw condition curves plus fitted condition curves.
    p_conditions <- ggplot() +
        geom_vline(data = landmarks, aes(xintercept = t),
                   linetype = "dashed", color = "grey60", linewidth = 0.3) +
        geom_hline(yintercept = 0, linetype = "dotted",
                   color = "grey55", linewidth = 0.3) +
        geom_line(data = raw_pp,
                  aes(time, eeg,
                      group = interaction(participant, condition),
                      color = condition),
                  alpha = 0.22, linewidth = 0.25) +
        geom_line(data = raw_grand,
                  aes(time, eeg, color = condition),
                  linewidth = 0.85, linetype = "dotted") +
        geom_ribbon(data = epred_summ,
                    aes(time, ymin = lo, ymax = hi, fill = condition),
                    alpha = 0.18) +
        geom_line(data = epred_summ,
                  aes(time, mean, color = condition),
                  linewidth = 0.9) +
        facet_wrap(~ channel, ncol = 3) +
        scale_color_manual(values = cond_cols) +
        scale_fill_manual(values = cond_cols) +
        labs(x = NULL, y = "EEG voltage (uV)",
             title = "Raw means and BGAM predictions",
             subtitle = "Thin: participant means; dotted: raw grand mean; solid + ribbon: posterior mean + 95% CrI") +
        theme_minimal(base_size = 9) +
        theme(legend.position = "top",
              strip.text = element_text(face = "bold"))

    # Panel 2: raw and fitted condition difference, cond2 - cond1.
    p_diff <- ggplot() +
        geom_vline(data = landmarks, aes(xintercept = t),
                   linetype = "dashed", color = "grey60", linewidth = 0.3) +
        geom_hline(yintercept = 0, linetype = "dotted",
                   color = "grey55", linewidth = 0.3) +
        geom_line(data = raw_diff_pp,
                  aes(time, diff, group = participant),
                  color = "grey50", alpha = 0.28, linewidth = 0.25) +
        geom_line(data = raw_diff_grand,
                  aes(time, diff),
                  color = "grey20", linewidth = 0.85, linetype = "dotted") +
        geom_ribbon(data = epred_diff,
                    aes(time, ymin = lo, ymax = hi),
                    fill = "black", alpha = 0.16) +
        geom_line(data = epred_diff,
                  aes(time, mean),
                  color = "black", linewidth = 0.9) +
        facet_wrap(~ channel, ncol = 3) +
        labs(x = "Time (s)", y = "Difference (cond2 - cond1, uV)",
             title = "Condition difference",
             subtitle = "Thin: participant differences; dotted: raw grand mean; solid + ribbon: posterior mean + 95% CrI") +
        theme_minimal(base_size = 9) +
        theme(strip.text = element_text(face = "bold"))

    p_conditions / p_diff +
        plot_annotation(
            title = paste0("Predictions vs raw: ",
                           tools::file_path_sans_ext(basename(inputs$fit_file))),
            subtitle = sprintf("P=20, R=20, k=(%d,%d)", k_xy, k_t),
            theme = theme(plot.title = element_text(face = "bold", size = 12))
        )
}

epred <- posterior_predictions(sim_model, fit)
p_fit <- plot_fit(
    raw_pp = participant_mean,
    epred_summ = summarise_prediction(epred),
    epred_diff = summarise_prediction_difference(epred)
)

ggsave(file.path(inputs$out_dir, paste0("pred_vs_raw_", slide_tag, ".pdf")),
       p_fit, width = 10.5, height = 9.5, device = cairo_pdf)
