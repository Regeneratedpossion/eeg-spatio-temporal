###########################################################
# 15_plot_posterior_probability.R
#
# Slide figure:
#   Posterior probability per electrode, P = 20, k = (20,20)
###########################################################

library(brms)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tidybayes)

source(file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()), "code/utils.R"))

inputs <- slide_inputs()
dir.create(inputs$out_dir, showWarnings = FALSE, recursive = TRUE)
slide_tag <- Sys.getenv("SLIDE_TAG", "n20_hard_sem")

basis <- basis_settings()
k_xy <<- basis$k_xy
k_t <<- basis$k_t
full <<- basis$full

sim_data <- readRDS(inputs$sim_file)
fit <- readRDS(inputs$fit_file)
sim_model <- aggregate_trials(sim_data)

# Plot settings
channels <- intersect(
    c("CZ", "FCZ", "CPZ", "C2", "PZ", "FP1"),
    unique(as.character(sim_model$channel))
)
cond_levels <- c("cond1", "cond2")
landmarks <- data.frame(t = c(0.160, 0.250, 0.342))
p_threshold <- ODDS_THRESHOLD / (ODDS_THRESHOLD + 1)

# Predict both conditions on selected channels, then compute P(cond2 - cond1 > 0).
pred_grid <- expand_grid(
    channel = channels,
    time = sort(unique(sim_model$time)),
    condition = factor(cond_levels, levels = cond_levels)
) %>%
    left_join(
        sim_model %>%
            filter(channel %in% channels) %>%
            select(channel, xproj, yproj) %>%
            distinct(),
        by = "channel"
    ) %>%
    mutate(eeg_sem = 1)

# For each electrode and time point, compute posterior probability
# that the condition difference is positive: P(cond2 - cond1 > 0).
post_prob <- add_epred_draws(pred_grid, fit, re_formula = NA) %>%
    ungroup() %>%
    select(channel, time, .draw, condition, .epred) %>%
    pivot_wider(names_from = condition, values_from = .epred) %>%
    mutate(diff = cond2 - cond1) %>%
    group_by(channel, time) %>%
    summarise(p_sig = mean(diff > 0), .groups = "drop") %>%
    mutate(
        channel = factor(channel, levels = channels),
        odds = p_sig / (1 - p_sig),
        above_threshold = odds > ODDS_THRESHOLD
    )

# Green bands mark continuous time ranges where posterior odds exceed threshold.
green_bands <- post_prob %>%
    arrange(channel, time) %>%
    group_by(channel) %>%
    mutate(run = cumsum(c(TRUE, diff(above_threshold) != 0))) %>%
    filter(above_threshold) %>%
    group_by(channel, run) %>%
    summarise(xmin = min(time), xmax = max(time), .groups = "drop")

p <- ggplot(post_prob, aes(x = time)) +
    geom_rect(
        data = green_bands,
        aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
        fill = "#1ec27d", alpha = 0.25, inherit.aes = FALSE
    ) +
    geom_vline(data = landmarks, aes(xintercept = t),
               linetype = "dashed", color = "grey60", linewidth = 0.4) +
    geom_hline(yintercept = 0.5, linetype = "dashed",
               color = "grey30", linewidth = 0.5) +
    geom_hline(yintercept = p_threshold, linetype = "dashed",
               color = "darkgreen", linewidth = 0.5) +
    geom_line(aes(y = p_sig), color = "black", linewidth = 0.8) +
    facet_wrap(~ channel, ncol = 3) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(
        title = sprintf("Posterior probability per electrode: P=20, R=20, k=(%d,%d)",
                        k_xy, k_t),
        subtitle = sprintf("Green band = odds > %s. Dashed verticals: true onset / peak / offset.",
                           ODDS_THRESHOLD),
        x = "Time (s)",
        y = expression(p(beta > 0 ~ "|" ~ data))
    ) +
    theme_minimal(base_size = 10) +
    theme(strip.text = element_text(face = "bold"))

ggsave(file.path(inputs$out_dir, paste0("posterior_p_sig_facets_", slide_tag, ".pdf")),
       p, width = 10, height = 6, device = cairo_pdf)
