###########################################################
# 08_slide_topomap_n20.R
#
# Slide figure:
#   Topomaps over the signal window, P = 20, R = 20
###########################################################

library(brms)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(tidybayes)

source(file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()), "code/utils.R"))

inputs <- slide_inputs()
dir.create(inputs$out_dir, showWarnings = FALSE, recursive = TRUE)

sim_data <- readRDS(inputs$sim_file)
gt <- readRDS(inputs$gt_file)
fit <- readRDS(inputs$fit_file)

basis <- basis_settings()
k_xy <<- basis$k_xy
k_t <<- basis$k_t
full <<- basis$full

sim_model <- aggregate_trials(sim_data)
sensors <- sim_model %>% select(channel, xproj, yproj) %>% distinct()
times <- sort(unique(sim_model$time))
target_times <- times[vapply(c(0.160, 0.250, 0.340),
                             function(x) which.min(abs(times - x)), integer(1))]

layer_levels <- c("0. Ground truth", "1. Empirical difference",
                  "2. BGAM estimate", "3. Residual")
cond_levels <- c("cond1", "cond2")
seed_xy <- sensors %>% filter(channel == gt$seed_electrode) %>% slice(1)

# Interpolation grid inside the scalp circle.
head_radius <- max(sqrt(sensors$xproj^2 + sensors$yproj^2)) * 1.05
base_grid <- expand.grid(
    xproj = seq(-head_radius, head_radius, length.out = 140),
    yproj = seq(-head_radius, head_radius, length.out = 140)
) %>%
    filter(sqrt(xproj^2 + yproj^2) <= head_radius)

temporal_weight <- function(t) {
    peak <- (gt$true_onset + gt$true_offset) / 2
    sigma_t <- (gt$true_offset - gt$true_onset) / 3
    z <- (t - peak) / sigma_t
    pmax(ifelse(abs(z) <= 1.5,
                (dnorm(z) - dnorm(1.5)) / (dnorm(0) - dnorm(1.5)),
                0), 0)
}

# Ground truth: hard-threshold spatial signal times temporal ERP weight.
gt_elec <- crossing(sensors, time = target_times) %>%
    mutate(
        d_to_seed = sqrt((xproj - seed_xy$xproj)^2 + (yproj - seed_xy$yproj)^2),
        spatial_weight = if_else(
            d_to_seed <= gt$sigma_space,
            exp(-d_to_seed^2 / (2 * gt$sigma_space^2)),
            0
        ),
        value = gt$bump_amplitude * temporal_weight(time) * spatial_weight,
        layer = "0. Ground truth"
    )

# Empirical difference: observed participant-aggregated cond2 - cond1.
emp_elec <- sim_model %>%
    filter(time %in% target_times) %>%
    group_by(channel, time, xproj, yproj, condition) %>%
    summarise(eeg = mean(eeg_mean), .groups = "drop") %>%
    pivot_wider(names_from = condition, values_from = eeg) %>%
    mutate(value = cond2 - cond1,
           layer = "1. Empirical difference")

# BGAM estimate: posterior mean of cond2 - cond1 on the same grid.
pred_grid <- expand_grid(
    channel = sort(unique(sim_model$channel)),
    time = target_times,
    condition = factor(cond_levels, levels = cond_levels)
) %>%
    left_join(sensors, by = "channel") %>%
    mutate(eeg_sem = 1)

model_elec <- add_epred_draws(pred_grid, fit, re_formula = NA) %>%
    ungroup() %>%
    select(channel, time, .draw, condition, .epred) %>%
    pivot_wider(names_from = condition, values_from = .epred) %>%
    mutate(diff = cond2 - cond1) %>%
    group_by(channel, time) %>%
    summarise(value = mean(diff), .groups = "drop") %>%
    left_join(sensors, by = "channel") %>%
    mutate(layer = "2. BGAM estimate")

# Residual: empirical difference minus BGAM estimate.
resid_elec <- emp_elec %>%
    select(channel, time, xproj, yproj, empirical = value) %>%
    left_join(model_elec %>% select(channel, time, model = value),
              by = c("channel", "time")) %>%
    mutate(value = empirical - model,
           layer = "3. Residual")

interp_layer <- function(elec_df) {
    bind_rows(lapply(split(elec_df, elec_df$time), function(e) {
        base_grid %>%
            mutate(
                value = idw_interp(., e %>% select(xproj, yproj), e$value,
                                   power = 2, k = 6),
                layer = unique(e$layer),
                time = unique(e$time)
            )
    }))
}

elec_layers <- bind_rows(gt_elec, emp_elec, model_elec, resid_elec) %>%
    mutate(layer = factor(layer, levels = layer_levels),
           active = channel %in% gt$active_electrodes)

plot_grid <- bind_rows(
    interp_layer(gt_elec),
    interp_layer(emp_elec),
    interp_layer(model_elec),
    interp_layer(resid_elec)
) %>%
    mutate(layer = factor(layer, levels = layer_levels))

shared_lim <- plot_grid %>%
    filter(layer != "3. Residual") %>%
    pull(value) %>%
    abs() %>%
    max(na.rm = TRUE)
decor <- scalp_decorations(head_radius)

plot_one <- function(layer_name) {
    ggplot(plot_grid %>% filter(layer == layer_name),
           aes(x = xproj, y = yproj)) +
        geom_raster(aes(fill = value), interpolate = TRUE) +
        geom_contour(aes(z = value), bins = 7, color = "grey25",
                     linewidth = 0.18, alpha = 0.65) +
        geom_path(data = decor$head_circle, aes(x, y), inherit.aes = FALSE,
                  linewidth = 0.45) +
        geom_path(data = decor$nose, aes(x, y), inherit.aes = FALSE,
                  linewidth = 0.45) +
        geom_path(data = decor$ear_left, aes(x, y), inherit.aes = FALSE,
                  linewidth = 0.45) +
        geom_path(data = decor$ear_right, aes(x, y), inherit.aes = FALSE,
                  linewidth = 0.45) +
        geom_point(data = elec_layers %>% filter(layer == layer_name),
                   aes(xproj, yproj, shape = active), inherit.aes = FALSE,
                   size = 0.85, color = "black", alpha = 0.8) +
        facet_wrap(~ sprintf("%.0f ms", time * 1000), nrow = 1) +
        coord_equal(xlim = c(-head_radius * 1.12, head_radius * 1.12),
                    ylim = c(-head_radius * 1.12, head_radius * 1.18),
                    expand = FALSE) +
        scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                             midpoint = 0, limits = c(-shared_lim, shared_lim),
                             name = "uV") +
        scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 21)) +
        guides(shape = "none") +
        labs(title = layer_name, x = NULL, y = NULL) +
        theme_void(base_size = 10) +
        theme(
            legend.position = "right",
            plot.background = element_rect(fill = "white", color = NA),
            panel.background = element_rect(fill = "white", color = NA),
            legend.background = element_rect(fill = "white", color = NA),
            legend.key = element_rect(fill = "white", color = NA),
            strip.text = element_text(face = "bold", size = 8),
            plot.title = element_text(face = "bold", size = 10)
        )
}

p <- plot_one("0. Ground truth") /
    plot_one("1. Empirical difference") /
    plot_one("2. BGAM estimate") /
    plot_one("3. Residual") +
    plot_annotation(
        title = "Topomaps over the signal window: P=20, R=20",
        subtitle = "Columns: onset / peak / offset. Residual = empirical - model.",
        theme = theme(
            plot.title = element_text(face = "bold", size = 12),
            plot.subtitle = element_text(color = "grey25", size = 9)
        )
    )

ggsave(file.path(inputs$out_dir, "topomap_n20_hard_sem_three_times.pdf"),
       p, width = 11.2, height = 9.2, device = cairo_pdf)
ggsave(file.path(inputs$out_dir, "topomap_n20_hard_sem_three_times.png"),
       p, width = 11.2, height = 9.2, dpi = 300,
       device = function(...) grDevices::png(..., type = "cairo"))
