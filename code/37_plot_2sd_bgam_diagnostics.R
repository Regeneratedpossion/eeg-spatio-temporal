diagnostic_channels <- function(available_channels) {
    requested <- c("CZ", "FCZ", "CPZ", "C2", "PZ", "FP1")
    missing <- setdiff(requested, available_channels)
    if (length(missing)) {
        stop(
            "Missing diagnostic channels: ",
            paste(missing, collapse = ", "),
            call. = FALSE
        )
    }
    requested
}

diagnostic_landmarks <- function(ground_truth, available_times) {
    if (is.null(ground_truth$true_onset) ||
        is.null(ground_truth$true_offset)) {
        stop("ground_truth must contain true_onset and true_offset",
             call. = FALSE)
    }
    targets <- c(
        ground_truth$true_onset,
        0.25,
        ground_truth$true_offset
    )
    unname(vapply(
        targets,
        function(target) {
            available_times[which.min(abs(available_times - target))]
        },
        numeric(1)
    ))
}

topomap_diagnostic_times <- function(available_times) {
    targets <- c(0.16, 0.25, 0.34)
    unname(vapply(
        targets,
        function(target) {
            available_times[which.min(abs(available_times - target))]
        },
        numeric(1)
    ))
}

topomap_colour_limit <- function(interpolated_data) {
    if (!"value" %in% names(interpolated_data)) {
        stop("interpolated_data must contain value", call. = FALSE)
    }
    limit <- max(abs(interpolated_data$value), na.rm = TRUE)
    if (!is.finite(limit) || limit <= 0) {
        stop("Topomap colour limit must be positive and finite",
             call. = FALSE)
    }
    limit
}

variable_contour_data <- function(grid_data, breaks = NULL) {
    if (!all(c("time_label", "value") %in% names(grid_data))) {
        stop("grid_data must contain time_label and value", call. = FALSE)
    }
    groups <- split(grid_data, grid_data$time_label, drop = TRUE)
    keep <- vapply(groups, function(group) {
        values <- group$value[is.finite(group$value)]
        if (length(values) <= 1L) {
            return(FALSE)
        }
        value_range <- range(values)
        if (is.null(breaks)) {
            return(diff(value_range) > 0)
        }
        any(breaks > value_range[1] & breaks < value_range[2])
    }, logical(1))
    if (!any(keep)) {
        return(grid_data[FALSE, , drop = FALSE])
    }
    output <- do.call(rbind, groups[keep])
    rownames(output) <- NULL
    output
}

summarise_difference_draws <- function(
        diff_draws, grid, sesoi = 0.10, odds_threshold = 10) {
    diff_draws <- as.matrix(diff_draws)
    if (ncol(diff_draws) != nrow(grid)) {
        stop("Each grid row must correspond to one draw-matrix column",
             call. = FALSE)
    }
    if (length(sesoi) != 1L || !is.finite(sesoi)) {
        stop("sesoi must be one finite number", call. = FALSE)
    }
    if (length(odds_threshold) != 1L || !is.finite(odds_threshold) ||
        odds_threshold <= 0) {
        stop("odds_threshold must be one positive finite number",
             call. = FALSE)
    }
    probability_threshold <- odds_threshold / (1 + odds_threshold)
    output <- data.frame(
        grid,
        mean = colMeans(diff_draws),
        lo = unname(apply(diff_draws, 2, quantile, probs = 0.025)),
        hi = unname(apply(diff_draws, 2, quantile, probs = 0.975)),
        p_sig = colMeans(diff_draws > sesoi),
        stringsAsFactors = FALSE
    )
    output$above_threshold <- output$p_sig > probability_threshold
    attr(output, "sesoi") <- sesoi
    attr(output, "odds_threshold") <- odds_threshold
    attr(output, "probability_threshold") <- probability_threshold
    output
}

prepare_six_channel_draws <- function(sim_model, fit, channels) {
    sensors <- sim_model %>%
        filter(channel %in% channels) %>%
        distinct(channel, xproj, yproj)
    grid <- expand_grid(
        channel = channels,
        time = sort(unique(sim_model$time))
    ) %>%
        left_join(sensors, by = "channel") %>%
        mutate(eeg_sem = 1)
    condition_1 <- grid %>%
        mutate(condition = factor("cond1", levels = c("cond1", "cond2")))
    condition_2 <- grid %>%
        mutate(condition = factor("cond2", levels = c("cond1", "cond2")))
    difference_draws <- posterior_epred(
        fit, newdata = condition_2, re_formula = NA
    ) - posterior_epred(
        fit, newdata = condition_1, re_formula = NA
    )
    list(grid = grid, difference_draws = difference_draws)
}

participant_difference_data <- function(sim_model, channels) {
    participant <- sim_model %>%
        filter(channel %in% channels) %>%
        select(participant, channel, time, condition, eeg_mean) %>%
        pivot_wider(names_from = condition, values_from = eeg_mean) %>%
        mutate(diff = cond2 - cond1)
    grand <- participant %>%
        group_by(channel, time) %>%
        summarise(diff = mean(diff), .groups = "drop")
    list(participant = participant, grand = grand)
}

probability_bands <- function(probability_data) {
    probability_data %>%
        arrange(channel, time) %>%
        group_by(channel) %>%
        mutate(run = cumsum(c(TRUE, diff(above_threshold) != 0))) %>%
        filter(above_threshold) %>%
        group_by(channel, run) %>%
        summarise(
            xmin = min(time) - 0.005,
            xmax = max(time) + 0.005,
            .groups = "drop"
        )
}

plot_bgam_fit <- function(participant, grand, posterior_summary,
                          population_truth, channels, landmarks) {
    channel_levels <- channels
    participant$channel <- factor(participant$channel, levels = channel_levels)
    grand$channel <- factor(grand$channel, levels = channel_levels)
    posterior_summary$channel <- factor(
        posterior_summary$channel, levels = channel_levels
    )
    population_truth$channel <- factor(
        population_truth$channel, levels = channel_levels
    )
    landmark_data <- data.frame(time = landmarks)

    ggplot() +
        geom_vline(
            data = landmark_data,
            aes(xintercept = time),
            linetype = "dashed", color = "grey70", linewidth = 0.3
        ) +
        geom_hline(
            yintercept = 0, linetype = "dotted",
            color = "grey55", linewidth = 0.3
        ) +
        geom_line(
            data = participant,
            aes(time, diff, group = participant),
            color = "grey65", alpha = 0.28, linewidth = 0.25
        ) +
        geom_ribbon(
            data = posterior_summary,
            aes(time, ymin = lo, ymax = hi, fill = "95% CrI"),
            alpha = 0.20
        ) +
        geom_line(
            data = grand,
            aes(time, diff, color = "Empirical grand mean"),
            linewidth = 0.75, linetype = "dotted"
        ) +
        geom_line(
            data = population_truth,
            aes(time, effect, color = "Superpopulation truth"),
            linewidth = 0.95, linetype = "longdash"
        ) +
        geom_line(
            data = posterior_summary,
            aes(time, mean, color = "BGAM posterior mean"),
            linewidth = 0.95
        ) +
        facet_wrap(~ channel, ncol = 3) +
        scale_color_manual(values = c(
            "Empirical grand mean" = "#333333",
            "Superpopulation truth" = "#2166ac",
            "BGAM posterior mean" = "#b2182b"
        )) +
        scale_fill_manual(values = c("95% CrI" = "#ef8a62")) +
        labs(
            title = "2 SD BGAM fit against superpopulation truth",
            subtitle = paste(
                "Thin grey: participant differences.",
                "Vertical lines: 200 / 250 / 300 ms."
            ),
            x = "Time (s)",
            y = "Condition difference (uV)",
            color = NULL,
            fill = NULL
        ) +
        theme_minimal(base_size = 10) +
        theme(
            legend.position = "top",
            strip.text = element_text(face = "bold"),
            plot.title = element_text(face = "bold")
        )
}

plot_posterior_probability <- function(probability_data, channels,
                                       landmarks, probability_threshold) {
    probability_data$channel <- factor(
        probability_data$channel, levels = channels
    )
    green_bands <- probability_bands(probability_data)
    green_bands$channel <- factor(green_bands$channel, levels = channels)
    landmark_data <- data.frame(time = landmarks)

    ggplot(probability_data, aes(time, p_sig)) +
        geom_rect(
            data = green_bands,
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "#1ec27d", alpha = 0.22
        ) +
        geom_vline(
            data = landmark_data,
            aes(xintercept = time),
            linetype = "dashed", color = "grey65", linewidth = 0.3
        ) +
        geom_hline(
            yintercept = 0.5, linetype = "dotted",
            color = "grey40", linewidth = 0.4
        ) +
        geom_hline(
            yintercept = probability_threshold,
            linetype = "dashed", color = "darkgreen", linewidth = 0.55
        ) +
        geom_line(color = "black", linewidth = 0.8) +
        facet_wrap(~ channel, ncol = 3) +
        coord_cartesian(ylim = c(0, 1)) +
        labs(
            title = "Posterior probability of a practically relevant effect",
            subtitle = sprintf(
                "P(cond2 - cond1 > %.1f uV | data); green = odds > %g (p > %.4f)",
                ROPE_UV, ODDS_THRESHOLD, probability_threshold
            ),
            x = "Time (s)",
            y = "Posterior probability"
        ) +
        theme_minimal(base_size = 10) +
        theme(
            strip.text = element_text(face = "bold"),
            plot.title = element_text(face = "bold")
        )
}

prepare_topomap_data <- function(sim_model, posterior_map, population_truth,
                                 ground_truth, landmarks) {
    sensors <- sim_model %>% distinct(channel, xproj, yproj)
    truth <- population_truth %>%
        filter(time %in% landmarks) %>%
        left_join(sensors, by = "channel") %>%
        transmute(channel, time, xproj, yproj, value = effect,
                  layer = "1. Ground truth (superpopulation)")
    empirical <- sim_model %>%
        filter(time %in% landmarks) %>%
        group_by(channel, time, xproj, yproj, condition) %>%
        summarise(value = mean(eeg_mean), .groups = "drop") %>%
        pivot_wider(names_from = condition, values_from = value) %>%
        transmute(channel, time, xproj, yproj, value = cond2 - cond1,
                  layer = "2. Empirical difference")
    model <- posterior_map %>%
        filter(time %in% landmarks) %>%
        transmute(channel, time, xproj, yproj, value = mean_diff,
                  layer = "4. BGAM estimate")
    empirical_minus_truth <- empirical %>%
        select(channel, time, xproj, yproj, empirical = value) %>%
        left_join(
            truth %>% select(channel, time, truth = value),
            by = c("channel", "time")
        ) %>%
        transmute(
            channel, time, xproj, yproj,
            value = empirical - truth,
            layer = "3. Noise residual (Empirical - Ground truth)"
        )
    empirical_minus_bgam <- empirical %>%
        select(channel, time, xproj, yproj, empirical = value) %>%
        left_join(
            model %>% select(channel, time, model = value),
            by = c("channel", "time")
        ) %>%
        transmute(
            channel, time, xproj, yproj,
            value = empirical - model,
            layer = "5. Empirical - BGAM"
        )
    bind_rows(
        truth, empirical, empirical_minus_truth, model, empirical_minus_bgam
    ) %>%
        mutate(
            active = channel %in% ground_truth$active_electrodes,
            time_label = factor(
                sprintf("%.0f ms", time * 1000),
                levels = sprintf("%.0f ms", landmarks * 1000)
            )
        )
}

interpolate_topomap_layers <- function(electrode_data, base_grid) {
    groups <- split(
        electrode_data,
        interaction(electrode_data$layer, electrode_data$time, drop = TRUE)
    )
    bind_rows(lapply(groups, function(electrodes) {
        base_grid %>%
            mutate(
                value = idw_interp(
                    ., electrodes %>% select(xproj, yproj), electrodes$value,
                    power = 2, k = 6
                ),
                layer = unique(electrodes$layer),
                time = unique(electrodes$time),
                time_label = unique(electrodes$time_label)
            )
    }))
}

plot_topomap_diagnostics <- function(electrode_data, landmarks) {
    head_radius <- max(
        sqrt(electrode_data$xproj^2 + electrode_data$yproj^2)
    ) * 1.05
    base_grid <- expand.grid(
        xproj = seq(-head_radius, head_radius, length.out = 140),
        yproj = seq(-head_radius, head_radius, length.out = 140)
    ) %>%
        filter(sqrt(xproj^2 + yproj^2) <= head_radius)
    interpolated <- interpolate_topomap_layers(electrode_data, base_grid)
    shared_limit <- topomap_colour_limit(interpolated)
    contour_breaks <- seq(
        -shared_limit, shared_limit, length.out = 9
    )[2:8]
    decorations <- scalp_decorations(head_radius)

    plot_layer <- function(layer_name) {
        grid_layer <- interpolated %>% filter(layer == layer_name)
        contour_layer <- variable_contour_data(
            grid_layer, breaks = contour_breaks
        )
        electrode_layer <- electrode_data %>% filter(layer == layer_name)
        ggplot(grid_layer, aes(xproj, yproj)) +
            geom_raster(aes(fill = value), interpolate = TRUE) +
            geom_contour(
                data = contour_layer,
                aes(z = value), breaks = contour_breaks, color = "grey25",
                linewidth = 0.18, alpha = 0.65
            ) +
            geom_path(
                data = decorations$head_circle, aes(x, y),
                inherit.aes = FALSE, linewidth = 0.45
            ) +
            geom_path(
                data = decorations$nose, aes(x, y),
                inherit.aes = FALSE, linewidth = 0.45
            ) +
            geom_path(
                data = decorations$ear_left, aes(x, y),
                inherit.aes = FALSE, linewidth = 0.45
            ) +
            geom_path(
                data = decorations$ear_right, aes(x, y),
                inherit.aes = FALSE, linewidth = 0.45
            ) +
            geom_point(
                data = electrode_layer,
                aes(xproj, yproj, shape = active),
                inherit.aes = FALSE, size = 0.8,
                color = "black", fill = "white"
            ) +
            facet_wrap(~ time_label, nrow = 1) +
            coord_equal(
                xlim = c(-head_radius * 1.12, head_radius * 1.12),
                ylim = c(-head_radius * 1.12, head_radius * 1.18),
                expand = FALSE
            ) +
            scale_fill_gradient2(
                low = "#2166ac", mid = "white", high = "#b2182b",
                midpoint = 0, limits = c(-shared_limit, shared_limit),
                name = "uV"
            ) +
            scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 21)) +
            guides(shape = "none") +
            labs(title = layer_name, x = NULL, y = NULL) +
            theme_void(base_size = 9) +
            theme(
                plot.title = element_text(face = "bold", size = 10),
                strip.text = element_text(face = "bold", size = 8),
                legend.position = "right"
            )
    }

    combined_plot <- plot_layer("1. Ground truth (superpopulation)") /
        plot_layer("2. Empirical difference") /
        plot_layer("3. Noise residual (Empirical - Ground truth)") /
        plot_layer("4. BGAM estimate") /
        plot_layer("5. Empirical - BGAM") +
        plot_annotation(
            title = "2 SD spatial diagnostics",
            subtitle = paste(
                "Columns: 160 / 250 / 340 ms.",
                "All panels use one symmetric colour scale."
            ),
            theme = theme(
                plot.title = element_text(face = "bold", size = 13),
                plot.subtitle = element_text(size = 9)
            )
        )
    attr(combined_plot, "topomap_colour_limit") <- shared_limit
    combined_plot
}

run_2sd_bgam_diagnostic_figures <- function() {
    root <- project_root()
    result_dir <- file.path(root, "results", "paired_radius_bgam_pilot", "2sd")
    output_dir <- file.path(result_dir, "figures")
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

    sim_data <- readRDS(data_path("montecarlo_r2sd", "sim_001.rds"))
    ground_truth <- readRDS(data_path("montecarlo_r2sd", "gt_001.rds"))
    fit <- readRDS(file.path(result_dir, "fit_seed001.rds"))
    posterior_map <- readRDS(file.path(result_dir, "posterior_map_seed001.rds"))
    sim_model <- aggregate_trials(sim_data)
    channels <- diagnostic_channels(unique(as.character(sim_model$channel)))
    times <- sort(unique(sim_model$time))
    landmarks <- diagnostic_landmarks(ground_truth, times)
    topomap_times <- topomap_diagnostic_times(times)

    population_truth <- paired_radius_superpopulation_truth(
        ground_truth,
        template_path = file.path(
            root, "resources", "brms_meeg", "code", "erp_template.R"
        )
    )
    draw_bundle <- prepare_six_channel_draws(sim_model, fit, channels)
    posterior_summary <- summarise_difference_draws(
        draw_bundle$difference_draws,
        draw_bundle$grid %>% select(channel, time),
        sesoi = ROPE_UV,
        odds_threshold = ODDS_THRESHOLD
    )
    probability_threshold <- attr(
        posterior_summary, "probability_threshold"
    )
    empirical <- participant_difference_data(sim_model, channels)
    truth_six <- population_truth %>%
        filter(channel %in% channels) %>%
        arrange(match(channel, channels), time)

    fit_plot <- plot_bgam_fit(
        empirical$participant, empirical$grand, posterior_summary,
        truth_six, channels, landmarks
    )
    probability_plot <- plot_posterior_probability(
        posterior_summary, channels, landmarks, probability_threshold
    )
    topomap_electrodes <- prepare_topomap_data(
        sim_model, posterior_map, population_truth, ground_truth, topomap_times
    )
    topomap_plot <- plot_topomap_diagnostics(topomap_electrodes, topomap_times)
    shared_limit <- attr(topomap_plot, "topomap_colour_limit")

    fit_path <- file.path(output_dir, "bgam_fit_superpopulation.png")
    probability_path <- file.path(
        output_dir, "posterior_probability_sesoi.png"
    )
    topomap_path <- file.path(
        output_dir, "topomap_diagnostics_superpopulation.png"
    )
    data_path_out <- file.path(output_dir, "diagnostic_plot_data.rds")
    ggsave(fit_path, fit_plot, width = 11, height = 7.5, dpi = 300)
    ggsave(
        probability_path, probability_plot,
        width = 11, height = 7.0, dpi = 300
    )
    ggsave(topomap_path, topomap_plot, width = 11.5, height = 13, dpi = 300)
    saveRDS(
        list(
            channels = channels,
            landmarks = landmarks,
            topomap_times = topomap_times,
            topomap_colour_limit = shared_limit,
            sesoi = ROPE_UV,
            odds_threshold = ODDS_THRESHOLD,
            probability_threshold = probability_threshold,
            participant_difference = empirical$participant,
            empirical_grand = empirical$grand,
            posterior_summary = posterior_summary,
            superpopulation_truth = population_truth,
            topomap_electrodes = topomap_electrodes
        ),
        data_path_out
    )
    invisible(c(fit_path, probability_path, topomap_path, data_path_out))
}

if (!interactive() && sys.nframe() == 0) {
    suppressPackageStartupMessages({
        library(brms)
        library(dplyr)
        library(tidyr)
        library(ggplot2)
        library(patchwork)
    })
    source(file.path(
        Sys.getenv("PROJECT_ROOT", unset = getwd()), "code", "utils.R"
    ))
    source(file.path(
        project_root(), "code/paired_radius_superpopulation_truth.R"
    ))
    run_2sd_bgam_diagnostic_figures()
}
