suppressPackageStartupMessages({
    library(dplyr)
    library(ggplot2)
})

source(file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()), "code/utils.R"))

truth_topomap_electrodes <- function(ground_truth, sensors,
                                     target_time = 0.25) {
    truth_times <- sort(unique(ground_truth$continuous_effect$time))
    selected_time <- pick_time(target_time, truth_times)
    continuous_truth_cells(
        ground_truth,
        channels = sensors$channel,
        times = selected_time
    ) %>%
        rename(value = effect) %>%
        left_join(sensors, by = "channel") %>%
        mutate(active = channel %in% ground_truth$active_electrodes)
}

truth_topomap_plot <- function(electrodes, title, fill_limit) {
    head_radius <- max(sqrt(
        electrodes$xproj^2 + electrodes$yproj^2
    )) * 1.05
    grid <- expand.grid(
        xproj = seq(-head_radius, head_radius, length.out = 140),
        yproj = seq(-head_radius, head_radius, length.out = 140)
    ) %>%
        filter(sqrt(xproj^2 + yproj^2) <= head_radius)
    grid$value <- idw_interp(
        grid,
        electrodes %>% select(xproj, yproj),
        electrodes$value,
        power = 2,
        k = 6
    )
    decor <- scalp_decorations(head_radius)

    ggplot(grid, aes(xproj, yproj)) +
        geom_raster(aes(fill = value), interpolate = TRUE) +
        geom_contour(aes(z = value), bins = 7, color = "grey25",
                     linewidth = 0.2, alpha = 0.65) +
        geom_path(data = decor$head_circle, aes(x, y), inherit.aes = FALSE,
                  linewidth = 0.5) +
        geom_path(data = decor$nose, aes(x, y), inherit.aes = FALSE,
                  linewidth = 0.5) +
        geom_path(data = decor$ear_left, aes(x, y), inherit.aes = FALSE,
                  linewidth = 0.5) +
        geom_path(data = decor$ear_right, aes(x, y), inherit.aes = FALSE,
                  linewidth = 0.5) +
        geom_point(
            data = electrodes,
            aes(xproj, yproj, shape = active),
            inherit.aes = FALSE,
            size = 1.5,
            color = "black"
        ) +
        coord_equal(
            xlim = c(-head_radius * 1.12, head_radius * 1.12),
            ylim = c(-head_radius * 1.12, head_radius * 1.18),
            expand = FALSE
        ) +
        scale_fill_gradient2(
            low = "#2166ac", mid = "white", high = "#b2182b",
            midpoint = 0, limits = c(-fill_limit, fill_limit), name = "uV"
        ) +
        scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 21)) +
        guides(shape = "none") +
        labs(title = title, subtitle = "Ground truth at 250 ms",
             x = NULL, y = NULL) +
        theme_void(base_size = 11) +
        theme(
            plot.title = element_text(face = "bold", hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            legend.position = "right"
        )
}

plot_paired_radius_truth_pilot <- function() {
    specs <- data.frame(
        radius_rule = "half-max",
        label = "Half-max (9 active)",
        data_subdir = "montecarlo_rhalfmax",
        stringsAsFactors = FALSE
    )
    sensors <- readRDS(data_path("sim_data_pilot.rds")) %>%
        distinct(channel, xproj, yproj)
    electrodes <- setNames(lapply(seq_len(nrow(specs)), function(i) {
        ground_truth <- readRDS(data_path(
            specs$data_subdir[i], "gt_001.rds"
        ))
        truth_topomap_electrodes(ground_truth, sensors, 0.25)
    }), specs$radius_rule)
    fill_limit <- max(abs(unlist(lapply(electrodes, `[[`, "value"))))
    plots <- setNames(lapply(seq_len(nrow(specs)), function(i) {
        truth_topomap_plot(
            electrodes[[specs$radius_rule[i]]],
            specs$label[i],
            fill_limit
        )
    }), specs$radius_rule)

    out_dir <- file.path(project_root(), "results/paired_radius_pilot")
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    for (rule in names(plots)) {
        ggsave(
            file.path(out_dir, paste0("truth_250ms_", rule, ".png")),
            plots[[rule]], width = 4.2, height = 4.4, dpi = 300
        )
        ggsave(
            file.path(out_dir, paste0("truth_250ms_", rule, ".pdf")),
            plots[[rule]], width = 4.2, height = 4.4
        )
    }
    invisible(plots)
}

if (!interactive() && sys.nframe() == 0) {
    plot_paired_radius_truth_pilot()
}
