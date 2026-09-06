###########################################################
# 13_plot_simulated_eeg.R
#
# Slide figure:
#   Participant-level simulated condition differences at CZ
#   under low and high heterogeneity.
###########################################################

library(dplyr)
library(ggplot2)

project_root <- Sys.getenv("PROJECT_ROOT", unset = getwd())
scenarios <- c("low", "high")
channels <- c("CZ", "FCZ", "CPZ", "C2", "PZ", "FP1")

read_scenario <- function(scenario) {
    data_dir <- file.path(
        project_root,
        "data",
        paste0("montecarlo_p20_", scenario, "_heterogeneity")
    )
    sim_file <- file.path(data_dir, "sim_001.rds")
    gt_file <- file.path(data_dir, "gt_001.rds")

    if (!file.exists(sim_file) || !file.exists(gt_file)) {
        stop(
            "Missing seed-1 heterogeneity data for scenario '", scenario,
            "'. Expected: ", sim_file, " and ", gt_file
        )
    }

    sim_data <- readRDS(sim_file)
    ground_truth <- readRDS(gt_file)
    if (is.null(ground_truth$continuous_effect)) {
        stop("Missing ground_truth$continuous_effect in ", gt_file)
    }

    participant_by_electrode <- sim_data %>%
        filter(channel %in% channels) %>%
        mutate(channel = factor(channel, levels = channels)) %>%
        group_by(participant, channel, time) %>%
        summarise(
            difference = mean(eeg[condition == "cond2"]) -
                mean(eeg[condition == "cond1"]),
            .groups = "drop"
        ) %>%
        mutate(scenario = scenario)

    group_by_electrode <- participant_by_electrode %>%
        group_by(scenario, channel, time) %>%
        summarise(difference = mean(difference), .groups = "drop")

    truth_by_electrode <- ground_truth$continuous_effect %>%
        filter(channel %in% channels) %>%
        transmute(
            scenario = scenario,
            channel = factor(channel, levels = channels),
            time,
            difference = effect
        )

    list(
        participant = participant_by_electrode %>% filter(channel == "CZ"),
        group = group_by_electrode %>% filter(channel == "CZ"),
        truth = truth_by_electrode %>% filter(channel == "CZ"),
        electrode_participant = participant_by_electrode,
        electrode_group = group_by_electrode,
        electrode_truth = truth_by_electrode
    )
}

scenario_data <- lapply(scenarios, read_scenario)
participant_diff <- bind_rows(lapply(scenario_data, `[[`, "participant"))
group_mean <- bind_rows(lapply(scenario_data, `[[`, "group"))
truth <- bind_rows(lapply(scenario_data, `[[`, "truth"))
electrode_participant <- bind_rows(
    lapply(scenario_data, `[[`, "electrode_participant")
)
electrode_group <- bind_rows(lapply(scenario_data, `[[`, "electrode_group"))
electrode_truth <- bind_rows(lapply(scenario_data, `[[`, "electrode_truth"))

scenario_labels <- c(
    low = "Low",
    high = "High"
)
participant_diff$scenario <- factor(
    participant_diff$scenario,
    levels = scenarios,
    labels = scenario_labels
)
group_mean$scenario <- factor(
    group_mean$scenario,
    levels = scenarios,
    labels = scenario_labels
)
truth$scenario <- factor(
    truth$scenario,
    levels = scenarios,
    labels = scenario_labels
)

p <- ggplot() +
    geom_hline(yintercept = 0, color = "grey70", linewidth = 0.3) +
    geom_line(
        data = participant_diff,
        aes(time, difference, group = interaction(scenario, participant)),
        color = "grey68", alpha = 0.55, linewidth = 0.28
    ) +
    geom_line(
        data = group_mean,
        aes(time, difference, color = "Observed group mean"),
        linewidth = 0.9
    ) +
    geom_line(
        data = truth,
        aes(time, difference, color = "Saved continuous truth"),
        linewidth = 1.0, linetype = "22"
    ) +
    facet_wrap(~scenario, nrow = 1) +
    scale_color_manual(
        values = c(
            "Observed group mean" = "#0072B2",
            "Saved continuous truth" = "#D55E00"
        )
    ) +
    labs(
        x = "Time (s)",
        y = expression("Condition difference (" * mu * "V)"),
        color = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
        legend.position = "top",
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", size = 11),
        plot.margin = margin(2, 4, 2, 2)
    )

out_file <- Sys.getenv(
    "SLIDE_HETEROGENEITY_FIG",
    unset = file.path(project_root, "slides", "simulated_heterogeneity_cz_seed1.pdf")
)
dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
ggsave(out_file, p, width = 10, height = 5, device = cairo_pdf)
message("Saved: ", out_file)

# Appendix figures: use the same electrodes and y-axis scale in all scenarios.
y_limits <- range(
    c(electrode_participant$difference, electrode_truth$difference),
    finite = TRUE
)
y_padding <- 0.04 * diff(y_limits)
y_limits <- y_limits + c(-y_padding, y_padding)
out_dir <- dirname(out_file)

for (scenario in scenarios) {
    p_electrodes <- ggplot() +
        geom_hline(yintercept = 0, color = "grey70", linewidth = 0.25) +
        geom_line(
            data = electrode_participant %>%
                filter(.data$scenario == .env$scenario),
            aes(time, difference, group = participant),
            color = "grey68", alpha = 0.5, linewidth = 0.22
        ) +
        geom_line(
            data = electrode_group %>%
                filter(.data$scenario == .env$scenario),
            aes(time, difference, color = "Observed group mean"),
            linewidth = 0.75
        ) +
        geom_line(
            data = electrode_truth %>%
                filter(.data$scenario == .env$scenario),
            aes(time, difference, color = "Saved continuous truth"),
            linewidth = 0.85, linetype = "22"
        ) +
        facet_wrap(~channel, ncol = 3) +
        coord_cartesian(ylim = y_limits) +
        scale_color_manual(
            values = c(
                "Observed group mean" = "#0072B2",
                "Saved continuous truth" = "#D55E00"
            )
        ) +
        labs(
            x = "Time (s)",
            y = expression("Condition difference (" * mu * "V)"),
            color = NULL
        ) +
        theme_minimal(base_size = 10) +
        theme(
            legend.position = "top",
            panel.grid.minor = element_blank(),
            strip.text = element_text(face = "bold"),
            plot.margin = margin(2, 4, 2, 2)
        )

    scenario_file <- file.path(
        out_dir,
        paste0("simulated_electrodes_", scenario, "_seed1.pdf")
    )
    ggsave(scenario_file, p_electrodes, width = 10, height = 5.6,
           device = cairo_pdf)
    message("Saved: ", scenario_file)
}
