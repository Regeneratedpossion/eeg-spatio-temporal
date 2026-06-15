###########################################################
# 05_slide_simulated_eeg_n20.R
#
# Slide figure:
#   Simulated EEG data, P = 20, R = 20
###########################################################

library(dplyr)
library(ggplot2)

source(file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()), "code/utils.R"))

inputs <- slide_inputs()
dir.create(inputs$out_dir, showWarnings = FALSE, recursive = TRUE)
sim_data <- readRDS(inputs$sim_file)

# Plot settings
channels <- intersect(
    c("CZ", "FCZ", "CPZ", "C2", "PZ", "FP1"),
    unique(as.character(sim_data$channel))
)
cond_levels <- c("cond1", "cond2")
landmarks <- data.frame(t = c(0.160, 0.250, 0.342))
cond_cols <- c(cond1 = "#3498db", cond2 = "#e74c3c")

# Average trials within each participant before plotting group summaries.
participant_mean <- sim_data %>%
    filter(channel %in% channels) %>%
    mutate(condition = factor(condition, levels = cond_levels),
           channel = factor(channel, levels = channels)) %>%
    group_by(participant, channel, condition, time) %>%
    summarise(eeg = mean(eeg), .groups = "drop")

# Plot group mean EEG curves; ribbon = mean +/- 1 SE across participants.
overview <- participant_mean %>%
    group_by(channel, condition, time) %>%
    summarise(mean = mean(eeg),
              se = sd(eeg) / sqrt(dplyr::n()),
              .groups = "drop")

p_overview <- ggplot(overview, aes(time, mean,
                                   color = condition, fill = condition)) +
    geom_vline(data = landmarks, aes(xintercept = t),
               inherit.aes = FALSE, linetype = "dashed",
               color = "grey60", linewidth = 0.3) +
    geom_hline(yintercept = 0, linetype = "dotted",
               color = "grey55", linewidth = 0.3) +
    # Shaded uncertainty band around the group mean.
    geom_ribbon(aes(ymin = mean - se, ymax = mean + se),
                alpha = 0.18, color = NA) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~ channel, ncol = 3) +
    scale_color_manual(values = cond_cols) +
    scale_fill_manual(values = cond_cols) +
    labs(title = "Simulated EEG data: P=20, R=20",
         subtitle = "Group mean +/- 1 SE across participants. Dashed verticals: true onset / peak / offset.",
         x = "Time (s)", y = "EEG voltage (uV)",
         color = "Condition", fill = "Condition") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "top",
          strip.text = element_text(face = "bold"))

ggsave(file.path(inputs$out_dir, "sim_data_overview_n20_hard_sem.pdf"),
       p_overview, width = 10, height = 5.7, device = cairo_pdf)
