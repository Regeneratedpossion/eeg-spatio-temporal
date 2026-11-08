###########################################################
# 05_aggregate_fixed_grid.R
#
# Aggregates the 100 per-seed fixed_grid CSVs produced by
# 05_extract_fixed_grid.R into three diagnostic plots for
# the single_node_cond<COND>_ds_cr_seeds001_100 batch:
#
#   1. bgam_fit_aggregate.png     - 6-channel fit curve
#        (mean across converged fits' posterior means, with a band
#        showing the 2.5-97.5% spread ACROSS SEEDS - i.e. Monte
#        Carlo recovery variability, not a within-seed CrI)
#   2. posterior_probability_aggregate.png - mean posterior
#        probability of exceeding SESOI, across converged fits
#   3. topomap_aggregate.png      - topo maps at 160/250/340ms,
#        averaged model estimate vs empirical vs truth
#
# Ground truth is identical across all 100 seeds (verified), so
# it is computed once from seed 1.
###########################################################

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(patchwork)
})
options(bitmapType = "cairo")

source(file.path(Sys.getenv("PROJECT_ROOT", unset = here::here()), "code/R/utils.R"))
source(file.path(project_root(), "code/R/plot_diagnostics.R"))
source(file.path(project_root(), "code/R/superpopulation_truth.R"))

base_dir <- file.path(project_root(), "results",
    fixed_grid_batch())
plot_label <- Sys.getenv("FIXED_GRID_LABEL", "cr basis")
out_dir <- file.path(base_dir, "figures_aggregate_passed")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cond_id <- sub(".*single_node_cond([A-Z])_.*", "\\1", base_dir)
man <- read.table(file.path(project_root(), "resources/passed_fits.txt"), col.names = c("cond", "seed"))
seed_dirs <- sprintf("seed_%03d", sort(man$seed[man$cond == cond_id]))
channels6 <- c("CZ", "FCZ", "CPZ", "C2", "PZ", "FP1")

# -----------------------------------------------------------
# Load and stack all per-seed CSVs
# -----------------------------------------------------------
read_seed_csv <- function(sub, file) {
    bind_rows(lapply(seed_dirs, function(sd) {
        p <- file.path(base_dir, sd, "fixed_grid", file)
        if (!file.exists(p)) return(NULL)
        read.csv(p, stringsAsFactors = FALSE)
    }))
}

six_all <- read_seed_csv(seed_dirs, "six_channel.csv")
topo_all <- read_seed_csv(seed_dirs, "topomap_3t.csv")
sesoi_all <- read_seed_csv(seed_dirs, "sesoi.csv")

stopifnot(length(unique(six_all$seed)) == length(seed_dirs))
stopifnot(length(unique(topo_all$seed)) == length(seed_dirs))

mean_sesoi <- mean(sesoi_all$global_sesoi)
cat("Mean global SESOI across converged fits:", mean_sesoi, "\n")
odds_threshold <- ODDS_THRESHOLD
probability_threshold <- odds_threshold / (1 + odds_threshold)

# -----------------------------------------------------------
# Ground truth / population truth (seed-invariant, use seed 1)
# -----------------------------------------------------------
config1 <- readRDS(file.path(base_dir, "seed_001", "config.rds"))
ground_truth <- readRDS(
    data_path(config1$data_subdir, "gt_001.rds")
)
population_truth <- paired_radius_superpopulation_truth(
    ground_truth,
    template_path = file.path(
        project_root(), "resources", "brms_meeg", "code", "erp_template.R"
    )
)
all_times <- sort(unique(six_all$time))
landmarks <- diagnostic_landmarks(ground_truth, all_times)
topomap_times <- topomap_diagnostic_times(all_times)

# -----------------------------------------------------------
# 1 & 2: six-channel fit curve + posterior probability
# Aggregate across converged fits' posterior point estimates.
# -----------------------------------------------------------
six_summary <- six_all %>%
    group_by(channel, time) %>%
    summarise(
        mean = mean(model_mean_diff),
        lo = quantile(model_mean_diff, 0.025),
        hi = quantile(model_mean_diff, 0.975),
        p_sig = mean(model_p_sig),
        empirical_diff = mean(empirical_diff),
        .groups = "drop"
    )

posterior_summary <- six_summary %>%
    select(channel, time, mean, lo, hi) %>%
    as.data.frame()

grand <- six_summary %>%
    transmute(channel, time, diff = empirical_diff) %>%
    as.data.frame()

# Empirical grand mean overlaps the superpopulation truth almost exactly
# at this sample size, so it is dropped from the plot to reduce clutter.
grand_empty <- data.frame(channel = character(), time = numeric(), diff = numeric())

# no participant-level traces retained from the extraction step
participant_empty <- data.frame(
    channel = character(), time = numeric(),
    diff = numeric(), participant = character()
)

truth6 <- population_truth %>%
    filter(channel %in% channels6) %>%
    arrange(match(channel, channels6), time) %>%
    as.data.frame()

fit_plot <- plot_bgam_fit(
    participant_empty, grand_empty, posterior_summary,
    truth6, channels6, landmarks
) + labs(
    title = paste0(plot_label, ": BGAM fit vs population-average truth (aggregate over converged fits)"),
    subtitle = NULL
) + scale_fill_manual(
    values = c("95% CrI" = "#ef8a62"),
    labels = "2.5-97.5% range of posterior means across seeds"
)

probability_data <- six_summary %>%
    transmute(channel, time, p_sig, above_threshold = p_sig > probability_threshold) %>%
    as.data.frame()
attr(probability_data, "sesoi") <- mean_sesoi

probability_plot <- plot_posterior_probability(
    probability_data, channels6, landmarks,
    probability_threshold, sesoi = mean_sesoi
) + labs(
    title = paste0(plot_label, ": mean posterior probability across converged fits"),
    subtitle = paste0(
        "Mean global SESOI (avg over converged fits) = ",
        sprintf("%.4f", mean_sesoi), " uV; ",
        "y = mean_seed P(cond2-cond1 > SESOI | data)"
    )
)

# -----------------------------------------------------------
# 3: topomap at 3 fixed times, averaged model vs empirical vs truth
# -----------------------------------------------------------
topo_summary <- topo_all %>%
    group_by(channel, time, xproj, yproj) %>%
    summarise(
        model_mean_diff = mean(model_mean_diff),
        empirical_diff = mean(empirical_diff),
        .groups = "drop"
    )

sensors <- topo_summary %>% distinct(channel, xproj, yproj)

truth_topo_raw <- population_truth %>%
    filter(time %in% topomap_times) %>%
    transmute(channel, time, truth = effect) %>%
    left_join(sensors, by = "channel")

topo_wide <- topo_summary %>%
    left_join(truth_topo_raw %>% select(channel, time, truth), by = c("channel", "time"))

truth_topo <- topo_wide %>%
    transmute(channel, time, xproj, yproj, value = truth,
              layer = "1. Ground truth (population average)")
empirical_topo <- topo_wide %>%
    transmute(channel, time, xproj, yproj, value = empirical_diff,
              layer = "2. Empirical difference")
residual_truth_topo <- topo_wide %>%
    transmute(channel, time, xproj, yproj, value = empirical_diff - truth,
              layer = "3. Noise residual (Empirical - Ground truth)")
model_topo <- topo_wide %>%
    transmute(channel, time, xproj, yproj, value = model_mean_diff,
              layer = "4. BGAM estimate")
residual_model_topo <- topo_wide %>%
    transmute(channel, time, xproj, yproj, value = empirical_diff - model_mean_diff,
              layer = "5. Empirical - BGAM")

electrode_data <- bind_rows(
    truth_topo, empirical_topo, residual_truth_topo,
    model_topo, residual_model_topo
) %>%
    mutate(
        active = channel %in% ground_truth$active_electrodes,
        time_label = factor(
            sprintf("%.0f ms", time * 1000),
            levels = sprintf("%.0f ms", topomap_times * 1000)
        )
    )

topomap_plot_base <- plot_topomap_diagnostics(electrode_data, topomap_times)
shared_limit <- attr(topomap_plot_base, "topomap_colour_limit")
topomap_plot <- topomap_plot_base + patchwork::plot_annotation(
    title = paste0(plot_label, ": spatial diagnostics, aggregate over converged fits"),
    subtitle = paste(
        "Columns:", paste(sprintf("%.0f ms", topomap_times * 1000), collapse = " / "),
        "- all panels share one symmetric colour scale."
    ),
    theme = theme(
        plot.title = element_text(face = "bold", size = 13),
        plot.subtitle = element_text(size = 9)
    )
)

# -----------------------------------------------------------
# Save
# -----------------------------------------------------------
ggsave(file.path(out_dir, "bgam_fit_aggregate.png"), fit_plot,
       width = 11, height = 7.5, dpi = 300)
ggsave(file.path(out_dir, "posterior_probability_aggregate.png"), probability_plot,
       width = 11, height = 7, dpi = 300)
ggsave(file.path(out_dir, "topomap_aggregate.png"), topomap_plot,
       width = 11.5, height = 13, dpi = 300)
ggsave(file.path(out_dir, "bgam_fit_aggregate.pdf"), fit_plot,
       width = 11, height = 7.5)
ggsave(file.path(out_dir, "posterior_probability_aggregate.pdf"), probability_plot,
       width = 11, height = 7)
ggsave(file.path(out_dir, "topomap_aggregate.pdf"), topomap_plot,
       width = 11.5, height = 13)

saveRDS(
    list(
        six_summary = six_summary,
        topo_summary = topo_summary,
        mean_sesoi = mean_sesoi,
        population_truth = population_truth,
        landmarks = landmarks,
        topomap_times = topomap_times
    ),
    file.path(out_dir, "aggregate_plot_data.rds")
)

cat("Saved:\n")
cat(" ", file.path(out_dir, "bgam_fit_aggregate.png"), "\n")
cat(" ", file.path(out_dir, "posterior_probability_aggregate.png"), "\n")
cat(" ", file.path(out_dir, "topomap_aggregate.png"), "\n")
