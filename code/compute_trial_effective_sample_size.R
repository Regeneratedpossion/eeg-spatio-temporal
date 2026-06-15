###########################################################
# compute_trial_effective_sample_size.R
#
# Estimate the effective number of independent trials (R_eff)
# from empirical trial-to-trial correlations.
#
# For each participant x channel x condition group, trials are treated
# as repeated time series. We compute all pairwise trial correlations
# within the requested time/channel region, average correlations with
# Fisher z, then use:
#
#   R_eff = R / (1 + (R - 1) * rho)
#
# where R is the nominal number of trials and rho is the average
# trial-to-trial correlation.
###########################################################

library(dplyr)

source(file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()), "code/utils.R"))

dataset <- Sys.getenv("DATASET", "n20_hard")

files <- switch(dataset,
    n20_hard = list(
        sim = data_path("montecarlo_n20_hard", "sim_001.rds"),
        gt = data_path("montecarlo_n20_hard", "gt_001.rds"),
        out = diag_path("trial_effective_sample_size_n20_hard.csv")
    ),
    pilot = list(
        sim = data_path("sim_data_pilot.rds"),
        gt = data_path("ground_truth_pilot.rds"),
        out = diag_path("trial_effective_sample_size_pilot.csv")
    ),
    stop("Unknown DATASET: ", dataset, " (use n20_hard or pilot)")
)

sim_data <- readRDS(files$sim)
ground_truth <- readRDS(files$gt)
dir.create(dirname(files$out), showWarnings = FALSE, recursive = TRUE)

trial_correlations <- function(dat) {
    wide <- dat %>%
        select(trial, time, eeg) %>%
        tidyr::pivot_wider(names_from = time, values_from = eeg) %>%
        arrange(trial)

    mat <- as.matrix(wide %>% select(-trial))
    if (nrow(mat) < 2 || ncol(mat) < 2) return(numeric())

    cm <- suppressWarnings(cor(t(mat), use = "pairwise.complete.obs"))
    cm[upper.tri(cm)]
}

summarise_region <- function(region_name, dat) {
    cors <- dat %>%
        group_by(participant, channel, condition) %>%
        group_modify(~ data.frame(rho = trial_correlations(.x))) %>%
        ungroup() %>%
        filter(is.finite(rho), abs(rho) < 1)

    R <- dat %>%
        summarise(R = max(dplyr::n_distinct(trial)), .by = c(participant, channel, condition)) %>%
        pull(R) %>%
        max()

    mean_rho <- tanh(mean(atanh(cors$rho)))
    r_eff <- R / (1 + (R - 1) * mean_rho)
    r_eff <- min(R, r_eff)

    data.frame(
        region = region_name,
        n_groups = dplyr::n_distinct(
            paste(cors$participant, cors$channel, cors$condition, sep = "_")
        ),
        n_pairwise_corr = nrow(cors),
        R = R,
        mean_rho_fisher = mean_rho,
        median_rho = median(cors$rho),
        q25_rho = unname(quantile(cors$rho, 0.25)),
        q75_rho = unname(quantile(cors$rho, 0.75)),
        R_eff = r_eff
    )
}

baseline <- sim_data %>%
    filter(time < ground_truth$true_onset)

active_window <- sim_data %>%
    filter(time >= ground_truth$true_onset,
           time <= ground_truth$true_offset)

active_window_active_channels <- active_window %>%
    filter(channel %in% ground_truth$active_electrodes)

active_window_seed <- active_window %>%
    filter(channel == ground_truth$seed_electrode)

out <- bind_rows(
    summarise_region("baseline_all_channels", baseline),
    summarise_region("active_window_active_channels", active_window_active_channels),
    summarise_region("active_window_seed_CZ", active_window_seed)
)

write.csv(out, files$out, row.names = FALSE)
print(out)
cat("Saved:", files$out, "\n")
