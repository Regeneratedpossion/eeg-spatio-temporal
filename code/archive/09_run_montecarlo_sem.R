###########################################################
# 09_run_montecarlo_sem.R
#
# One Monte Carlo run:
# load one simulation, fit BGAM, compare methods, save results.
###########################################################

source(file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()), "code/utils.R"))
source(file.path(project_root(), "code/pipeline/fit_bgam.R"))
source(file.path(project_root(), "code/pipeline/compare_methods.R"))

# Parameters
# Seed comes from the SLURM array task ID.
seed <- as.integer(c(commandArgs(trailingOnly = TRUE), 1)[1])
ALPHA <- 0.1
basis <- basis_settings()

paths <- mc_paths()
files <- montecarlo_files(seed, paths)
dir.create(paths$compare_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(paths$summary_dir, showWarnings = FALSE, recursive = TRUE)

# Load dataset, fit BGAM, summarise posterior, compare methods, and save CSV.
sim_data     <- readRDS(files$sim)
ground_truth <- readRDS(files$gt)

fit <- fit_bgam(sim_data, seed = seed, k_xy = basis$k_xy, k_t = basis$k_t)

np    <- nuts_params(fit)
n_div <- sum(np$Value[np$Parameter == "divergent__"])

pilot_model <- aggregate_trials(sim_data)
elec        <- pilot_model %>% select(channel, xproj, yproj) %>% distinct()
channels    <- sort(unique(pilot_model$channel))
times       <- sort(unique(pilot_model$time))
cond_levels <- c("cond1", "cond2")

newd_full <- expand_grid(
    channel   = channels,
    time      = times,
    condition = factor(cond_levels, levels = cond_levels)
) %>%
    left_join(elec, by = "channel") %>%
    mutate(eeg_sem = 1)

epred_full <- add_epred_draws(newd_full, fit, re_formula = NA)
n_dr       <- length(unique(epred_full$.draw))
bgam_full  <- epred_full %>% ungroup() %>%
    select(channel, time, condition, .draw, .epred) %>%
    pivot_wider(names_from = condition, values_from = .epred) %>%
    mutate(diff = cond2 - cond1) %>%
    group_by(channel, time) %>%
    summarise(mean_diff = mean(diff),
              p_sig     = mean(diff > ROPE_UV),
              .groups   = "drop") %>%
    arrange(channel, time)

saveRDS(list(bgam_full   = bgam_full,
             ROPE_UV     = ROPE_UV,
             n_draws     = n_dr,
             seed        = seed,
             active_elec = ground_truth$active_electrodes),
        files$summary)

adj_chan <- readRDS(adjacency_path("adjacency_xproj_yproj_61ch.rds"))

result <- compare_methods_2d(
    sim_data     = sim_data,
    ground_truth = ground_truth,
    adj_chan     = adj_chan,
    bgam_full    = bgam_full,
    mne          = mne,
    alpha        = ALPHA
) %>%
    mutate(sim_id = seed,
           n_div  = n_div)

write.csv(result, files$csv, row.names = FALSE)
