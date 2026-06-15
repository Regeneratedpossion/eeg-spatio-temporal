###########################################################
# run_montecarlo.R
#
# One iteration of the Monte Carlo pipeline:
#   seed -> generate -> fit BGAM -> compare 5 methods -> save 1 csv
#
# Designed to be invoked as a SLURM array task:
#   Rscript code/run_montecarlo.R <seed>
# where <seed> is the SLURM_ARRAY_TASK_ID (1..N_SIMS).
#
# Output: results/comparisons/montecarlo/sim_<seed>.csv (one row per
# method, columns = bias_onset_ms, bias_offset_ms, etc.).
#
# Optional environment overrides:
#   MC_ALPHA=0.05
#   MC_DATA_SUBDIR=montecarlo_n20
#   MC_COMPARE_SUBDIR=montecarlo_n20
#   MC_SUMMARY_SUBDIR=montecarlo_n20
###########################################################

# ---------------------------------------------------------------
# Resolve seed from command line
# ---------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
seed <- if (length(args) >= 1) as.integer(args[[1]]) else 1L
stopifnot(!is.na(seed) && seed >= 1)

cat(sprintf("=== run_montecarlo: seed = %d ===\n", seed))

# ---------------------------------------------------------------
# Source helpers (utils, generators, fitter, MNE bridge)
# ---------------------------------------------------------------
source(file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()), "code/utils.R"))

suppressPackageStartupMessages({
    library(reticulate)
    library(purrr)
})

# Sourcing 01 / 03 only defines functions (fit_bgam(), compare_methods()).
# Each script's own driver block is guarded by sys.nframe() so nothing
# else fires. 00_generate_data.R is NOT sourced — datasets are pre-built
# locally and synced to the cluster (see code/00_generate_montecarlo.R).
source(file.path(project_root(), "code/01_fit_bgam.R"))
source(file.path(project_root(), "code/03_compare_methods.R"))

# 03 already imports MNE into a local `mne` object; if absent it sets
# it to NULL. Ensure logging is muted for run_montecarlo's own usage.
if (!is.null(mne)) mne$set_log_level("ERROR")

ALPHA <- as.numeric(Sys.getenv("MC_ALPHA", unset = "0.1"))
stopifnot(!is.na(ALPHA), ALPHA > 0, ALPHA < 1)
cat(sprintf("Using alpha = %.3f\n", ALPHA))

# Output directories
compare_subdir <- Sys.getenv("MC_COMPARE_SUBDIR", unset = "montecarlo")
data_subdir    <- Sys.getenv("MC_DATA_SUBDIR",    unset = "montecarlo")
summary_subdir <- Sys.getenv("MC_SUMMARY_SUBDIR", unset = "montecarlo")

out_dir   <- file.path(project_root(), "results/comparisons", compare_subdir)
data_dir  <- file.path(project_root(), "data", data_subdir)
sumry_dir <- file.path(project_root(), "results", summary_subdir)
for (d in c(out_dir, data_dir, sumry_dir)) {
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

sim_path  <- file.path(data_dir, sprintf("sim_%03d.rds", seed))
gt_path   <- file.path(data_dir, sprintf("gt_%03d.rds",  seed))
sumry_path<- file.path(sumry_dir, sprintf("bgam_summary_%03d.rds", seed))
csv_path  <- file.path(out_dir,  sprintf("sim_%03d.csv", seed))

# ---------------------------------------------------------------
# 1) Load pre-generated dataset from disk (built locally on the Mac
#    via code/00_generate_montecarlo.R and rsync'd to the cluster).
#    Bail out clearly if the seed's data is missing — the cluster
#    should not silently regenerate.
# ---------------------------------------------------------------
cat("Step 1/3: load dataset...\n")
if (!file.exists(sim_path) || !file.exists(gt_path)) {
    stop(sprintf(
        "Missing input data for seed %d:\n  %s\n  %s\nRun 00_generate_montecarlo.R with matching MC_DATA_SUBDIR first.",
        seed, sim_path, gt_path
    ))
}
sim_data     <- readRDS(sim_path)
ground_truth <- readRDS(gt_path)

true_on  <- ground_truth$true_onset
true_off <- ground_truth$true_offset
seed_ch  <- ground_truth$seed_electrode

# ---------------------------------------------------------------
# 2) Fit BGAM (slow step). The fit object itself is NOT saved
#    (~65 MB each is too large for 50 sims); instead, we save a
#    compact (time, channel)-level summary of posterior odds below.
# ---------------------------------------------------------------
cat("Step 2/3: fit BGAM...\n")
fit <- fit_bgam(sim_data, seed = seed)

# Quick convergence guard — bail out if the fit is unhealthy so the
# aggregator can drop this simulation cleanly.
np    <- nuts_params(fit)
n_div <- sum(np$Value[np$Parameter == "divergent__"])
if (n_div > 0) {
    cat(sprintf("WARNING: %d divergent transitions; marking sim as failed.\n",
                n_div))
}

# ---------------------------------------------------------------
# 3a) Save BGAM posterior summary on the full (channel x time) grid.
#     This lets us recompute odds / cluster decisions later under a
#     different ROPE or odds threshold without re-fitting.
# ---------------------------------------------------------------
cat("Step 3a/3: BGAM posterior summary on full grid...\n")
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
        sumry_path)

# ---------------------------------------------------------------
# 3b) Compare 6 methods on the full spatio-temporal grid.
#     onset/offset come from the largest connected component, so the
#     2D BGAM model is actually exercised in 2D (not collapsed to CZ).
#     compare_methods_2d() comes from 03_compare_methods.R.
# ---------------------------------------------------------------
cat("Step 3b/3: compare methods (2D)...\n")

# Load channel adjacency (61-channel template, ordered to match
# unique(sim_data$channel) inside compare_methods_2d).
adj_chan <- readRDS(adjacency_path("adjacency_xproj_yproj_61ch.rds"))

result <- compare_methods_2d(
    sim_data     = sim_data,
    ground_truth = ground_truth,
    adj_chan     = adj_chan,
    bgam_full    = bgam_full,       # reuse the Step 3a precomputation
    mne          = mne,
    alpha        = ALPHA
) %>%
    mutate(sim_id = seed,
           n_div  = n_div)

write.csv(result, csv_path, row.names = FALSE)
cat(sprintf("Saved:\n  %s\n  %s\n  %s\n  %s\n",
            sim_path, gt_path, sumry_path, csv_path))
