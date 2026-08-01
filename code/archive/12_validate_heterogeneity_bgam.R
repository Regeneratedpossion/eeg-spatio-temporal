###########################################################
# Validate one BGAM fit under one heterogeneity scenario.
#
# Outputs:
#   - complete brms fit
#   - posterior effect map on the channel x time grid
#   - convergence and recovery summary
#   - posterior-probability heatmap
###########################################################

args <- commandArgs(trailingOnly = TRUE)
scenario <- if (length(args) >= 1L) args[[1]] else "low"
seed <- if (length(args) >= 2L) as.integer(args[[2]]) else 1L

Sys.setenv(MC_SCENARIO = scenario)
array_job_id <- Sys.getenv(
    "SLURM_ARRAY_JOB_ID",
    unset = Sys.getenv("SLURM_JOB_ID", unset = "local")
)
array_task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = scenario)
Sys.setenv(
    CMDSTAN_CACHE_DIR = file.path(
        Sys.getenv("PROJECT_ROOT", unset = getwd()),
        "cmdstan_cache",
        paste0("heterogeneity_", array_job_id, "_", array_task_id)
    )
)

source(file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()), "code/utils.R"))
source(file.path(project_root(), "code/03_fit_bgam.R"))
source(file.path(project_root(), "code/04_compare_methods.R"))

library(ggplot2)

basis <- basis_settings()
bs_space <- Sys.getenv("VALIDATION_BS_SPACE", unset = "ps")
bs_time <- Sys.getenv("VALIDATION_BS_TIME", unset = "ps")
model_tag <- sprintf("%s_%s_k%d%d", bs_space, bs_time,
                     basis$k_xy, basis$k_t)
max_treedepth <- as.integer(Sys.getenv("VALIDATION_MAX_TREEDEPTH", "10"))
adapt_delta <- as.numeric(Sys.getenv("VALIDATION_ADAPT_DELTA", "0.80"))

paths <- mc_paths()
files <- montecarlo_files(seed, paths)
out_dir <- file.path(project_root(), "results", "heterogeneity_bgam_validation",
                     model_tag, scenario)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(files$sim) || !file.exists(files$gt)) {
    stop("Missing simulation or ground truth for scenario ", scenario,
         ", seed ", seed, ": ", files$sim)
}

sim_data <- readRDS(files$sim)
ground_truth <- readRDS(files$gt)
if (!identical(ground_truth$scenario, scenario)) {
    stop("Ground-truth scenario does not match requested scenario.")
}

start_time <- Sys.time()
fit <- fit_bgam(
    sim_data,
    seed = seed,
    k_xy = basis$k_xy,
    k_t = basis$k_t,
    bs_space = bs_space,
    bs_time = bs_time,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth
)
elapsed_hours <- as.numeric(Sys.time() - start_time, units = "hours")

fit_file <- file.path(
    out_dir,
    sprintf("fit_%s_%s_seed%03d.rds", model_tag, scenario, seed)
)
saveRDS(fit, fit_file)

np <- nuts_params(fit)
draw_summary <- posterior::summarise_draws(as_draws_df(fit))
n_divergent <- sum(np$Value[np$Parameter == "divergent__"])
n_treedepth <- sum(np$Value[np$Parameter == "treedepth__"] >= max_treedepth)
max_rhat <- max(draw_summary$rhat, na.rm = TRUE)
min_ess_bulk <- min(draw_summary$ess_bulk, na.rm = TRUE)
min_ess_tail <- min(draw_summary$ess_tail, na.rm = TRUE)

model_data <- aggregate_trials(sim_data)
electrodes <- model_data %>%
    select(channel, xproj, yproj) %>%
    distinct() %>%
    arrange(channel)
channels <- electrodes$channel
times <- sort(unique(model_data$time))
condition_levels <- c("cond1", "cond2")

prediction_grid <- expand_grid(channel = channels, time = times) %>%
    left_join(electrodes, by = "channel") %>%
    mutate(eeg_sem = 1)
newdata_c1 <- prediction_grid %>%
    mutate(condition = factor("cond1", levels = condition_levels))
newdata_c2 <- prediction_grid %>%
    mutate(condition = factor("cond2", levels = condition_levels))

epred_c1 <- posterior_epred(fit, newdata = newdata_c1, re_formula = NA)
epred_c2 <- posterior_epred(fit, newdata = newdata_c2, re_formula = NA)
effect_draws <- epred_c2 - epred_c1

posterior_map <- prediction_grid %>%
    transmute(
        channel,
        time,
        mean_diff = colMeans(effect_draws),
        p_sig = colMeans(effect_draws > ROPE_UV)
    ) %>%
    mutate(
        odds = pmax(p_sig, 1 / nrow(effect_draws)) /
            pmax(1 - p_sig, 1 / nrow(effect_draws))
    )

odds_matrix <- posterior_map %>%
    select(channel, time, odds) %>%
    pivot_wider(names_from = channel, values_from = odds) %>%
    arrange(time) %>%
    select(all_of(channels)) %>%
    as.matrix()

adjacency <- readRDS(adjacency_path("adjacency_xproj_yproj_61ch.rds"))
adjacency <- adjacency[channels, channels]
predicted_mask <- .filter_components(
    odds_mask(odds_matrix, threshold = ODDS_THRESHOLD),
    adjacency,
    min_size = MIN_CLUSTER_SIZE
)
timing <- .extract_onset_offset_2d(
    predicted_mask,
    times,
    adjacency,
    min_size = MIN_CLUSTER_SIZE
)

truth_mask <- outer(
    times %in% ground_truth$temporal_active_times,
    channels %in% ground_truth$active_electrodes,
    "&"
)
tp <- sum(predicted_mask & truth_mask)
fp <- sum(predicted_mask & !truth_mask)
fn <- sum(!predicted_mask & truth_mask)
iou <- if ((tp + fp + fn) > 0) tp / (tp + fp + fn) else NA_real_

convergence_pass <- n_divergent == 0L && max_rhat <= 1.01 &&
    min_ess_bulk >= 400 && min_ess_tail >= 400

validation_summary <- data.frame(
    scenario = scenario,
    seed = seed,
    n_participants = length(unique(sim_data$participant)),
    k_xy = basis$k_xy,
    k_t = basis$k_t,
    bs_space = bs_space,
    bs_time = bs_time,
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    elapsed_hours = elapsed_hours,
    n_divergent = n_divergent,
    n_treedepth = n_treedepth,
    max_rhat = max_rhat,
    min_ess_bulk = min_ess_bulk,
    min_ess_tail = min_ess_tail,
    convergence_pass = convergence_pass,
    true_onset = ground_truth$true_onset,
    estimated_onset = timing$onset,
    onset_error_ms = (timing$onset - ground_truth$true_onset) * 1000,
    true_offset = ground_truth$true_offset,
    estimated_offset = timing$offset,
    offset_error_ms = (timing$offset - ground_truth$true_offset) * 1000,
    tp = tp,
    fp = fp,
    fn = fn,
    iou = iou,
    stringsAsFactors = FALSE
)

map_file <- file.path(out_dir,
                      sprintf("posterior_map_%s_%s_seed%03d.rds",
                              model_tag, scenario, seed))
summary_file <- file.path(
    out_dir,
    sprintf("validation_summary_%s_%s_seed%03d.csv",
            model_tag, scenario, seed)
)
saveRDS(posterior_map, map_file)
write.csv(validation_summary, summary_file, row.names = FALSE)

plot_data <- posterior_map %>%
    mutate(channel = factor(channel, levels = rev(channels)))
probability_plot <- ggplot(plot_data, aes(time, channel, fill = p_sig)) +
    geom_tile() +
    geom_vline(xintercept = c(ground_truth$true_onset,
                              ground_truth$true_offset),
               linetype = "dashed", colour = "white") +
    scale_fill_viridis_c(limits = c(0, 1)) +
    labs(x = "Time (s)", y = "Electrode",
         fill = "Pr(effect > 0.1)",
         title = paste("BGAM posterior map:", model_tag, scenario)) +
    theme_minimal(base_size = 9)
ggsave(
    file.path(out_dir,
              sprintf("posterior_map_%s_%s_seed%03d.pdf",
                      model_tag, scenario, seed)),
    probability_plot,
    width = 9,
    height = 11
)

print(validation_summary)
