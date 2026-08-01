###########################################################
# 23_basis_screening.R
#
# Full-BGAM k=(12,12) baseline models:
#   A_k1212_tp = tensor-product thin-plate basis
#   B_k1212_ps = tensor-product P-spline basis
###########################################################

source(file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()), "code/utils.R"))

setup_cmdstan_cache()

# Select one of the two predefined full-BGAM baselines.
# Slurm array task 1 -> A_k1212_tp; task 2 -> B_k1212_ps.
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
model_label <- c("A_k1212_tp", "B_k1212_ps")[task_id]
space_basis <- c("tp", "ps")[task_id]
time_basis <- space_basis

# Long-chain baseline settings.
data_subdir <- "montecarlo_p5_gaussian_hardlabel"
out_subdir <- "basis_screening_full_ab"
seed <- 1L
iter <- 3000L
warmup <- 1000L
space_k <- 12L
time_k <- 12L
adapt_delta <- 0.99
max_treedepth <- 13L
refresh <- 100L

# Output folder.
out_dir <- file.path(project_root(), "results", out_subdir)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Simulation data.
sim_file <- data_path(data_subdir, "sim_001.rds")
sim_data <- readRDS(sim_file)

# Participant x sensor x time means and SEMs.
sim_data_model <- aggregate_trials(sim_data)
contrasts(sim_data_model$condition) <- c(-0.5, 0.5)

start_time <- Sys.time()

fit <- brm(
    formula = bf(
        eeg_mean | se(eeg_sem) ~
            condition +
            (1 + condition | participant) +
            t2(xproj, yproj, time,
               by = condition,
               bs = c(space_basis, time_basis),
               k = c(space_k, time_k),
               d = c(2, 1),
               full = FALSE)
    ),
    data = sim_data_model,
    family = gaussian(),
    chains = 4,
    cores = 4,
    iter = iter,
    warmup = warmup,
    seed = seed,
    threads = threading(threads = 4),
    backend = "cmdstanr",
    refresh = refresh,
    control = list(adapt_delta = adapt_delta,
                   max_treedepth = max_treedepth)
)

elapsed_sec <- as.numeric(Sys.time() - start_time, units = "secs")
nuts <- nuts_params(fit)
draw_summary <- posterior::summarise_draws(as_draws_df(fit))

saveRDS(fit, file.path(out_dir, paste0("fit_", model_label, ".rds")))

diagnostics <- data.frame(
    label = model_label,
    elapsed_hours = round(elapsed_sec / 3600, 2),
    n_divergent = sum(nuts$Value[nuts$Parameter == "divergent__"]),
    n_treedepth = sum(nuts$Value[nuts$Parameter == "treedepth__"] >= max_treedepth),
    max_rhat = max(draw_summary$rhat, na.rm = TRUE),
    min_ess_bulk = min(draw_summary$ess_bulk, na.rm = TRUE),
    min_ess_tail = min(draw_summary$ess_tail, na.rm = TRUE)
)

write.csv(diagnostics,
          file.path(out_dir, paste0("summary_", model_label, ".csv")),
          row.names = FALSE)
print(diagnostics)
