###########################################################
# 11_basis_screening.R
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

if (task_id == 1) {
    model_label <- "A_k1212_tp"
    space_basis <- "tp"
    time_basis <- "tp"
} else if (task_id == 2) {
    model_label <- "B_k1212_ps"
    space_basis <- "ps"
    time_basis <- "ps"
} else {
    stop("Run as a Slurm array with SLURM_ARRAY_TASK_ID=1 or 2.")
}

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
save_fit <- TRUE

# Output folder.
out_dir <- file.path(project_root(), "results", out_subdir)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Simulation data.
sim_file <- data_path(data_subdir, "sim_001.rds")
sim_data <- readRDS(sim_file)

# Participant x sensor x time means and SEMs.
sim_data_model <- aggregate_trials(sim_data)
contrasts(sim_data_model$condition) <- c(-0.5, 0.5)

bgam_formula <- bf(as.formula(sprintf(
    paste(
        "eeg_mean | se(eeg_sem) ~",
        "condition +",
        "(1 + condition | participant) +",
        "t2(xproj, yproj, time,",
        "by = condition,",
        "bs = c('%s', '%s'),",
        "k = c(%d, %d),",
        "d = c(2, 1),",
        "full = FALSE)"
    ),
    space_basis, time_basis, space_k, time_k
)))

start_time <- Sys.time()

result <- tryCatch({
    fit <- brm(
        formula = bgam_formula,
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

    fit_file <- NA_character_
    if (save_fit) {
        fit_file <- file.path(out_dir, sprintf("fit_%s_iter%d_warmup%d.rds",
                                               model_label, iter, warmup))
        saveRDS(fit, fit_file)
    }

    elapsed_sec <- as.numeric(Sys.time() - start_time, units = "secs")
    np <- nuts_params(fit)
    n_div <- sum(np$Value[np$Parameter == "divergent__"])
    max_td <- fit$fit@stan_args[[1]]$control$max_treedepth
    n_tree <- sum(np$Value[np$Parameter == "treedepth__"] >= max_td)
    summ <- posterior::summarise_draws(as_draws_df(fit))

    data.frame(
        label = model_label,
        status = "ok",
        data_subdir = data_subdir,
        out_subdir = out_subdir,
        seed = seed,
        k_xy = space_k,
        k_t = time_k,
        bs_space = space_basis,
        bs_time = time_basis,
        iter = iter,
        warmup = warmup,
        adapt_delta = adapt_delta,
        max_treedepth = max_treedepth,
        elapsed_sec = elapsed_sec,
        n_divergent = n_div,
        n_treedepth = n_tree,
        max_rhat = max(summ$rhat, na.rm = TRUE),
        n_rhat_gt_1_01 = sum(summ$rhat > 1.01, na.rm = TRUE),
        min_ess_bulk = min(summ$ess_bulk, na.rm = TRUE),
        mean_ess_bulk = mean(summ$ess_bulk, na.rm = TRUE),
        min_ess_tail = min(summ$ess_tail, na.rm = TRUE),
        fit_file = fit_file,
        error = NA_character_
    )
}, error = function(e) {
    data.frame(
        label = model_label,
        status = "error",
        data_subdir = data_subdir,
        out_subdir = out_subdir,
        seed = seed,
        k_xy = space_k,
        k_t = time_k,
        bs_space = space_basis,
        bs_time = time_basis,
        iter = iter,
        warmup = warmup,
        adapt_delta = adapt_delta,
        max_treedepth = max_treedepth,
        elapsed_sec = as.numeric(Sys.time() - start_time, units = "secs"),
        n_divergent = NA_integer_,
        n_treedepth = NA_integer_,
        max_rhat = NA_real_,
        n_rhat_gt_1_01 = NA_integer_,
        min_ess_bulk = NA_real_,
        mean_ess_bulk = NA_real_,
        min_ess_tail = NA_real_,
        fit_file = NA_character_,
        error = conditionMessage(e)
    )
})

out_file <- file.path(out_dir, sprintf("basis_screen_%s_iter%d_warmup%d.csv",
                                       model_label, iter, warmup))
write.csv(result, out_file, row.names = FALSE)
print(result)
