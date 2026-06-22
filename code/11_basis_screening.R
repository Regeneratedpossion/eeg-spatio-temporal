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

configs <- list(
    list(label = "A_k1212_tp", k = c(12L, 12L), bs = c("tp", "tp")),
    list(label = "B_k1212_ps", k = c(12L, 12L), bs = c("ps", "ps"))
)
cfg <- configs[[task_id]]
if (is.null(cfg)) {
    stop("Run as a Slurm array with SLURM_ARRAY_TASK_ID=1 or 2.")
}

# Runtime settings. Defaults reproduce the long-chain k=(12,12) baseline.
data_subdir <- Sys.getenv("SCREEN_DATA_SUBDIR", "montecarlo_p5_gaussian_hardlabel")
out_subdir <- Sys.getenv("SCREEN_OUT_SUBDIR", "basis_screening_full_ab")
seed <- as.integer(Sys.getenv("SCREEN_SEED", "1"))
iter <- as.integer(Sys.getenv("SCREEN_ITER", "3000"))
warmup <- as.integer(Sys.getenv("SCREEN_WARMUP", "1000"))
adapt_delta <- as.numeric(Sys.getenv("SCREEN_ADAPT_DELTA", "0.99"))
max_treedepth <- as.integer(Sys.getenv("SCREEN_MAX_TREEDEPTH", "13"))
refresh <- as.integer(Sys.getenv("SCREEN_REFRESH", "100"))
save_fit <- tolower(Sys.getenv("SCREEN_SAVE_FIT", "true")) %in%
    c("true", "1", "yes")

# Load and aggregate trial-level simulation data to participant x sensor x time cells.
out_dir <- file.path(project_root(), "results", out_subdir)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sim_data <- readRDS(data_path(data_subdir, sprintf("sim_%03d.rds", seed)))
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
    cfg$bs[1], cfg$bs[2], cfg$k[1], cfg$k[2]
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
                                               cfg$label, iter, warmup))
        saveRDS(fit, fit_file)
    }

    elapsed_sec <- as.numeric(Sys.time() - start_time, units = "secs")
    np <- nuts_params(fit)
    n_div <- sum(np$Value[np$Parameter == "divergent__"])
    max_td <- fit$fit@stan_args[[1]]$control$max_treedepth
    n_tree <- sum(np$Value[np$Parameter == "treedepth__"] >= max_td)
    summ <- posterior::summarise_draws(as_draws_df(fit))

    data.frame(
        label = cfg$label,
        status = "ok",
        data_subdir = data_subdir,
        out_subdir = out_subdir,
        seed = seed,
        k_xy = cfg$k[1],
        k_t = cfg$k[2],
        bs_space = cfg$bs[1],
        bs_time = cfg$bs[2],
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
        label = cfg$label,
        status = "error",
        data_subdir = data_subdir,
        out_subdir = out_subdir,
        seed = seed,
        k_xy = cfg$k[1],
        k_t = cfg$k[2],
        bs_space = cfg$bs[1],
        bs_time = cfg$bs[2],
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
                                       cfg$label, iter, warmup))
write.csv(result, out_file, row.names = FALSE)
print(result)
