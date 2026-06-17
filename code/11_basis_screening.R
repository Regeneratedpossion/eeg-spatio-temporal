###########################################################
# 11_basis_screening.R
#
# Short feasibility screen for basis size and basis type.
# Uses P=5, seed=1, SEM, Gaussian signal + hard-threshold labels.
###########################################################

source(file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()), "code/utils.R"))

setup_cmdstan_cache()

task_id <- as.integer(Sys.getenv(
    "SLURM_ARRAY_TASK_ID",
    unset = as.character(as.integer(c(commandArgs(trailingOnly = TRUE), 1)[1]))
))

configs <- data.frame(
    task_id = 1:4,
    label = c("A_k1212_tp", "B_k1212_ps", "C_k2020_tp", "D_k2020_ps"),
    k_xy = c(12L, 12L, 20L, 20L),
    k_t = c(12L, 12L, 20L, 20L),
    bs_space = c("tp", "ps", "tp", "ps"),
    bs_time = c("tp", "ps", "tp", "ps")
)

cfg <- configs[configs$task_id == task_id, ]
stopifnot(nrow(cfg) == 1)

data_subdir <- Sys.getenv("SCREEN_DATA_SUBDIR",
                          "montecarlo_p5_gaussian_hardlabel")
seed <- as.integer(Sys.getenv("SCREEN_SEED", "1"))
iter <- as.integer(Sys.getenv("SCREEN_ITER", "200"))
warmup <- as.integer(Sys.getenv("SCREEN_WARMUP", "100"))

out_dir <- file.path(project_root(), "results", "basis_screening")
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
    cfg$bs_space, cfg$bs_time, cfg$k_xy, cfg$k_t
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
        refresh = 10,
        control = list(adapt_delta = 0.99, max_treedepth = 13)
    )

    elapsed_sec <- as.numeric(Sys.time() - start_time, units = "secs")
    np <- nuts_params(fit)
    n_div <- sum(np$Value[np$Parameter == "divergent__"])
    max_td <- fit$fit@stan_args[[1]]$control$max_treedepth
    n_tree <- sum(np$Value[np$Parameter == "treedepth__"] >= max_td)
    summ <- posterior::summarise_draws(as_draws_df(fit))

    data.frame(
        label = cfg$label,
        status = "ok",
        seed = seed,
        k_xy = cfg$k_xy,
        k_t = cfg$k_t,
        bs_space = cfg$bs_space,
        bs_time = cfg$bs_time,
        iter = iter,
        warmup = warmup,
        elapsed_sec = elapsed_sec,
        n_divergent = n_div,
        n_treedepth = n_tree,
        max_rhat = max(summ$rhat, na.rm = TRUE),
        min_ess_bulk = min(summ$ess_bulk, na.rm = TRUE),
        mean_ess_bulk = mean(summ$ess_bulk, na.rm = TRUE),
        error = NA_character_
    )
}, error = function(e) {
    data.frame(
        label = cfg$label,
        status = "error",
        seed = seed,
        k_xy = cfg$k_xy,
        k_t = cfg$k_t,
        bs_space = cfg$bs_space,
        bs_time = cfg$bs_time,
        iter = iter,
        warmup = warmup,
        elapsed_sec = as.numeric(Sys.time() - start_time, units = "secs"),
        n_divergent = NA_integer_,
        n_treedepth = NA_integer_,
        max_rhat = NA_real_,
        min_ess_bulk = NA_real_,
        mean_ess_bulk = NA_real_,
        error = conditionMessage(e)
    )
})

out_file <- file.path(out_dir, sprintf("basis_screen_%s.csv", cfg$label))
write.csv(result, out_file, row.names = FALSE)
print(result)
