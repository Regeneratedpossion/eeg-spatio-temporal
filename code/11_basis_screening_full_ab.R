###########################################################
# Full diagnostic rerun for basis screening A/B.
# This is the single full-BGAM k=(12,12) A/B script.
#
# Defaults reproduce the original full run:
#   SCREEN_DATA_SUBDIR=montecarlo_p5_gaussian_hardlabel
#   SCREEN_OUT_SUBDIR=basis_screening_full_ab
#   SCREEN_ADAPT_DELTA=0.99
#   SCREEN_MAX_TREEDEPTH=13
#
# The same script can run sampler-control sensitivity checks by
# overriding these environment variables.
###########################################################

source(file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()), "code/utils.R"))

setup_cmdstan_cache()

task_id <- as.integer(Sys.getenv(
    "SLURM_ARRAY_TASK_ID",
    unset = as.character(as.integer(c(commandArgs(trailingOnly = TRUE), 1)[1]))
))

configs <- data.frame(
    task_id = 1:2,
    label = c("A_k1212_tp", "B_k1212_ps"),
    k_xy = c(12L, 12L),
    k_t = c(12L, 12L),
    bs_space = c("tp", "ps"),
    bs_time = c("tp", "ps")
)

cfg <- configs[configs$task_id == task_id, ]
stopifnot(nrow(cfg) == 1)

data_subdir <- Sys.getenv("SCREEN_DATA_SUBDIR",
                          "montecarlo_p5_gaussian_hardlabel")
out_subdir <- Sys.getenv("SCREEN_OUT_SUBDIR",
                         "basis_screening_full_ab")
seed <- as.integer(Sys.getenv("SCREEN_SEED", "1"))
iter <- as.integer(Sys.getenv("SCREEN_ITER", "3000"))
warmup <- as.integer(Sys.getenv("SCREEN_WARMUP", "1000"))
adapt_delta <- as.numeric(Sys.getenv("SCREEN_ADAPT_DELTA", "0.99"))
max_treedepth <- as.integer(Sys.getenv("SCREEN_MAX_TREEDEPTH", "13"))

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
        refresh = 100,
        control = list(adapt_delta = adapt_delta,
                       max_treedepth = max_treedepth)
    )

    fit_file <- file.path(out_dir, sprintf("fit_%s_iter%d_warmup%d.rds",
                                           cfg$label, iter, warmup))
    saveRDS(fit, fit_file)

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
        k_xy = cfg$k_xy,
        k_t = cfg$k_t,
        bs_space = cfg$bs_space,
        bs_time = cfg$bs_time,
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
        k_xy = cfg$k_xy,
        k_t = cfg$k_t,
        bs_space = cfg$bs_space,
        bs_time = cfg$bs_time,
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

out_file <- file.path(out_dir, sprintf("basis_screen_full_%s_iter%d_warmup%d.csv",
                                       cfg$label, iter, warmup))
write.csv(result, out_file, row.names = FALSE)
print(result)
