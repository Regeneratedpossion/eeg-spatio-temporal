###########################################################
# 24_basis_screening_full.R
#
# Full A/B basis screening with configurable sampler control and sigma=TRUE.
# Used for sensitivity checks requested by Paul:
#   adapt_delta = 0.80, max_treedepth = 10
#   eeg_mean | se(eeg_sem, sigma = TRUE)
###########################################################

source(file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()), "code/utils.R"))

setup_cmdstan_cache()

task_id <- as.integer(Sys.getenv(
    "SLURM_ARRAY_TASK_ID",
    unset = as.character(as.integer(c(commandArgs(trailingOnly = TRUE), 1)[1]))
))

configs <- data.frame(
    task_id = 1:2,
    label_prefix = c("A", "B"),
    basis_suffix = c("tp", "ps"),
    bs_space = c("tp", "ps"),
    bs_time = c("tp", "ps")
)

cfg <- configs[configs$task_id == task_id, ]
stopifnot(nrow(cfg) == 1)

data_subdir <- Sys.getenv("SCREEN_DATA_SUBDIR",
                          "montecarlo_p5_gaussian_hardlabel")
seed <- as.integer(Sys.getenv("SCREEN_SEED", "1"))
iter <- as.integer(Sys.getenv("SCREEN_ITER", "3000"))
warmup <- as.integer(Sys.getenv("SCREEN_WARMUP", "1000"))
k_xy <- as.integer(Sys.getenv("SCREEN_K_XY", "12"))
k_t <- as.integer(Sys.getenv("SCREEN_K_T", "12"))
adapt_delta <- as.numeric(Sys.getenv("SCREEN_ADAPT_DELTA", "0.8"))
max_treedepth <- as.integer(Sys.getenv("SCREEN_MAX_TREEDEPTH", "10"))
out_subdir <- Sys.getenv("SCREEN_OUT_SUBDIR",
                         "basis_screening_full_ab_sigma_true_default_control")

cfg$k_xy <- k_xy
cfg$k_t <- k_t
cfg$label <- sprintf("%s_k%02d%02d_%s",
                     cfg$label_prefix, cfg$k_xy, cfg$k_t, cfg$basis_suffix)

adapt_tag <- gsub("\\.", "p", as.character(adapt_delta))
label_tag <- sprintf("%s_%s_adapt%s_td%d",
                     basename(data_subdir), cfg$label, adapt_tag, max_treedepth)

out_dir <- file.path(project_root(), "results", out_subdir)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

sim_data <- readRDS(data_path(data_subdir, sprintf("sim_%03d.rds", seed)))
sim_data_model <- aggregate_trials(sim_data)
contrasts(sim_data_model$condition) <- c(-0.5, 0.5)

bgam_formula <- bf(as.formula(sprintf(
    paste(
        "eeg_mean | se(eeg_sem, sigma = TRUE) ~",
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
                                           label_tag, iter, warmup))
    saveRDS(fit, fit_file)

    elapsed_sec <- as.numeric(Sys.time() - start_time, units = "secs")
    np <- nuts_params(fit)
    n_div <- sum(np$Value[np$Parameter == "divergent__"])
    n_tree <- sum(np$Value[np$Parameter == "treedepth__"] >= max_treedepth)
    summ <- posterior::summarise_draws(as_draws_df(fit))

    data.frame(
        label = cfg$label,
        label_tag = label_tag,
        status = "ok",
        data_subdir = data_subdir,
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
        label_tag = label_tag,
        status = "error",
        data_subdir = data_subdir,
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
                                       label_tag, iter, warmup))
write.csv(result, out_file, row.names = FALSE)
print(result)
