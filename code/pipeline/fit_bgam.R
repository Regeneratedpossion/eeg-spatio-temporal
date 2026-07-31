###########################################################
# Fit spatiotemporal BGAM
#
###########################################################

source(file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()), "code/utils.R"))

setup_cmdstan_cache()

fit_bgam <- function(sim_data,
                     seed = 42,
                     k_xy = 16,
                     k_t  = 12,
                     bs_space = "ps",
                     bs_time = "ps",
                     full = FALSE,
                     adapt_delta = 0.80,
                     max_treedepth = 10,
                     chains = 4,
                     cores = 4,
                     iter = 3000,
                     warmup = 1000,
                     threads = 4) {

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
            "full = %s)"
        ),
        bs_space, bs_time, k_xy, k_t, as.character(full)
    )))

    bgam_fit <- brm(
        formula = bgam_formula,
        data = sim_data_model,
        family = gaussian(),
        chains = chains,
        cores = cores,
        iter = iter,
        warmup = warmup,
        seed = seed,
        threads = threading(threads = threads),
        backend = "cmdstanr",
        refresh = 100,
        control = list(adapt_delta = adapt_delta,
                       max_treedepth = max_treedepth)
    )

    bgam_fit
}
