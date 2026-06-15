###########################################################
# Fit spatiotemporal BGAM
#
###########################################################

source(file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()), "code/utils.R"))

setup_cmdstan_cache()

fit_bgam <- function(sim_data,
                     seed = 42,
                     k_xy = 20,
                     k_t  = 20,
                     full = FALSE) {

    sim_data_model <- aggregate_trials(sim_data)
    contrasts(sim_data_model$condition) <- c(-0.5, 0.5)

    bgam_fit <- brm(
        formula = eeg_mean | se(eeg_sem) ~
        condition +
        (1 + condition | participant) +
        t2(xproj, yproj, time,
        by = condition,
        bs = c("tp", "tp"),
        k = c(k_xy, k_t),
        d = c(2, 1),
        full = full),
        data = sim_data_model,
        family = gaussian(),
        chains = 4,
        cores = 4,
        iter = 3000,
        warmup = 1000,
        seed = seed,
        threads = threading(threads = 4),
        backend = "cmdstanr",
        control = list(adapt_delta = 0.99, max_treedepth = 13)
    )

    bgam_fit
}
