###############################################################
# 17_prior_predictive_check.R
#
# Prior predictive check for hard-threshold + SEM BGAM.
###############################################################

library(brms)
library(ggplot2)

source(file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()), "code/utils.R"))

inputs <- ppc_inputs()
basis <- basis_settings()
dir.create(inputs$fig_dir, showWarnings = FALSE, recursive = TRUE)

sim_data_model <- aggregate_trials(readRDS(inputs$sim_file))
contrasts(sim_data_model$condition) <- c(-0.5, 0.5)

# brms needs a proper fixed-effect prior when sample_prior = "only".
prior_fit <- brm(
    formula = eeg_mean | se(eeg_sem) ~
        condition +
        (1 + condition | participant) +
        t2(xproj, yproj, time,
           by = condition,
           bs = c("tp", "tp"),
           k = c(basis$k_xy, basis$k_t),
           d = c(2, 1),
           full = basis$full),
    data = sim_data_model,
    family = gaussian(),
    prior = prior(normal(0, 5), class = b),
    sample_prior = "only",
    chains = 2,
    cores = 2,
    iter = 1000,
    backend = "cmdstanr",
    seed = 42
)

capture.output(prior_summary(prior_fit),
               file = file.path(inputs$fig_dir, "ppc_prior_summary.txt"))

save_ppc_plot <- function(filename, title, type = NULL, stat = NULL) {
    p <- pp_check(prior_fit, type = type, stat = stat, ndraws = 100) +
        ggtitle(title) +
        theme_minimal()

    ggsave(file.path(inputs$fig_dir, filename),
           plot = p, width = 8, height = 5, dpi = 150, bg = "white")
}

save_ppc_plot("ppc_prior_overlay.png",
              "Prior predictive: density of simulated eeg_mean")
save_ppc_plot("ppc_prior_max.png",
              "Prior predictive: max of simulated eeg_mean per draw",
              type = "stat", stat = "max")
save_ppc_plot("ppc_prior_min.png",
              "Prior predictive: min of simulated eeg_mean per draw",
              type = "stat", stat = "min")
save_ppc_plot("ppc_prior_sd.png",
              "Prior predictive: SD of simulated eeg_mean per draw",
              type = "stat", stat = "sd")
