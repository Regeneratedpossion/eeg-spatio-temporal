script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[[1]])
test_project_root <- normalizePath(file.path(dirname(script_path), ".."))
Sys.setenv(PROJECT_ROOT = test_project_root)

source(file.path(test_project_root, "code", "pipeline", "fit_paired_radius.R"))
source(file.path(
    test_project_root, "code", "diagnostics",
    "35_diagnose_paired_radius_bgam_pilot.R"
))

config <- paired_radius_bgam_config()
fit_arguments <- names(formals(fit_bgam))
paired_fit_arguments <- names(formals(run_paired_radius_bgam_fit))
stopifnot(all(c("chains", "cores", "iter", "warmup", "threads") %in%
              fit_arguments))
stopifnot(identical(paired_fit_arguments, c("radius_rule", "seed")))
stopifnot(config$k_xy == 16L)
stopifnot(config$k_t == 12L)
stopifnot(identical(config$bs_space, "ps"))
stopifnot(identical(config$bs_time, "ps"))
stopifnot(config$chains == 4L)
stopifnot(config$iter == 3000L)
stopifnot(config$warmup == 1000L)
stopifnot(config$adapt_delta == 0.80)
stopifnot(config$max_treedepth == 10L)
stopifnot(identical(config$se_variable, "eeg_sem"))
stopifnot(isTRUE(config$sigma_true))
stopifnot(identical(config$prior_mode, "brms_default"))

specs <- paired_radius_bgam_specs()
stopifnot(identical(specs$radius_rule, c("1sd", "half-max", "2sd")))
stopifnot(identical(
    specs$data_subdir,
    c("montecarlo_r1sd", "montecarlo_rhalfmax", "montecarlo_r2sd")
))

theta <- c(2, 1, 0)
theta_hat <- c(1, 1, 0.5)
stopifnot(projection_recovery(theta_hat, theta) == 0.6)
stopifnot(isTRUE(all.equal(
    spatial_leakage(theta_hat, c(TRUE, TRUE, FALSE)),
    0.25 / 2.25,
    tolerance = 0
)))

metrics <- spatial_fp_metrics(
    predicted = c("B", "C"),
    active = c("A", "B"),
    all_channels = c("A", "B", "C", "D", "E")
)
stopifnot(metrics$TP == 1L, metrics$FP == 1L, metrics$FN == 1L)
stopifnot(metrics$n_inactive == 3L, metrics$FPR == 1 / 3)

sensors <- data.frame(
    channel = c("A", "B", "C", "D"),
    xproj = c(0, 1, 2, 0),
    yproj = c(0, 0, 0, 2)
)
distances <- fp_nearest_active_distance(
    fp_channels = c("C", "D"),
    active_channels = c("A", "B"),
    sensors = sensors
)
stopifnot(distances$nearest_active_distance[distances$channel == "C"] == 1)
stopifnot(distances$nearest_active_distance[distances$channel == "D"] == 2)

posterior_fixture <- expand.grid(
    channel = c("A", "B"),
    time = c(0, 0.01),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
)
posterior_fixture$p_sig <- 0.99
adjacency_fixture <- matrix(
    c(0, 1, 1, 0),
    nrow = 2,
    dimnames = list(c("A", "B"), c("A", "B"))
)
bgam_only_masks <- paired_radius_method_masks(
    posterior_map = posterior_fixture,
    channels = c("A", "B"),
    times = c(0, 0.01),
    adjacency = adjacency_fixture,
    bgam_only = TRUE
)
stopifnot(identical(names(bgam_only_masks), "BGAM_2D"))
stopifnot(all(bgam_only_masks$BGAM_2D))

cat("paired-radius BGAM pilot tests passed\n")
