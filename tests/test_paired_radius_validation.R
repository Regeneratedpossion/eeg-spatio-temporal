script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[[1]])
project_root <- normalizePath(file.path(dirname(script_path), ".."))
Sys.setenv(PROJECT_ROOT = project_root)

source(file.path(
    project_root, "code", "pilot", "00_validate_paired_radius_pilot.R"
))

sim <- expand.grid(
    participant = c("p1", "p2"),
    channel = c("CZ", "PZ"),
    condition = c("cond1", "cond2"),
    trial = 1:2,
    time = 0.25,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
)
sim$xproj <- ifelse(sim$channel == "CZ", 0, 1)
sim$yproj <- 0
sim$condition <- factor(sim$condition, levels = c("cond1", "cond2"))
sim$eeg <- ifelse(
    sim$condition == "cond2" & sim$channel == "CZ",
    2,
    0
)

truth <- list(
    signal_mode = "hard_signal",
    active_radius = 0.5,
    radius_rule = "1sd",
    active_electrodes = "CZ",
    continuous_effect = data.frame(
        channel = c("CZ", "PZ"),
        time = c(0.25, 0.25),
        spatial_distance = c(0, 1),
        spatial_weight_raw = c(1, exp(-2)),
        spatial_weight = c(1, 0),
        temporal_weight = c(1, 1),
        effect = c(2, 0)
    )
)

residual <- empirical_truth_residual(sim, truth)
stopifnot(all(residual$residual == 0))

summary <- validate_radius_objects(sim, truth, expected_active = 1L)
stopifnot(summary$n_active == 1L)
stopifnot(summary$max_abs_inactive_effect == 0)
stopifnot(summary$n_nonzero_inactive_effect == 0L)
stopifnot(summary$outermost_active_weight == 1)

machine_only <- compare_residual_vectors(
    c(0.1, -0.2),
    c(0.1 + 2e-16, -0.2 - 2e-16)
)
stopifnot(machine_only$equal_at_14_digits)
stopifnot(!machine_only$exactly_identical)

real_mismatch <- compare_residual_vectors(
    c(0.1, -0.2),
    c(0.1 + 1e-8, -0.2)
)
stopifnot(!real_mismatch$equal_at_14_digits)

bad_truth <- truth
bad_truth$continuous_effect$effect[2] <- 1e-8
stopifnot(inherits(
    try(validate_radius_objects(sim, bad_truth, 1L), silent = TRUE),
    "try-error"
))

cat("paired radius validation tests passed\n")
