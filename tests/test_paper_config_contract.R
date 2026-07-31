script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[[1]])
project_root <- normalizePath(file.path(dirname(script_path), ".."))
Sys.setenv(PROJECT_ROOT = project_root)

source(file.path(project_root, "code", "config", "paper_config.R"))

cfg <- paper_analysis_config()
stopifnot(
    identical(cfg$simulation$scenario, "high"),
    identical(cfg$simulation$n_ppt, 20L),
    identical(cfg$simulation$n_trials, 20L),
    identical(cfg$simulation$n_reps, 20L),
    identical(cfg$simulation$signal_mode, "hard_signal"),
    identical(cfg$simulation$true_onset, 0.160),
    identical(cfg$simulation$true_offset, 0.342),
    identical(cfg$model$k_xy, 16L),
    identical(cfg$model$k_t, 12L),
    identical(cfg$model$bs_space, "ps"),
    identical(cfg$model$bs_time, "ps"),
    identical(cfg$model$sigma_true, TRUE),
    identical(cfg$radii$expected_active, c(5L, 9L, 21L))
)

cat("paper configuration contract tests passed\n")
