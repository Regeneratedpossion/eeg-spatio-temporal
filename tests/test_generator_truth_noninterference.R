script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[[1]])
project_root <- normalizePath(file.path(dirname(script_path), ".."))
generator_path <- file.path(project_root, "code", "pipeline", "simulation.R")
montecarlo_path <- file.path(
    project_root, "code", "pipeline", "generate_montecarlo.R"
)

generator_source <- paste(readLines(generator_path), collapse = "\n")
montecarlo_source <- paste(readLines(montecarlo_path), collapse = "\n")

# Saving diagnostic truth must retain an explicit historical replay mode.
stopifnot(grepl("legacy", generator_source, fixed = TRUE))
stopifnot(grepl("continuous_effect = continuous_ground_truth",
                generator_source, fixed = TRUE))
stopifnot(grepl("participant_onsets_ms", generator_source, fixed = TRUE))

# The Monte Carlo wrapper must make the scenario explicit.
stopifnot(grepl("MC_SCENARIO", montecarlo_source, fixed = TRUE))
stopifnot(grepl("scenario_manifest", montecarlo_source, fixed = TRUE))

cat("generator truth non-interference regression tests passed\n")
