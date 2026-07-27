script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[[1]])
project_root <- normalizePath(file.path(dirname(script_path), ".."))
generator_path <- file.path(project_root, "code", "00_generate_data.R")
montecarlo_path <- file.path(project_root, "code", "00_generate_montecarlo.R")

generator_source <- paste(readLines(generator_path), collapse = "\n")
montecarlo_source <- paste(readLines(montecarlo_path), collapse = "\n")

# Saving diagnostic truth must not alter the historical signal or noise model.
forbidden_generator_patterns <- c(
    "heterogeneity_scenarios",
    "draw_participant_parameters",
    "participant_amplitude",
    "participant\\$noise_scale",
    "omega_A",
    "omega_sigma"
)
stopifnot(!any(vapply(
    forbidden_generator_patterns,
    grepl,
    logical(1),
    x = generator_source
)))

stopifnot(grepl("continuous_effect = continuous_ground_truth",
                generator_source, fixed = TRUE))
stopifnot(grepl("participant_onsets_ms", generator_source, fixed = TRUE))

# The Monte Carlo wrapper must retain its historical inputs and not introduce
# a new simulation scenario while adding the saved truth field.
stopifnot(!grepl("MC_SCENARIO", montecarlo_source, fixed = TRUE))
stopifnot(!grepl("scenario_manifest", montecarlo_source, fixed = TRUE))

cat("generator truth non-interference regression tests passed\n")
