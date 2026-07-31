test_project_root <- normalizePath(
    Sys.getenv("PROJECT_ROOT", unset = getwd()),
    mustWork = TRUE
)
source(file.path(test_project_root, "code", "utils.R"))

expect_error <- function(expr, pattern) {
    message_text <- tryCatch({
        force(expr)
        NA_character_
    }, error = function(e) conditionMessage(e))

    stopifnot(
        !is.na(message_text),
        grepl(pattern, message_text, fixed = TRUE)
    )
}

# Historical P20 truth reconstructed once from the actual 93-point ERP template,
# the 150--170 ms discrete onset distribution, and the 100 Hz analysis grid.
legacy_times <- seq(0.16, 0.34, by = 0.01)
legacy_truth <- c(
    0.07107146, 0.23111883, 0.43876655, 0.65501926, 0.86920727,
    1.06902382, 1.24163778, 1.37503274, 1.45938272, 1.48824700,
    1.45938272, 1.37503274, 1.24163778, 1.06902382, 0.86920727,
    0.65501926, 0.43876655, 0.23111883, 0.07107146
)

truth <- list(
    true_onset = 0.16,
    true_offset = 0.34,
    continuous_effect = data.frame(
        channel = "CZ",
        time = legacy_times,
        effect = legacy_truth
    )
)

selected <- continuous_truth_cells(truth, "CZ", c(0.16, 0.25, 0.34))
stopifnot(isTRUE(all.equal(
    selected$effect,
    legacy_truth[c(1, 10, 19)],
    tolerance = 1e-12
)))

truth_summary <- continuous_truth_summary(
    truth,
    channel = "CZ",
    peak_time = 0.25,
    window_times = legacy_times
)
stopifnot(isTRUE(all.equal(
    truth_summary$window_mean,
    0.8583562,
    tolerance = 1e-7
)))

# Regression values from the saved k=(16,12) and k=(12,12) posterior summaries.
k1612_recovery <- 100 * 0.796082710022911 / truth_summary$window_mean
k1212_recovery <- 100 * 0.726330090492104 / truth_summary$window_mean
stopifnot(abs(k1612_recovery - 92.7) < 0.1)
stopifnot(abs(k1212_recovery - 84.6) < 0.1)

expect_error(
    continuous_truth_cells(list(), "CZ", 0.25),
    "Exact continuous ground truth is unavailable"
)

malformed <- truth
malformed$continuous_effect$effect <- NULL
expect_error(
    continuous_truth_cells(malformed, "CZ", 0.25),
    "continuous_effect is missing required columns"
)

incomplete <- truth
incomplete$continuous_effect <- subset(
    incomplete$continuous_effect,
    time != 0.25
)
expect_error(
    continuous_truth_cells(incomplete, "CZ", 0.25),
    "Missing continuous ground-truth cells"
)

duplicated_truth <- truth
duplicated_truth$continuous_effect <- rbind(
    duplicated_truth$continuous_effect,
    duplicated_truth$continuous_effect[1, ]
)
expect_error(
    continuous_truth_cells(duplicated_truth, "CZ", 0.16),
    "Duplicated continuous ground-truth cells"
)

truth_matrix <- continuous_truth_matrix(truth, legacy_times, "CZ")
stopifnot(
    identical(dim(truth_matrix), c(length(legacy_times), 1L)),
    isTRUE(all.equal(as.vector(truth_matrix), legacy_truth, tolerance = 1e-12))
)

two_channel_truth <- list(
    continuous_effect = data.frame(
        channel = c("CZ", "FZ", "CZ", "FZ"),
        time = c(0.10, 0.10, 0.20, 0.20),
        effect = c(1, 10, 2, 20)
    )
)
two_channel_matrix <- continuous_truth_matrix(
    two_channel_truth,
    times = c(0.10, 0.20),
    channels = c("CZ", "FZ")
)
stopifnot(isTRUE(all.equal(
    two_channel_matrix,
    matrix(c(1, 2, 10, 20), nrow = 2,
           dimnames = list(c("0.1", "0.2"), c("CZ", "FZ")))
)))

generator_source <- paste(
    readLines(file.path(test_project_root, "code", "pipeline", "simulation.R")),
    collapse = "\n"
)
diagnostic_source <- paste(
    readLines(file.path(
        test_project_root, "code", "diagnostics", "02_diagnostics.R"
    )),
    collapse = "\n"
)
topomap_source <- paste(
    readLines(file.path(
        test_project_root, "code", "reporting", "08_slide_topomap_n20.R"
    )),
    collapse = "\n"
)
laplace_config_source <- paste(
    readLines(file.path(
        test_project_root, "code", "sensitivity", "laplacian", "00_config.R"
    )),
    collapse = "\n"
)
laplace_hard_source <- paste(
    readLines(file.path(
        test_project_root, "code", "sensitivity", "laplacian",
        "03_diagnose_hard.R"
    )),
    collapse = "\n"
)
laplace_continuous_source <- paste(
    readLines(file.path(
        test_project_root, "code", "sensitivity", "laplacian",
        "04_diagnose_continuous.R"
    )),
    collapse = "\n"
)
sesoi_plot_source <- paste(
    readLines(file.path(
        test_project_root, "code", "reporting",
        "31_plot_all_electrode_sesoi_posterior.R"
    )),
    collapse = "\n"
)

stopifnot(grepl("continuous_effect = continuous_ground_truth", generator_source,
                fixed = TRUE))
stopifnot(grepl("continuous_truth_summary", diagnostic_source, fixed = TRUE))
stopifnot(!grepl("truth_window_mean <- function", diagnostic_source, fixed = TRUE))
stopifnot(!grepl("sigma_t", diagnostic_source, fixed = TRUE))
stopifnot(!grepl("ground_truth$bump_amplitude", diagnostic_source, fixed = TRUE))
stopifnot(grepl("continuous_truth_cells", topomap_source, fixed = TRUE))
stopifnot(!grepl("temporal_weight <- function", topomap_source, fixed = TRUE))
stopifnot(!grepl("truth_temporal_bump <- function", laplace_config_source,
                 fixed = TRUE))
stopifnot(grepl("continuous_truth_summary", laplace_hard_source, fixed = TRUE))
stopifnot(!grepl("truth_temporal_bump", laplace_hard_source, fixed = TRUE))
stopifnot(!grepl("gt$bump_amplitude", laplace_hard_source, fixed = TRUE))
stopifnot(grepl("continuous_truth_matrix", laplace_continuous_source,
                fixed = TRUE))
stopifnot(!grepl("ground_truth$sigma_space", laplace_continuous_source,
                 fixed = TRUE))
stopifnot(!grepl("truth_temporal_bump", laplace_continuous_source, fixed = TRUE))
stopifnot(grepl("continuous_truth_cells", sesoi_plot_source, fixed = TRUE))

cat("continuous ground-truth regression tests passed\n")
