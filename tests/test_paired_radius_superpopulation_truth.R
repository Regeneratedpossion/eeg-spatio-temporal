test_project_root <- Sys.getenv("PROJECT_ROOT", unset = getwd())

source(file.path(
    test_project_root, "code", "paired_radius_superpopulation_truth.R"
))

template_path <- file.path(
    test_project_root, "resources", "brms_meeg", "code", "erp_template.R"
)
ground_truth <- readRDS(file.path(
    test_project_root, "data", "montecarlo_r2sd", "gt_001.rds"
))

population <- paired_radius_superpopulation_truth(
    ground_truth, template_path = template_path
)
cz_250 <- subset(
    population,
    channel == "CZ" & abs(time - 0.25) < 1e-12
)$effect

stopifnot(length(cz_250) == 1L)
stopifnot(isTRUE(all.equal(cz_250, 1.457449, tolerance = 1e-6)))
stopifnot(nrow(population) == 61L * 51L)
stopifnot(all(subset(
    population,
    !channel %in% ground_truth$active_electrodes
)$effect == 0))

changed_realisation <- ground_truth
changed_realisation$participant_onsets_ms[] <- 180
changed_realisation$participant_parameters$onset_ms[] <- 180
changed_realisation$participant_parameters$amplitude[] <- 100
changed_realisation$continuous_effect$effect[] <- 999
changed_realisation$population_mean_waveform$temporal_weight[] <- 999
population_after_change <- paired_radius_superpopulation_truth(
    changed_realisation, template_path = template_path
)
stopifnot(identical(population, population_after_change))

selected <- paired_radius_superpopulation_truth(
    ground_truth,
    template_path = template_path,
    channels = c("CZ", "AF1"),
    times = 0.25
)
stopifnot(identical(selected$channel, c("AF1", "CZ")))
stopifnot(all(selected$time == 0.25))

expect_error <- function(expression, pattern) {
    message <- tryCatch(
        {
            force(expression)
            NA_character_
        },
        error = function(error) conditionMessage(error)
    )
    stopifnot(!is.na(message), grepl(pattern, message, fixed = TRUE))
}

missing_support <- ground_truth
missing_support$onset_support_ms <- NULL
expect_error(
    paired_radius_superpopulation_truth(
        missing_support, template_path = template_path
    ),
    "ground_truth is missing required fields: onset_support_ms"
)

diagnostic_source <- paste(
    readLines(file.path(
        test_project_root,
        "code",
        "35_diagnose_paired_radius_bgam_pilot.R"
    )),
    collapse = "\n"
)
stopifnot(grepl(
    'source(file.path(project_root(), "code/paired_radius_superpopulation_truth.R"))',
    diagnostic_source,
    fixed = TRUE
))
stopifnot(grepl(
    "paired_radius_superpopulation_truth(",
    diagnostic_source,
    fixed = TRUE
))
stopifnot(grepl(
    'truth_estimand = "superpopulation"',
    diagnostic_source,
    fixed = TRUE
))
stopifnot(grepl(
    "realized_truth_cz_250 = realized_truth_peak",
    diagnostic_source,
    fixed = TRUE
))

cat("paired-radius superpopulation truth tests passed\n")
