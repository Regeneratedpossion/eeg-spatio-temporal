suppressPackageStartupMessages({
    library(dplyr)
})

source(file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()), "code/utils.R"))

paired_radius_validation_specs <- function() {
    data.frame(
        radius_rule = "half-max",
        data_subdir = "montecarlo_rhalfmax",
        expected_active = 9L,
        stringsAsFactors = FALSE
    )
}

validate_radius_objects <- function(sim_data, ground_truth, expected_active) {
    if (!identical(ground_truth$signal_mode, "hard_signal")) {
        stop("Ground truth is not hard_signal", call. = FALSE)
    }
    if (is.null(ground_truth$active_radius) ||
        is.null(ground_truth$radius_rule)) {
        stop("Missing radius metadata", call. = FALSE)
    }
    n_active <- length(ground_truth$active_electrodes)
    if (n_active != expected_active) {
        stop("Expected ", expected_active, " active electrodes, found ",
             n_active, call. = FALSE)
    }

    ce <- ground_truth$continuous_effect
    inactive <- !(ce$channel %in% ground_truth$active_electrodes)
    max_abs_inactive <- max(abs(ce$effect[inactive]))
    n_nonzero_inactive <- sum(ce$effect[inactive] != 0)
    if (max_abs_inactive != 0 || n_nonzero_inactive != 0L) {
        stop("Inactive ground-truth effects are not exactly zero", call. = FALSE)
    }

    active_spatial <- ce %>%
        filter(channel %in% ground_truth$active_electrodes) %>%
        distinct(channel, spatial_distance, spatial_weight_raw) %>%
        arrange(desc(spatial_distance))
    outer_distance <- active_spatial$spatial_distance[1]
    outer <- active_spatial %>%
        filter(spatial_distance == outer_distance)

    data.frame(
        radius_rule = ground_truth$radius_rule,
        active_radius = ground_truth$active_radius,
        n_active = n_active,
        max_abs_inactive_effect = max_abs_inactive,
        n_nonzero_inactive_effect = n_nonzero_inactive,
        outermost_active_electrode = paste(sort(outer$channel), collapse = ","),
        outermost_active_distance = outer_distance,
        outermost_active_weight = outer$spatial_weight_raw[1],
        stringsAsFactors = FALSE
    )
}

validate_paired_radius_pilot <- function(seeds = 1:5) {
    specs <- paired_radius_validation_specs()
    result_dir <- file.path(project_root(), "results/paired_radius_pilot")
    dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)
    summaries <- list()

    for (seed in seeds) {
        for (i in seq_len(nrow(specs))) {
            sim_path <- data_path(
                specs$data_subdir[i], sprintf("sim_%03d.rds", seed)
            )
            gt_path <- data_path(
                specs$data_subdir[i], sprintf("gt_%03d.rds", seed)
            )
            sim_data <- readRDS(sim_path)
            ground_truth <- readRDS(gt_path)
            summary <- validate_radius_objects(
                sim_data, ground_truth, specs$expected_active[i]
            )
            summary$seed <- seed
            summaries[[length(summaries) + 1L]] <- summary
        }
    }

    summary_table <- bind_rows(summaries) %>%
        select(seed, everything())
    write.csv(
        summary_table,
        file.path(result_dir, "radius_validation.csv"),
        row.names = FALSE
    )
    print(summary_table)
    invisible(summary_table)
}

if (!interactive() && sys.nframe() == 0) {
    validate_paired_radius_pilot()
}
