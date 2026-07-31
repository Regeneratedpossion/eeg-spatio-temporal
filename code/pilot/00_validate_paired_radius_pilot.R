suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
})

source(file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()), "code/utils.R"))

paired_radius_validation_specs <- function() {
    data.frame(
        radius_rule = c("1sd", "half-max", "2sd"),
        data_subdir = c(
            "montecarlo_r1sd",
            "montecarlo_rhalfmax",
            "montecarlo_r2sd"
        ),
        expected_active = c(5L, 9L, 21L),
        stringsAsFactors = FALSE
    )
}

empirical_truth_residual <- function(sim_data, ground_truth) {
    empirical <- aggregate_trials(sim_data) %>%
        group_by(channel, time, condition) %>%
        summarise(value = mean(eeg_mean), .groups = "drop") %>%
        pivot_wider(names_from = condition, values_from = value) %>%
        mutate(empirical = cond2 - cond1) %>%
        select(channel, time, empirical)
    truth <- continuous_truth_cells(ground_truth)

    empirical %>%
        left_join(truth, by = c("channel", "time")) %>%
        mutate(residual = empirical - effect) %>%
        arrange(channel, time) %>%
        select(channel, time, residual)
}

compare_residual_vectors <- function(left, right) {
    stopifnot(length(left) == length(right))
    delta <- left - right
    list(
        exactly_identical = identical(left, right),
        equal_at_14_digits = identical(round(left, 14), round(right, 14)),
        max_abs_difference = max(abs(delta))
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
    residual_checks <- list()

    for (seed in seeds) {
        residuals <- list()
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
            residuals[[specs$radius_rule[i]]] <- empirical_truth_residual(
                sim_data, ground_truth
            )
        }

        pairs <- combn(specs$radius_rule, 2, simplify = FALSE)
        for (pair in pairs) {
            left <- residuals[[pair[1]]]
            right <- residuals[[pair[2]]]
            if (!identical(left[c("channel", "time")],
                           right[c("channel", "time")])) {
                stop("Residual grids differ for seed ", seed, call. = FALSE)
            }
            comparison <- compare_residual_vectors(
                left$residual, right$residual
            )
            residual_checks[[length(residual_checks) + 1L]] <- data.frame(
                seed = seed,
                left_rule = pair[1],
                right_rule = pair[2],
                exactly_identical = comparison$exactly_identical,
                equal_at_14_digits = comparison$equal_at_14_digits,
                max_abs_difference = comparison$max_abs_difference,
                stringsAsFactors = FALSE
            )
        }
    }

    summary_table <- bind_rows(summaries) %>%
        select(seed, everything())
    residual_table <- bind_rows(residual_checks)
    write.csv(
        summary_table,
        file.path(result_dir, "radius_validation.csv"),
        row.names = FALSE
    )
    write.csv(
        residual_table,
        file.path(result_dir, "paired_residual_validation.csv"),
        row.names = FALSE
    )
    print(summary_table)
    print(residual_table)
    if (!all(residual_table$equal_at_14_digits)) {
        stop("Empirical-minus-truth residuals differ beyond machine rounding",
             call. = FALSE)
    }
    invisible(list(summary = summary_table, residual = residual_table))
}

if (!interactive() && sys.nframe() == 0) {
    validate_paired_radius_pilot()
}
