###########################################################
# Significance-threshold sensitivity on heterogeneity data.
#
# Computes each frequentist p-value map once, then applies alpha values
# 0.05, 0.01, 0.005, and 0.001. Evaluation uses the same exact continuous
# truth and spatio-temporal cluster rules as the BGAM odds sensitivity.
###########################################################

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(readr)
})

project_dir <- Sys.getenv("PROJECT_ROOT", unset = getwd())
source(file.path(
    project_dir, "code", "21_posterior_odds_sensitivity.R"
))

compute_alpha_sensitivity <- function(
        pvalue_cells,
        truth_cells,
        adjacency,
        alphas = c(0.05, 0.01, 0.005, 0.001),
        min_cluster_size = 2L,
        min_onset_channels = 3L,
        min_onset_run = 2L) {
    bind_rows(lapply(unique(pvalue_cells$method), function(method_name) {
        method_cells <- pvalue_cells %>%
            filter(method == method_name) %>%
            mutate(odds = (1 - p_value) / p_value)
        equivalent_odds <- (1 - alphas) / alphas
        compute_odds_sensitivity(
            cells = method_cells,
            truth_cells = truth_cells,
            adjacency = adjacency,
            thresholds = equivalent_odds,
            min_cluster_size = min_cluster_size,
            min_onset_channels = min_onset_channels,
            min_onset_run = min_onset_run
        ) %>%
            mutate(method = method_name, alpha = alphas, .before = 1) %>%
            select(-odds_threshold, -posterior_probability)
    }))
}

.participant_difference_array <- function(sim_data) {
    cell_means <- sim_data %>%
        group_by(participant, channel, xproj, yproj, condition, time) %>%
        summarise(eeg_mean = mean(eeg), .groups = "drop") %>%
        select(participant, channel, time, condition, eeg_mean) %>%
        pivot_wider(names_from = condition, values_from = eeg_mean) %>%
        mutate(diff = cond2 - cond1)

    participants <- sort(unique(cell_means$participant))
    times <- sort(unique(cell_means$time))
    channels <- sort(unique(as.character(cell_means$channel)))
    cell_means <- cell_means %>%
        mutate(
            participant_index = match(participant, participants),
            time_index = match(time, times),
            channel_index = match(channel, channels)
        )
    differences <- array(
        NA_real_,
        dim = c(length(participants), length(times), length(channels)),
        dimnames = list(participants, as.character(times), channels)
    )
    differences[cbind(
        cell_means$participant_index,
        cell_means$time_index,
        cell_means$channel_index
    )] <- cell_means$diff

    list(array = differences, times = times, channels = channels)
}

.cell_t_test <- function(values) stats::t.test(values)$p.value

.mne_pvalue_map <- function(differences, adjacency, threshold = NULL,
                            permutations = 1024L, seed = 42L) {
    reticulate <- asNamespace("reticulate")
    mne <- reticulate$import("mne")
    scipy_sparse <- reticulate$import("scipy.sparse")
    adjacency_sparse <- scipy_sparse$csr_matrix(
        matrix(as.double(adjacency), nrow = nrow(adjacency))
    )
    arguments <- list(
        differences,
        adjacency = adjacency_sparse,
        n_permutations = as.integer(permutations),
        n_jobs = 1L,
        seed = as.integer(seed),
        out_type = "mask",
        threshold = threshold
    )
    result <- do.call(
        mne$stats$spatio_temporal_cluster_1samp_test,
        arguments
    )
    clusters <- result[[2]]
    cluster_pvalues <- result[[3]]
    pvalues <- matrix(1, dim(differences)[2], dim(differences)[3])
    for (cluster_index in seq_along(clusters)) {
        indices <- which(as.matrix(clusters[[cluster_index]]), arr.ind = TRUE)
        pvalues[indices] <- pmin(
            pvalues[indices], cluster_pvalues[cluster_index]
        )
    }
    pvalues
}

.matrix_to_pvalue_cells <- function(pvalues, method, channels, times,
                                    practical_null) {
    expand_grid(time = times, channel = channels) %>%
        mutate(
            method = method,
            p_value = as.vector(t(pvalues)),
            practical_null = practical_null[match(channel, channels)]
        ) %>%
        select(method, channel, time, p_value, practical_null)
}

run_heterogeneity_alpha_sensitivity <- function(
        project_root = project_dir,
        scenarios = c("low", "high"),
        alphas = c(0.05, 0.01, 0.005, 0.001),
        permutations = 1024L) {
    adjacency <- readRDS(file.path(
        project_root, "resources", "adjacency",
        "adjacency_xproj_yproj_61ch.rds"
    ))
    results <- bind_rows(lapply(scenarios, function(scenario) {
        message("Computing p-value maps for ", scenario, " heterogeneity...")
        sim_data <- readRDS(file.path(
            project_root, "data",
            paste0("montecarlo_p20_", scenario, "_heterogeneity"),
            "sim_001.rds"
        ))
        ground_truth <- readRDS(file.path(
            project_root, "data",
            paste0("montecarlo_p20_", scenario,
                   "_heterogeneity_exacttruth"),
            "gt_001.rds"
        ))
        global_sesoi <- .global_sesoi(project_root, scenario)

        bundle <- .participant_difference_array(sim_data)
        channels <- bundle$channels
        times <- bundle$times
        adjacency_scenario <- adjacency[channels, channels, drop = FALSE]
        raw_p <- apply(bundle$array, c(2, 3), .cell_t_test)
        bh_p <- matrix(
            p.adjust(raw_p, method = "BH"), nrow(raw_p), ncol(raw_p)
        )
        holm_p <- matrix(
            p.adjust(raw_p, method = "holm"), nrow(raw_p), ncol(raw_p)
        )
        cluster_p <- .mne_pvalue_map(
            bundle$array, adjacency_scenario,
            permutations = permutations
        )
        tfce_p <- .mne_pvalue_map(
            bundle$array, adjacency_scenario,
            threshold = reticulate::dict(start = 0, step = 0.2),
            permutations = permutations
        )

        sensors <- sim_data %>%
            select(channel, xproj, yproj) %>%
            distinct() %>%
            arrange(match(channel, channels))
        seed_sensor <- sensors %>%
            filter(channel == ground_truth$seed_electrode)
        practical_null <- sqrt(
            (sensors$xproj - seed_sensor$xproj)^2 +
                (sensors$yproj - seed_sensor$yproj)^2
        ) > 2 * ground_truth$sigma_space

        pvalue_cells <- bind_rows(
            .matrix_to_pvalue_cells(
                raw_p, "Raw p", channels, times, practical_null
            ),
            .matrix_to_pvalue_cells(
                bh_p, "FDR (BH)", channels, times, practical_null
            ),
            .matrix_to_pvalue_cells(
                holm_p, "Holm", channels, times, practical_null
            ),
            .matrix_to_pvalue_cells(
                cluster_p, "Cluster", channels, times, practical_null
            ),
            .matrix_to_pvalue_cells(
                tfce_p, "TFCE", channels, times, practical_null
            )
        )
        truth_cells <- paired_radius_superpopulation_truth(
            ground_truth,
            file.path(project_dir, "resources", "brms_meeg", "code", "erp_template.R")
        ) %>%
            transmute(channel, time, true_active = effect > global_sesoi)

        compute_alpha_sensitivity(
            pvalue_cells = pvalue_cells,
            truth_cells = truth_cells,
            adjacency = adjacency_scenario,
            alphas = alphas
        ) %>%
            mutate(
                scenario = scenario,
                global_sesoi = global_sesoi,
                permutations = permutations,
                .before = 1
            )
    }))

    table_file <- file.path(
        project_root, "results", "tables",
        "heterogeneity_alpha_sensitivity_exacttruth.csv"
    )
    write_csv(results, table_file)
    results
}

if (identical(Sys.getenv("RUN_ALPHA_SENSITIVITY"), "1")) {
    print(run_heterogeneity_alpha_sensitivity())
}
