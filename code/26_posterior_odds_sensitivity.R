###########################################################
# Posterior-odds sensitivity for heterogeneity BGAM results.
#
# Re-thresholds saved global-SESOI posterior cell probabilities. The fitted
# model, SESOI, continuous ground truth, cluster definition, and metric
# definitions are held fixed.
###########################################################

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(readr)
})

project_dir <- Sys.getenv("PROJECT_ROOT", unset = getwd())

.sensitivity_components <- function(mask, adjacency) {
    active_cells <- which(mask, arr.ind = TRUE)
    n_active <- nrow(active_cells)
    if (n_active == 0L) return(NULL)

    node_id <- matrix(0L, nrow(mask), ncol(mask))
    node_id[active_cells] <- seq_len(n_active)
    membership <- integer(n_active)
    component_id <- 0L

    for (start in seq_len(n_active)) {
        if (membership[start] != 0L) next
        component_id <- component_id + 1L
        membership[start] <- component_id
        queue <- start
        queue_position <- 1L

        while (queue_position <= length(queue)) {
            node <- queue[queue_position]
            queue_position <- queue_position + 1L
            time_index <- active_cells[node, 1]
            channel_index <- active_cells[node, 2]

            time_neighbours <- c(time_index - 1L, time_index + 1L)
            time_neighbours <- time_neighbours[
                time_neighbours >= 1L & time_neighbours <= nrow(mask)
            ]
            neighbour_nodes <- node_id[time_neighbours, channel_index]

            channel_neighbours <- which(adjacency[channel_index, ] == 1L)
            neighbour_nodes <- c(
                neighbour_nodes,
                node_id[time_index, channel_neighbours]
            )
            neighbour_nodes <- neighbour_nodes[neighbour_nodes > 0L]
            unseen <- neighbour_nodes[membership[neighbour_nodes] == 0L]
            if (length(unseen) > 0L) {
                unseen <- unique(unseen)
                membership[unseen] <- component_id
                queue <- c(queue, unseen)
            }
        }
    }

    list(
        active_cells = active_cells,
        membership = membership,
        size = tabulate(membership, nbins = component_id)
    )
}

.sensitivity_filter_components <- function(mask, adjacency, min_size) {
    components <- .sensitivity_components(mask, adjacency)
    if (is.null(components)) return(mask)
    keep <- components$size[components$membership] >= min_size
    output <- matrix(FALSE, nrow(mask), ncol(mask), dimnames = dimnames(mask))
    output[components$active_cells[keep, , drop = FALSE]] <- TRUE
    output
}

.sensitivity_timing <- function(mask, times, adjacency, min_size,
                                min_channels, min_run) {
    components <- .sensitivity_components(mask, adjacency)
    if (is.null(components)) {
        return(list(onset = NA_real_, offset = NA_real_))
    }
    eligible <- which(components$size >= min_size)
    if (length(eligible) == 0L) {
        return(list(onset = NA_real_, offset = NA_real_))
    }

    largest <- eligible[which.max(components$size[eligible])]
    time_index <- components$active_cells[
        components$membership == largest, 1
    ]
    n_detected <- tabulate(time_index, nbins = length(times))
    runs <- rle(n_detected >= min_channels)
    good <- which(runs$values & runs$lengths >= min_run)
    if (length(good) == 0L) {
        return(list(onset = NA_real_, offset = NA_real_))
    }

    starts <- cumsum(c(1L, head(runs$lengths, -1L)))
    ends <- cumsum(runs$lengths)
    list(
        onset = times[starts[good[1]]],
        offset = times[ends[tail(good, 1)]]
    )
}

.iou <- function(tp, fp, fn) tp / (tp + fp + fn)
.f1 <- function(tp, fp, fn) 2L * tp / (2L * tp + fp + fn)

.long_to_matrix <- function(data, value, channels, times) {
    value <- rlang::ensym(value)
    data %>%
        select(channel, time, !!value) %>%
        complete(time = times, channel = channels) %>%
        arrange(time, match(channel, channels)) %>%
        pull(!!value) %>%
        matrix(nrow = length(times), ncol = length(channels), byrow = TRUE,
               dimnames = list(as.character(times), channels))
}

compute_odds_sensitivity <- function(
        cells,
        truth_cells,
        adjacency,
        thresholds = c(1, 3, 6, 10, 20, 50, 100),
        min_cluster_size = 2L,
        min_onset_channels = 3L,
        min_onset_run = 2L) {
    channels <- sort(unique(as.character(cells$channel)))
    times <- sort(unique(cells$time))
    adjacency <- adjacency[channels, channels, drop = FALSE]

    odds <- .long_to_matrix(cells, odds, channels, times)
    practical_null <- .long_to_matrix(
        cells, practical_null, channels, times
    )
    truth_mask <- .long_to_matrix(
        truth_cells, true_active, channels, times
    )
    practical_null <- !is.na(practical_null) & practical_null
    truth_mask <- !is.na(truth_mask) & truth_mask

    truth_time <- rowSums(truth_mask) > 0L
    truth_onset <- min(times[truth_time])
    truth_offset <- max(times[truth_time])
    spatial_truth <- colSums(truth_mask) > 0L

    bind_rows(lapply(thresholds, function(threshold) {
        raw_pred <- !is.na(odds) & odds > threshold
        pred <- .sensitivity_filter_components(
            raw_pred, adjacency, min_size = min_cluster_size
        )
        timing <- .sensitivity_timing(
            pred,
            times,
            adjacency,
            min_size = min_cluster_size,
            min_channels = min_onset_channels,
            min_run = min_onset_run
        )

        st_tp <- sum(pred & truth_mask)
        st_fp <- sum(pred & !truth_mask)
        st_fn <- sum(!pred & truth_mask)
        spatial_pred <- colSums(pred) > 0L
        spatial_tp <- sum(spatial_pred & spatial_truth)
        spatial_fp <- sum(spatial_pred & !spatial_truth)
        spatial_fn <- sum(!spatial_pred & spatial_truth)

        tibble(
            odds_threshold = threshold,
            posterior_probability = threshold / (1 + threshold),
            truth_onset = truth_onset,
            estimated_onset = timing$onset,
            onset_error_ms = (timing$onset - truth_onset) * 1000,
            truth_offset = truth_offset,
            estimated_offset = timing$offset,
            offset_error_ms = (timing$offset - truth_offset) * 1000,
            spatial_tp = spatial_tp,
            spatial_fp = spatial_fp,
            spatial_fn = spatial_fn,
            spatial_iou = .iou(spatial_tp, spatial_fp, spatial_fn),
            spatial_F1 = .f1(spatial_tp, spatial_fp, spatial_fn),
            st_tp = st_tp,
            st_fp = st_fp,
            st_fn = st_fn,
            st_iou = .iou(st_tp, st_fp, st_fn),
            st_F1 = .f1(st_tp, st_fp, st_fn),
            n_pred = sum(pred),
            n_truth = sum(truth_mask),
            practical_null_fpr = mean(pred[practical_null])
        )
    }))
}

.global_sesoi <- function(project_root, scenario) {
    tag <- paste0(scenario, "_heterogeneity_seed001")
    read_csv(file.path(
        project_root, "results", "diagnostics",
        paste0(tag, "_sesoi"),
        paste0("delta_sesoi_thresholds_", tag, ".csv")
    ), show_col_types = FALSE) %>%
        filter(baseline_end == 0.1,
               threshold_type == "global_delta_sesoi") %>%
        pull(sesoi) %>%
        first()
}

run_heterogeneity_odds_sensitivity <- function(
        project_root = project_dir,
        scenarios = c("low", "medium", "high"),
        thresholds = c(1, 3, 6, 10, 20, 50, 100)) {
    adjacency <- readRDS(file.path(
        project_root, "resources", "adjacency",
        "adjacency_xproj_yproj_61ch.rds"
    ))

    results <- bind_rows(lapply(scenarios, function(scenario) {
        diagnostic_dir <- file.path(
            project_root, "results", "diagnostics",
            paste0(scenario, "_heterogeneity_seed001_sesoi")
        )
        tag <- paste0(scenario, "_heterogeneity_seed001")
        cells <- read_csv(file.path(
            diagnostic_dir,
            paste0("delta_sesoi_posterior_cells_", tag, ".csv")
        ), show_col_types = FALSE) %>%
            filter(baseline_end == 0.1,
                   threshold_type == "global_delta_sesoi")
        global_sesoi <- .global_sesoi(project_root, scenario)

        ground_truth <- readRDS(file.path(
            project_root, "data",
            paste0("montecarlo_p20_", scenario,
                   "_heterogeneity_exacttruth"),
            "gt_001.rds"
        ))
        truth_cells <- ground_truth$continuous_effect %>%
            transmute(channel, time, true_active = effect > global_sesoi)

        compute_odds_sensitivity(
            cells = cells,
            truth_cells = truth_cells,
            adjacency = adjacency,
            thresholds = thresholds
        ) %>%
            mutate(
                scenario = scenario,
                global_sesoi = global_sesoi,
                .before = 1
            )
    }))

    table_file <- file.path(
        project_root, "results", "tables",
        "heterogeneity_odds_sensitivity_exacttruth.csv"
    )
    write_csv(results, table_file)
    results
}

if (identical(Sys.getenv("RUN_THRESHOLD_SENSITIVITY"), "1")) {
    print(run_heterogeneity_odds_sensitivity())
}
