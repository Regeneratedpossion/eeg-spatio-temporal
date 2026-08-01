suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(mgcv)
})

source(file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()), "code/utils.R"))
comparison_dependency <- file.path(
    project_root(), "code", "paired_radius_compare_methods_dependency.R"
)
if (!file.exists(comparison_dependency)) {
    comparison_dependency <- file.path(project_root(), "code/04_compare_methods.R")
}
source(comparison_dependency)
source(file.path(project_root(), "code/08_fit_paired_radius_bgam.R"))
source(file.path(project_root(), "code/ground_truth.R"))

projection_recovery <- function(theta_hat, theta) {
    sum(theta_hat * theta) / sum(theta^2)
}

spatial_leakage <- function(theta_hat, active) {
    sum(theta_hat[!active]^2) / sum(theta_hat^2)
}

spatial_fp_metrics <- function(predicted, active, all_channels) {
    predicted <- intersect(predicted, all_channels)
    active <- intersect(active, all_channels)
    tp <- length(intersect(predicted, active))
    fp <- length(setdiff(predicted, active))
    fn <- length(setdiff(active, predicted))
    n_inactive <- length(setdiff(all_channels, active))
    precision <- if ((tp + fp) > 0L) tp / (tp + fp) else NA_real_
    recall <- if ((tp + fn) > 0L) tp / (tp + fn) else NA_real_
    f1 <- if (is.na(precision) || is.na(recall) ||
              (precision + recall) == 0) {
        NA_real_
    } else {
        2 * precision * recall / (precision + recall)
    }
    list(
        TP = tp, FP = fp, FN = fn, n_inactive = n_inactive,
        FPR = fp / n_inactive, precision = precision,
        recall = recall, F1 = f1
    )
}

fp_nearest_active_distance <- function(fp_channels, active_channels, sensors) {
    if (!length(fp_channels)) {
        return(data.frame(
            channel = character(), nearest_active_channel = character(),
            nearest_active_distance = numeric()
        ))
    }
    active <- sensors %>% filter(channel %in% active_channels)
    bind_rows(lapply(fp_channels, function(fp_channel) {
        fp <- sensors %>% filter(channel == fp_channel)
        distances <- sqrt(
            (active$xproj - fp$xproj)^2 +
            (active$yproj - fp$yproj)^2
        )
        nearest <- which.min(distances)
        data.frame(
            channel = fp_channel,
            nearest_active_channel = active$channel[nearest],
            nearest_active_distance = distances[nearest]
        )
    }))
}

fit_edf_companion <- function(sim_model, config) {
    edf_data <- sim_model %>%
        mutate(
            participant = factor(participant),
            condition = factor(condition, levels = c("cond1", "cond2")),
            condition_numeric = if_else(condition == "cond2", 0.5, -0.5)
        )
    contrasts(edf_data$condition) <- c(-0.5, 0.5)
    edf_fit <- bam(
        eeg_mean ~ condition +
            s(participant, bs = "re") +
            s(participant, by = condition_numeric, bs = "re") +
            s(xproj, yproj, by = condition, bs = "tp", k = config$k_xy) +
            s(time, by = condition, bs = config$bs_time, k = config$k_t),
        data = edf_data,
        weights = 1 / pmax(edf_data$eeg_sem^2, 1e-8),
        method = "fREML",
        discrete = FALSE,
        nthreads = 4
    )
    smooth_table <- as.data.frame(summary(edf_fit)$s.table)
    smooth_table$smooth <- rownames(smooth_table)
    marginal_smooths <- Filter(
        function(x) startsWith(x$label, "s(xproj,yproj)") ||
            startsWith(x$label, "s(time)"),
        edf_fit$smooth
    )
    components <- bind_rows(lapply(marginal_smooths, function(smooth) {
        row <- smooth_table[smooth_table$smooth == smooth$label, , drop = FALSE]
        stopifnot(nrow(row) == 1L)
        data.frame(
            smooth = smooth$label,
            condition = sub("^.*:condition", "", smooth$label),
            edf = row$edf,
            upper = smooth$bs.dim - 1,
            direction = if (startsWith(smooth$label, "s(xproj,yproj)")) {
                "space"
            } else {
                "time"
            },
            stringsAsFactors = FALSE
        )
    })) %>%
        mutate(
            fraction_of_upper = edf / upper
        )
    summary <- components %>%
        group_by(direction) %>%
        slice_max(fraction_of_upper, n = 1, with_ties = FALSE) %>%
        transmute(
            direction,
            edf,
            upper,
            fraction_of_upper,
            near_upper = fraction_of_upper >= 0.90
        ) %>%
        ungroup()
    list(fit = edf_fit, components = components, summary = summary)
}

paired_radius_posterior_map <- function(fit, sim_model) {
    sensors <- sim_model %>%
        distinct(channel, xproj, yproj) %>%
        arrange(channel)
    channels <- sensors$channel
    times <- sort(unique(sim_model$time))
    grid <- expand_grid(channel = channels, time = times) %>%
        left_join(sensors, by = "channel") %>%
        mutate(eeg_sem = 1)
    c1 <- grid %>%
        mutate(condition = factor("cond1", levels = c("cond1", "cond2")))
    c2 <- grid %>%
        mutate(condition = factor("cond2", levels = c("cond1", "cond2")))
    draws <- posterior_epred(fit, newdata = c2, re_formula = NA) -
        posterior_epred(fit, newdata = c1, re_formula = NA)
    grid %>%
        transmute(
            channel, time, xproj, yproj,
            mean_diff = colMeans(draws),
            p_sig = colMeans(draws > ROPE_UV),
            odds = pmax(p_sig, 1 / nrow(draws)) /
                pmax(1 - p_sig, 1 / nrow(draws))
        )
}

paired_radius_method_masks <- function(
        posterior_map, channels, times, adjacency, bgam_only = FALSE,
        arr = NULL, mne = NULL) {
    if (isTRUE(bgam_only)) {
        odds_matrix <- .bgam_full_to_odds(posterior_map, channels, times)
        return(list(
            BGAM_2D = .filter_components(
                odds_mask(odds_matrix),
                adjacency,
                min_size = MIN_CLUSTER_SIZE
            )
        ))
    }
    if (is.null(arr)) {
        stop("arr is required unless bgam_only is TRUE", call. = FALSE)
    }
    build_method_masks(
        arr, channels, times, adjacency, posterior_map,
        mne = mne, alpha = ALPHA
    )
}

plot_empirical_minus_bgam <- function(sim_model, posterior_map, radius_rule,
                                      out_file) {
    target_time <- pick_time(0.25, sort(unique(sim_model$time)))
    empirical <- sim_model %>%
        filter(time == target_time) %>%
        group_by(channel, xproj, yproj, condition) %>%
        summarise(value = mean(eeg_mean), .groups = "drop") %>%
        pivot_wider(names_from = condition, values_from = value) %>%
        mutate(empirical = cond2 - cond1) %>%
        select(channel, xproj, yproj, empirical)
    residual <- posterior_map %>%
        filter(time == target_time) %>%
        select(channel, bgam = mean_diff) %>%
        left_join(empirical, by = "channel") %>%
        mutate(value = empirical - bgam)
    head_radius <- max(sqrt(residual$xproj^2 + residual$yproj^2)) * 1.05
    grid <- expand.grid(
        xproj = seq(-head_radius, head_radius, length.out = 140),
        yproj = seq(-head_radius, head_radius, length.out = 140)
    ) %>% filter(sqrt(xproj^2 + yproj^2) <= head_radius)
    grid$value <- idw_interp(
        grid, residual %>% select(xproj, yproj), residual$value,
        power = 2, k = 6
    )
    limit <- max(abs(grid$value))
    decor <- scalp_decorations(head_radius)
    plot <- ggplot(grid, aes(xproj, yproj)) +
        geom_raster(aes(fill = value), interpolate = TRUE) +
        geom_contour(aes(z = value), bins = 7, color = "grey25",
                     linewidth = 0.2) +
        geom_path(data = decor$head_circle, aes(x, y), inherit.aes = FALSE) +
        geom_path(data = decor$nose, aes(x, y), inherit.aes = FALSE) +
        geom_path(data = decor$ear_left, aes(x, y), inherit.aes = FALSE) +
        geom_path(data = decor$ear_right, aes(x, y), inherit.aes = FALSE) +
        geom_point(data = residual, aes(xproj, yproj), inherit.aes = FALSE,
                   size = 0.8) +
        coord_equal() +
        scale_fill_gradient2(
            low = "#2166ac", mid = "white", high = "#b2182b",
            midpoint = 0, limits = c(-limit, limit), name = "uV"
        ) +
        labs(title = paste("Empirical - BGAM:", radius_rule),
             subtitle = "250 ms", x = NULL, y = NULL) +
        theme_void(base_size = 11)
    ggsave(out_file, plot, width = 5, height = 5, dpi = 300)
    invisible(residual)
}

run_paired_radius_bgam_diagnostics <- function(radius_rule,
                                               bgam_only = FALSE) {
    config <- paired_radius_bgam_config()
    specs <- paired_radius_bgam_specs()
    spec <- specs[specs$radius_rule == radius_rule, , drop = FALSE]
    if (nrow(spec) != 1L) {
        stop("Unknown radius_rule: ", radius_rule, call. = FALSE)
    }
    out_dir <- paired_radius_output_dir(radius_rule)
    fit <- readRDS(file.path(out_dir, "fit_seed001.rds"))
    sim_data <- readRDS(data_path(spec$data_subdir, "sim_001.rds"))
    ground_truth <- readRDS(data_path(spec$data_subdir, "gt_001.rds"))
    stopifnot(identical(readRDS(file.path(out_dir, "model_config.rds")),
                        config))

    sim_model <- aggregate_trials(sim_data)
    sensors <- sim_model %>%
        distinct(channel, xproj, yproj) %>%
        arrange(channel)
    channels <- sensors$channel
    times <- sort(unique(sim_model$time))
    target_time <- pick_time(0.25, times)

    posterior_map <- paired_radius_posterior_map(fit, sim_model)
    saveRDS(posterior_map, file.path(out_dir, "posterior_map_seed001.rds"))

    adjacency <- readRDS(adjacency_path("adjacency_xproj_yproj_61ch.rds"))
    adjacency <- adjacency[channels, channels]
    bundle <- if (isTRUE(bgam_only)) NULL else build_diff_array(sim_model)
    method_masks <- paired_radius_method_masks(
        posterior_map = posterior_map,
        channels = channels,
        times = times,
        adjacency = adjacency,
        bgam_only = bgam_only,
        arr = if (is.null(bundle)) NULL else bundle$arr,
        mne = mne
    )
    time_index <- match(target_time, times)
    truth_at_peak <- paired_radius_superpopulation_truth(
        ground_truth,
        template_path = file.path(
            project_root(), "resources", "brms_meeg", "code",
            "erp_template.R"
        ),
        channels = channels,
        times = target_time
    ) %>% arrange(match(channel, channels))
    realized_truth_at_peak <- continuous_truth_cells(
        ground_truth, channels = channels, times = target_time
    ) %>% arrange(match(channel, channels))
    active_channels <- ground_truth$active_electrodes

    selected_methods <- if (isTRUE(bgam_only)) {
        c(BGAM = "BGAM_2D")
    } else {
        c(BGAM = "BGAM_2D", cluster_mass = "cluster", TFCE = "tfce")
    }
    method_rows <- list()
    fp_rows <- list()
    for (method in names(selected_methods)) {
        mask <- method_masks[[selected_methods[[method]]]][time_index, ]
        predicted <- channels[which(mask)]
        metrics <- spatial_fp_metrics(predicted, active_channels, channels)
        method_rows[[method]] <- data.frame(
            radius_rule = radius_rule,
            method = method,
            target_time = target_time,
            TP = metrics$TP, FP = metrics$FP, FN = metrics$FN,
            n_inactive = metrics$n_inactive, FPR = metrics$FPR,
            precision = metrics$precision, recall = metrics$recall,
            F1 = metrics$F1,
            stringsAsFactors = FALSE
        )
        fp_channels <- setdiff(predicted, active_channels)
        if (length(fp_channels)) {
            fp_detail <- fp_nearest_active_distance(
                fp_channels, active_channels, sensors
            ) %>%
                left_join(
                    truth_at_peak %>%
                        select(channel, true_effect = effect),
                    by = "channel"
                ) %>%
                mutate(radius_rule = radius_rule, method = method,
                       .before = 1)
            fp_rows[[method]] <- fp_detail
        }
    }
    method_metrics <- bind_rows(method_rows)
    fp_details <- bind_rows(fp_rows)
    write.csv(method_metrics, file.path(out_dir, "method_metrics_250ms.csv"),
              row.names = FALSE)
    write.csv(fp_details, file.path(out_dir, "fp_details_250ms.csv"),
              row.names = FALSE)

    bgam_peak <- posterior_map %>%
        filter(channel == "CZ", time == target_time) %>%
        pull(mean_diff)
    truth_peak <- truth_at_peak %>%
        filter(channel == "CZ") %>%
        pull(effect)
    realized_truth_peak <- realized_truth_at_peak %>%
        filter(channel == "CZ") %>%
        pull(effect)
    theta_hat <- posterior_map %>%
        filter(time == target_time) %>%
        arrange(match(channel, channels)) %>%
        pull(mean_diff)
    theta <- truth_at_peak$effect
    active <- channels %in% active_channels

    np <- nuts_params(fit)
    draws_summary <- posterior::summarise_draws(as_draws_df(fit))
    n_divergent <- sum(np$Value[np$Parameter == "divergent__"])
    n_treedepth <- sum(
        np$Value[np$Parameter == "treedepth__"] >= config$max_treedepth
    )

    edf <- fit_edf_companion(sim_model, config)
    write.csv(edf$components, file.path(out_dir, "edf_components.csv"),
              row.names = FALSE)
    write.csv(edf$summary, file.path(out_dir, "edf_summary.csv"),
              row.names = FALSE)
    direction_value <- function(direction, column) {
        value <- edf$summary[edf$summary$direction == direction, column]
        if (length(value)) value[[1]] else NA
    }
    bgam_metrics <- method_metrics %>% filter(method == "BGAM")
    summary <- data.frame(
        truth_estimand = "superpopulation",
        diagnostic_scope = if (isTRUE(bgam_only)) "BGAM_only" else "all_methods",
        radius_rule = radius_rule,
        jump_weight = spec$jump_weight,
        n_active = length(active_channels),
        n_inactive = sum(!active),
        spatial_edf = direction_value("space", "edf"),
        spatial_edf_upper = direction_value("space", "upper"),
        spatial_edf_fraction = direction_value("space", "fraction_of_upper"),
        spatial_edf_near_upper = direction_value("space", "near_upper"),
        temporal_edf = direction_value("time", "edf"),
        temporal_edf_upper = direction_value("time", "upper"),
        temporal_edf_fraction = direction_value("time", "fraction_of_upper"),
        temporal_edf_near_upper = direction_value("time", "near_upper"),
        truth_cz_250 = truth_peak,
        realized_truth_cz_250 = realized_truth_peak,
        bgam_cz_250 = bgam_peak,
        peak_ratio = bgam_peak / truth_peak,
        FP = bgam_metrics$FP,
        FPR = bgam_metrics$FPR,
        recovery = projection_recovery(theta_hat, theta),
        leakage = spatial_leakage(theta_hat, active),
        n_divergent = n_divergent,
        n_treedepth = n_treedepth,
        max_rhat = max(draws_summary$rhat, na.rm = TRUE),
        min_ess_bulk = min(draws_summary$ess_bulk, na.rm = TRUE),
        min_ess_tail = min(draws_summary$ess_tail, na.rm = TRUE),
        stringsAsFactors = FALSE
    )
    write.csv(summary, file.path(out_dir, "diagnostic_summary.csv"),
              row.names = FALSE)

    residual <- plot_empirical_minus_bgam(
        sim_model, posterior_map, radius_rule,
        file.path(out_dir, "empirical_minus_bgam_250ms.png")
    )
    write.csv(residual,
              file.path(out_dir, "empirical_minus_bgam_250ms.csv"),
              row.names = FALSE)
    print(summary)
    print(method_metrics)
    invisible(summary)
}

if (!interactive() && sys.nframe() == 0) {
    args <- commandArgs(trailingOnly = TRUE)
    radius_rule <- if (length(args)) args[[1]] else "1sd"
    mode <- if (length(args) >= 2L) args[[2]] else "all-methods"
    if (!mode %in% c("all-methods", "bgam-only")) {
        stop("Mode must be all-methods or bgam-only", call. = FALSE)
    }
    run_paired_radius_bgam_diagnostics(
        radius_rule,
        bgam_only = identical(mode, "bgam-only")
    )
}
