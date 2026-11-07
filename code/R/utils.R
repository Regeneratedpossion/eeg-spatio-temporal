###########################################################
# utils.R — shared helpers across the project
###########################################################

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(brms)
    library(tidybayes)
    library(posterior)
})

select <- dplyr::select
filter <- dplyr::filter

# Shared thresholds.
ROPE_UV           <- 0.10
ODDS_THRESHOLD    <- 10
BIAS_TOLERANCE_MS <- 20

# Prefer PROJECT_ROOT; otherwise use the folder that holds the .here marker file.
project_root <- function() {
    p <- Sys.getenv("PROJECT_ROOT", unset = "")
    if (nzchar(p)) return(normalizePath(p, mustWork = FALSE))
    normalizePath(here::here(), mustWork = FALSE)
}

data_path <- function(...)   file.path(project_root(), "data", ...)
fit_path  <- function(...)   file.path(project_root(), "results/fits", ...)
diag_path <- function(...)   file.path(project_root(), "results/diagnostics", ...)
fig_path  <- function(...)   file.path(project_root(), "results/figures", ...)
resource_path <- function(...) file.path(project_root(), "resources", ...)
brms_meeg_path <- function(...) {
    external_path <- Sys.getenv("BRMS_MEEG", unset = "")
    if (nzchar(external_path)) {
        return(file.path(external_path, ...))
    }
    resource_path("brms_meeg/code", ...)
}
adjacency_path <- function(...) resource_path("adjacency", ...)
source(file.path(project_root(), "code", "R", "superpopulation_truth.R"))

resolve_path <- function(base, path) {
    if (grepl("^(/|[A-Za-z]:)", path)) path else file.path(base, path)
}

basis_settings <- function(default_xy = 16L, default_t = 12L) {
    list(
        k_xy = as.integer(Sys.getenv("K_XY", default_xy)),
        k_t = as.integer(Sys.getenv("K_T", default_t)),
        full = FALSE
    )
}

slide_inputs <- function() {
    data_subdir <- Sys.getenv("SLIDE_DATA_SUBDIR", "montecarlo_n20_hard")
    seed <- as.integer(Sys.getenv("SLIDE_SEED", "1"))
    sim_default <- file.path(data_subdir, sprintf("sim_%03d.rds", seed))
    gt_default <- file.path(data_subdir, sprintf("gt_%03d.rds", seed))

    list(
        out_dir = fig_path(Sys.getenv("SLIDE_OUT_SUBDIR",
                                      "n20_hard_sem_slides")),
        sim_file = resolve_path(data_path(), Sys.getenv("SLIDE_SIM_FILE", sim_default)),
        gt_file = resolve_path(data_path(), Sys.getenv("SLIDE_GT_FILE", gt_default)),
        fit_file = resolve_path(
            fit_path(),
            Sys.getenv("SLIDE_FIT_FILE", "bgam_fit_v5_k20_sem.rds")
        )
    )
}

mc_paths <- function() {
    scenario <- Sys.getenv("MC_SCENARIO", unset = "low")
    default_subdir <- paste0("montecarlo_p20_", scenario, "_heterogeneity")
    data_subdir <- Sys.getenv("MC_DATA_SUBDIR", unset = default_subdir)

    list(
        prefix = Sys.getenv("MC_AGG_PREFIX", unset = data_subdir),
        data_dir = data_path(data_subdir),
        compare_dir = file.path(project_root(), "results/comparisons",
                                Sys.getenv("MC_COMPARE_SUBDIR",
                                           unset = paste0(data_subdir, "_sem"))),
        summary_dir = file.path(project_root(), "results",
                                Sys.getenv("MC_SUMMARY_SUBDIR",
                                           unset = paste0(data_subdir, "_sem"))),
        out_dir = file.path(project_root(), "results/comparisons")
    )
}

montecarlo_files <- function(seed, paths = mc_paths()) {
    list(
        sim = file.path(paths$data_dir, sprintf("sim_%03d.rds", seed)),
        gt = file.path(paths$data_dir, sprintf("gt_%03d.rds", seed)),
        summary = file.path(paths$summary_dir, sprintf("bgam_summary_%03d.rds", seed)),
        csv = file.path(paths$compare_dir, sprintf("sim_%03d.csv", seed))
    )
}

ppc_inputs <- function() {
    data_subdir <- Sys.getenv("PPC_DATA_SUBDIR", "montecarlo_n20_hard")
    seed <- as.integer(Sys.getenv("PPC_SEED", "1"))
    sim_default <- file.path(data_subdir, sprintf("sim_%03d.rds", seed))

    list(
        sim_file = resolve_path(data_path(), Sys.getenv("PPC_SIM_FILE", sim_default)),
        fig_dir = fig_path()
    )
}

setup_cmdstan_cache <- function() {
    cache_dir <- Sys.getenv(
        "CMDSTAN_CACHE_DIR",
        unset = file.path(project_root(), "cmdstan_cache")
    )
    dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
    options(cmdstanr_write_stan_file_dir = cache_dir)
    invisible(cache_dir)
}

# Aggregate trials to participant x channel x condition x time.
# For model fitting, use `eeg_mean | se(eeg_sem)`: the SEM is the
# measurement error of the cell mean. `eeg_sd` is retained for diagnostics.
aggregate_trials <- function(sim_data) {
    sim_data %>%
        group_by(participant, channel, xproj, yproj, condition, time) %>%
        summarise(eeg_mean = mean(eeg),
                  eeg_sd   = sd(eeg),
                  R        = dplyr::n(),
                  eeg_sem  = sd(eeg) / sqrt(dplyr::n()),
                  .groups  = "drop") %>%
        mutate(condition = factor(condition, levels = c("cond1", "cond2")))
}

# Theoretical population truth (thesis Eq. 19) for the requested cells.
continuous_truth_cells <- function(ground_truth, channels = NULL, times = NULL) {
    truth <- paired_radius_superpopulation_truth(
        ground_truth,
        template_path = brms_meeg_path("erp_template.R")
    )

    requested_channels <- if (is.null(channels)) {
        unique(as.character(truth$channel))
    } else {
        as.character(channels)
    }
    requested_times <- if (is.null(times)) unique(truth$time) else times

    requested <- expand.grid(
        channel = requested_channels,
        time = requested_times,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )
    truth$channel <- as.character(truth$channel)
    matched <- merge(
        requested,
        truth,
        by = c("channel", "time"),
        all.x = TRUE,
        sort = FALSE
    )
    if (anyNA(matched$effect)) {
        stop("Missing continuous ground-truth cells", call. = FALSE)
    }

    matched <- matched[order(
        match(matched$channel, requested_channels),
        match(matched$time, requested_times)
    ), c("channel", "time", "effect"), drop = FALSE]
    rownames(matched) <- NULL
    matched
}

continuous_truth_summary <- function(ground_truth, channel, peak_time,
                                     window_times) {
    cells <- continuous_truth_cells(
        ground_truth,
        channels = channel,
        times = unique(c(peak_time, window_times))
    )

    list(
        peak = cells$effect[cells$time == peak_time],
        window_mean = mean(cells$effect[cells$time %in% window_times])
    )
}

continuous_truth_matrix <- function(ground_truth, times, channels) {
    cells <- continuous_truth_cells(
        ground_truth,
        channels = channels,
        times = times
    )
    matrix(
        cells$effect,
        nrow = length(times),
        ncol = length(channels),
        dimnames = list(as.character(times), as.character(channels))
    )
}

# Inverse-distance interpolation for topomaps.
idw_interp <- function(grid_xy, src_xy, src_v, power = 2, k = 6) {
    n_g <- nrow(grid_xy)
    out <- numeric(n_g)
    for (i in seq_len(n_g)) {
        d   <- sqrt((src_xy$xproj - grid_xy$xproj[i])^2 +
                    (src_xy$yproj - grid_xy$yproj[i])^2)
        idx <- order(d)[1:min(k, length(d))]
        w   <- 1 / (d[idx]^power + 1e-9)
        out[i] <- sum(w * src_v[idx]) / sum(w)
    }
    out
}

pick_time <- function(target_t, all_times) {
    all_times[which.min(abs(all_times - target_t))]
}

# Find contiguous 1D significant intervals.
find_clusters <- function(time, value, threshold) {
    data.frame(time = time, above = value > threshold) %>%
        mutate(cluster_id = cumsum(above & !lag(above, default = FALSE))) %>%
        filter(above) %>%
        group_by(cluster_id) %>%
        summarise(cluster_onset  = first(time),
                  cluster_offset = last(time),
                  .groups = "drop") %>%
        data.frame()
}

# Head outline used by topomap plots.
scalp_decorations <- function(head_radius) {
    list(
        head_circle = data.frame(
            x = head_radius * cos(seq(0, 2 * pi, length.out = 200)),
            y = head_radius * sin(seq(0, 2 * pi, length.out = 200))
        ),
        nose = data.frame(
            x = c(-0.08, 0, 0.08) * head_radius,
            y = c( 1.00, 1.10, 1.00) * head_radius
        ),
        ear_left = {
            theta <- seq(-pi / 2, pi / 2, length.out = 40)
            data.frame(
                x = (-1.0 - 0.06 * cos(theta)) * head_radius,
                y = ( 0.12 * sin(theta))       * head_radius
            )
        },
        ear_right = {
            theta <- seq(-pi / 2, pi / 2, length.out = 40)
            data.frame(
                x =  (1.0 + 0.06 * cos(theta)) * head_radius,
                y = ( 0.12 * sin(theta))       * head_radius
            )
        }
    )
}
