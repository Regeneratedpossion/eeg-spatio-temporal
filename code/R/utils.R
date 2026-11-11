###########################################################
# utils.R — shared helpers across the project
###########################################################

suppressPackageStartupMessages(library(dplyr))

# Plotting reference threshold, in microvolts.
PLOT_SESOI_UV <- 0.10

# Prefer PROJECT_ROOT; otherwise use the folder that holds the .here marker file.
project_root <- function() {
    p <- Sys.getenv("PROJECT_ROOT", unset = "")
    if (nzchar(p)) return(normalizePath(p, mustWork = FALSE))
    normalizePath(here::here(), mustWork = FALSE)
}

# Results folder of one scenario: FIXED_GRID_BATCH, or scenario<COND> with COND = 1, 2, 3 or 4.
fixed_grid_batch <- function() {
    batch <- Sys.getenv("FIXED_GRID_BATCH", unset = "")
    cond <- Sys.getenv("COND", unset = "")
    if (nzchar(batch)) return(batch)
    if (!cond %in% c("1", "2", "3", "4")) stop("Set COND to 1, 2, 3 or 4", call. = FALSE)
    sprintf("scenario%s", cond)
}

data_path <- function(...)   file.path(project_root(), "data", ...)
resource_path <- function(...) file.path(project_root(), "resources", ...)
brms_meeg_path <- function(...) {
    external_path <- Sys.getenv("BRMS_MEEG", unset = "")
    if (nzchar(external_path)) {
        return(file.path(external_path, ...))
    }
    resource_path("brms_meeg/code", ...)
}

# Aggregate trials to participant x channel x condition x time.
# For model fitting, use `eeg_mean | se(eeg_sem)`: the SEM is the
# measurement error of the cell mean. `eeg_sd` is retained for diagnostics.
aggregate_trials <- function(sim_data) {
    sim_data %>%
        group_by(participant, channel, xproj, yproj, condition, time) %>%
        summarise(
            eeg_mean = mean(eeg),
            eeg_sd = sd(eeg),
            R = dplyr::n(),
            .groups = "drop"
        ) %>%
        mutate(
            eeg_sem = eeg_sd / sqrt(R),
            condition = factor(condition, levels = c("cond1", "cond2"))
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
