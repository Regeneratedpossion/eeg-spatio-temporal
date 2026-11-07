###########################################################
# 26_coverage_bgam.R
#
# Coverage of the 95% credible interval of the condition effect g at every
# electrode-time cell, against the population truth of the scenario.
#
# Run:
#   Rscript code/26_coverage_bgam.R <fit.rds> <sim.rds> <gt.rds> <out_dir>
#
# Writes to <out_dir>:
#   coverage_cells.csv    one row per cell: truth, posterior mean, 2.5% and 97.5% quantiles
#   coverage_summary.csv  coverage in three groups of cells
###########################################################

suppressPackageStartupMessages({
    library(brms)
})

project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = getwd()))

truth_file <- file.path(project_root, "code", c(
    "ground_truth.R", "paired_radius_superpopulation_truth.R"
))
source(truth_file[file.exists(truth_file)][1])

coverage_bgam_fit <- function(fit_file, sim_file, gt_file, out_dir,
                              true_onset = 0.160, true_offset = 0.340) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    fit <- readRDS(fit_file)
    gt <- readRDS(gt_file)
    sim <- readRDS(sim_file)

    sensors <- unique(sim[, c("channel", "xproj", "yproj")])
    sensors <- sensors[order(sensors$channel), ]
    channels <- sensors$channel
    times <- sort(unique(fit$data$time))

    # Same grid and population-level draws of g as in 23_evaluate_bgam_fit.R.
    grid <- expand.grid(time = times, channel = channels,
                        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    grid <- merge(grid, sensors, by = "channel", sort = FALSE)
    grid <- grid[order(match(grid$channel, channels), grid$time), ]
    grid$eeg_sem <- 1
    levs <- levels(fit$data$condition)
    nd1 <- transform(grid, condition = factor(levs[1], levels = levs))
    nd2 <- transform(grid, condition = factor(levs[2], levels = levs))
    draws <- posterior_epred(fit, newdata = nd2, re_formula = NA) -
        posterior_epred(fit, newdata = nd1, re_formula = NA)

    template <- file.path(project_root, "resources", "brms_meeg", "code",
                          "erp_template.R")
    truth_df <- paired_radius_superpopulation_truth(gt, template)
    truth <- truth_df$effect[match(paste(grid$channel, grid$time),
                                   paste(truth_df$channel, truth_df$time))]
    stopifnot(!anyNA(truth))

    q <- apply(draws, 2, quantile, probs = c(0.025, 0.975))
    cells <- data.frame(
        channel = grid$channel, time = grid$time, truth = truth,
        post_mean = colMeans(draws), lo = q[1, ], hi = q[2, ]
    )
    cells$covered <- cells$lo <= cells$truth & cells$truth <= cells$hi

    in_window <- cells$time >= true_onset - 1e-9 & cells$time <= true_offset + 1e-9
    active_ch <- cells$channel %in% gt$active_electrodes
    cells$group <- ifelse(active_ch & in_window, "active",
                   ifelse(!active_ch & in_window, "inactive_in_window", "outside_window"))
    write.csv(cells, file.path(out_dir, "coverage_cells.csv"), row.names = FALSE)

    summ <- do.call(rbind, lapply(c("active", "inactive_in_window", "outside_window", "all"), function(g) {
        sel <- if (g == "all") rep(TRUE, nrow(cells)) else cells$group == g
        data.frame(group = g, n_cells = sum(sel), coverage = mean(cells$covered[sel]),
                   mean_width = mean((cells$hi - cells$lo)[sel]))
    }))
    write.csv(summ, file.path(out_dir, "coverage_summary.csv"), row.names = FALSE)
    invisible(summ)
}

if (sys.nframe() == 0L) {
    a <- commandArgs(trailingOnly = TRUE)
    stopifnot(length(a) == 4L)
    coverage_bgam_fit(a[1], a[2], a[3], a[4])
}
