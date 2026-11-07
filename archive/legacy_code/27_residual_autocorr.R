###########################################################
# 27_residual_autocorr.R
#
# Residual correlation of one BGAM fit: lag-1 autocorrelation in time within
# each participant x electrode x condition series, and the correlation between
# neighbouring electrodes at the same time point.
#
# Run:
#   Rscript code/27_residual_autocorr.R <fit.rds> <sim.rds> <out_dir>
###########################################################

suppressPackageStartupMessages(library(brms))

project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = getwd()))

residual_correlation <- function(fit_file, sim_file, out_dir, n_draws = 400L) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    fit <- readRDS(fit_file)
    d <- fit$data
    stopifnot("eeg_mean" %in% names(d))

    # The fit data has coordinates only, so electrode names come from the simulated data.
    sensors <- unique(readRDS(sim_file)[, c("channel", "xproj", "yproj")])
    d$channel <- sensors$channel[match(paste(d$xproj, d$yproj), paste(sensors$xproj, sensors$yproj))]
    stopifnot(!anyNA(d$channel))

    # Posterior mean fit with the participant effects, from a subset of draws.
    est <- fitted(fit, summary = TRUE, ndraws = n_draws)[, "Estimate"]
    d$res <- d$eeg_mean - est
    d <- d[order(d$participant, d$channel, d$condition, d$time), ]

    # Lag-1 autocorrelation in time within each series.
    key <- interaction(d$participant, d$channel, d$condition, drop = TRUE)
    lag1 <- vapply(split(d$res, key), function(r) {
        if (length(r) < 3L || sd(r) == 0) return(NA_real_)
        cor(r[-1], r[-length(r)])
    }, numeric(1))

    # Correlation between neighbouring electrodes at the same time point.
    adj <- readRDS(file.path(project_root, "resources", "adjacency",
                             "adjacency_xproj_yproj_61ch.rds"))
    wide <- reshape(d[, c("participant", "condition", "time", "channel", "res")],
                    idvar = c("participant", "condition", "time"),
                    timevar = "channel", direction = "wide")
    colnames(wide) <- sub("^res\\.", "", colnames(wide))
    chs <- intersect(rownames(adj), colnames(wide))
    pairs <- which(upper.tri(adj[chs, chs]) & adj[chs, chs] == 1, arr.ind = TRUE)
    space <- vapply(seq_len(nrow(pairs)), function(i) {
        cor(wide[[chs[pairs[i, 1]]]], wide[[chs[pairs[i, 2]]]], use = "complete.obs")
    }, numeric(1))

    out <- data.frame(
        n_series = sum(!is.na(lag1)),
        lag1_mean = mean(lag1, na.rm = TRUE),
        lag1_median = median(lag1, na.rm = TRUE),
        lag1_q25 = unname(quantile(lag1, 0.25, na.rm = TRUE)),
        lag1_q75 = unname(quantile(lag1, 0.75, na.rm = TRUE)),
        n_pairs = length(space),
        space_mean = mean(space, na.rm = TRUE),
        space_median = median(space, na.rm = TRUE),
        res_sd = sd(d$res)
    )
    write.csv(out, file.path(out_dir, "residual_autocorr.csv"), row.names = FALSE)
    invisible(out)
}

if (sys.nframe() == 0L) {
    a <- commandArgs(trailingOnly = TRUE)
    stopifnot(length(a) == 3L)
    residual_correlation(a[1], a[2], a[3])
}
