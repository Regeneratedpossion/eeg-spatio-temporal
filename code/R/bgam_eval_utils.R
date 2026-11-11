# Input, grid and output helpers of 04_evaluate_bgam.R.
# The statistical definitions stay in the main script.


# 1. Inputs --------------------------------------------------

# Reads the fit, simulated data and ground truth of one seed, and
# returns the sensors, channels, time points and adjacency.
load_bgam_eval_inputs <- function(fit_file, sim_file, gt_file,
                                  project_root) {

    fit <- readRDS(fit_file)
    gt <- readRDS(gt_file)
    sim <- readRDS(sim_file)

    sensors <- unique(sim[, c("channel", "xproj", "yproj")])
    sensors <- sensors[order(sensors$channel), ]
    channels <- sensors$channel
    times <- sort(unique(fit$data$time))

    adj <- readRDS(file.path(
        project_root, "resources", "adjacency",
        "adjacency_xproj_yproj_61ch.rds"
    ))
    adj <- adj[
        match(channels, rownames(adj)),
        match(channels, rownames(adj))
    ]

    list(
        fit = fit, gt = gt, sensors = sensors,
        channels = channels, times = times, adj = adj
    )
}


# 2. Prediction grid -----------------------------------------

# Time x channel grid and the two newdata sets (condition 1 and 2)
# for posterior_epred(). eeg_sem is fixed at 1.
prepare_bgam_prediction_grid <- function(fit, sensors, channels, times) {

    grid <- expand.grid(
        time = times, channel = channels,
        KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
    )
    grid <- merge(grid, sensors, by = "channel", sort = FALSE)
    grid <- grid[order(match(grid$channel, channels), grid$time), ]
    grid$eeg_sem <- 1

    levs <- levels(fit$data$condition)

    list(
        grid = grid,
        nd1 = transform(
            grid, condition = factor(levs[1], levels = levs)
        ),
        nd2 = transform(
            grid, condition = factor(levs[2], levels = levs)
        )
    )
}


# 3. Outputs -------------------------------------------------

# Writes the five evaluation CSV files of one fit.
save_bgam_eval_outputs <- function(out_dir, conv_effect, sesoi,
                                   by_elec, sensitivity, main) {

    write.csv(
        conv_effect,
        file.path(out_dir, "convergence_effect.csv"),
        row.names = FALSE
    )

    write.csv(
        sesoi,
        file.path(out_dir, "sesoi.csv"),
        row.names = FALSE
    )

    write.csv(
        data.frame(
            channel = names(by_elec), rate = as.numeric(by_elec)
        ),
        file.path(out_dir, "detection_by_electrode.csv"),
        row.names = FALSE
    )

    write.csv(
        sensitivity,
        file.path(out_dir, "evaluation_sensitivity.csv"),
        row.names = FALSE
    )

    write.csv(
        main,
        file.path(out_dir, "evaluation_main.csv"),
        row.names = FALSE
    )

    invisible(NULL)
}
