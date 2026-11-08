###########################################################
# 04_evaluate_bgam.R
#
# Evaluate one BGAM fit with the thesis definitions (Sections 2.2.2 and 2.2.4).
#
# Run:
#   Rscript code/04_evaluate/04_evaluate_bgam.R <fit.rds> <sim.rds> <gt.rds> <out_dir>
#
# Writes to <out_dir>:
#   evaluation_main.csv         one row: main rule (global SESOI 0-0.10 s, odds > 10)
#   evaluation_sensitivity.csv  all thresholds x SESOI variants
#   sesoi.csv                   SESOI values used
#   convergence_effect.csv      R-hat and bulk/tail ESS of g at every cell
###########################################################

suppressPackageStartupMessages({
    library(brms)
    library(dplyr)
    library(igraph)
})

project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = here::here()))

source(file.path(project_root, "code", "R", "evaluation_functions.R"))

# ---------------------------------------------------------------
# Main
# ---------------------------------------------------------------

evaluate_bgam_fit <- function(fit_file, sim_file, gt_file, out_dir,
                              s = evaluation_settings()) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    fit <- readRDS(fit_file)
    gt <- readRDS(gt_file)
    sim <- readRDS(sim_file)

    sensors <- unique(sim[, c("channel", "xproj", "yproj")])
    sensors <- sensors[order(sensors$channel), ]
    channels <- sensors$channel
    times <- sort(unique(fit$data$time))
    adj <- readRDS(file.path(project_root, "resources", "adjacency",
                             "adjacency_xproj_yproj_61ch.rds"))
    adj <- adj[match(channels, rownames(adj)), match(channels, rownames(adj))]

    # Posterior draws of the condition difference on the time x channel grid.
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
    est <- colMeans(draws)

    # Convergence of the evaluated quantity g at every cell (Section 2.2.5).
    conv_g <- effect_convergence(draws, posterior::nchains(fit))
    write.csv(cbind(grid[, c("channel", "time")], conv_g),
              file.path(out_dir, "convergence_effect.csv"), row.names = FALSE)

    # Truth: theoretical population field and the binary active set.
    template <- file.path(project_root, "resources", "brms_meeg", "code",
                          "erp_template.R")
    truth_df <- paired_radius_superpopulation_truth(gt, template)
    truth <- truth_df$effect[match(paste(grid$channel, grid$time),
                                   paste(truth_df$channel, truth_df$time))]
    stopifnot(!anyNA(truth))
    in_window <- grid$time >= s$true_onset - 1e-9 & grid$time <= s$true_offset + 1e-9
    active_ch <- grid$channel %in% gt$active_electrodes
    in_A <- active_ch & in_window
    null_window <- !active_ch & in_window
    ring_chan <- electrode_rings(adj, channels, channels %in% gt$active_electrodes)
    ring_grid <- ring_chan[match(grid$channel, channels)]
    to_matrix <- function(x) matrix(x, nrow = length(times), ncol = length(channels))

    # SESOI from the posterior mean in the baseline window.
    base <- grid$time >= 0 & grid$time <= s$main_baseline_end + 1e-9
    sesoi <- rbind(
        data.frame(rule = "zero", baseline_end = NA, channel = channels,
                   delta = 0),
        data.frame(rule = "global", baseline_end = s$main_baseline_end,
                   channel = channels, delta = max(0, max(est[base])))
    )
    write.csv(sesoi, file.path(out_dir, "sesoi.csv"), row.names = FALSE)

    rule_grid <- unique(sesoi[, c("rule", "baseline_end")])
    topo <- abs(grid$time - s$topo_time) < 1e-9
    rows <- list()
    for (r in seq_len(nrow(rule_grid))) {
        rr <- rule_grid[r, ]
        sel <- sesoi$rule == rr$rule & (is.na(rr$baseline_end) & is.na(sesoi$baseline_end) |
                                        sesoi$baseline_end %in% rr$baseline_end)
        delta <- sesoi$delta[sel][match(grid$channel, sesoi$channel[sel])]
        p <- colMeans(sweep(draws, 2, delta, ">"))
        odds <- pmax(p, 1 / nrow(draws)) / pmax(1 - p, 1 / nrow(draws))
        for (o in s$odds_grid) {
            det <- odds > o
            if (isTRUE(rr$rule == "global") && isTRUE(abs(rr$baseline_end - s$main_baseline_end) < 1e-9) &&
                o == s$main_odds) elec_det <- det
            oo <- extract_onset_offset(to_matrix(det), times, adj, s)
            rows[[length(rows) + 1L]] <- data.frame(
                rule = rr$rule, baseline_end = rr$baseline_end, odds = o,
                t(classification(det, in_A)),
                t(setNames(classification(det[topo], active_ch[topo]),
                           paste0(c("TP", "FP", "FN", "TN", "precision", "recall",
                                    "F1", "IoU", "FPR"), "_250"))),
                FPR_window = mean(det[null_window]),
                FPR_ring1 = mean(det[in_window & ring_grid == 1L]),
                FPR_ring2 = mean(det[in_window & ring_grid == 2L]),
                FPR_ring3 = mean(det[in_window & ring_grid == 3L]),
                post_offset = mean(det[active_ch & grid$time > s$true_offset + 1e-9]),
                onset = oo[["onset"]], offset = oo[["offset"]],
                onset_error_ms = 1000 * (oo[["onset"]] - s$true_onset),
                offset_error_ms = 1000 * (oo[["offset"]] - s$true_offset),
                largest_tied = oo[["tied"]]
            )
        }
    }
    # Proportion of detected cells in the active window, per electrode (main rule).
    by_elec <- tapply(elec_det[in_window], grid$channel[in_window], mean)
    write.csv(data.frame(channel = names(by_elec), rate = as.numeric(by_elec)),
              file.path(out_dir, "detection_by_electrode.csv"), row.names = FALSE)
    sens <- do.call(rbind, rows)
    # FPR: inactive electrodes x active window (Section 2.2.4); FPR_cells: whole grid.
    sens$FPR_cells <- sens$FPR
    sens$FPR <- sens$FPR_window
    sens$FPR_window <- NULL
    write.csv(sens, file.path(out_dir, "evaluation_sensitivity.csv"), row.names = FALSE)

    # Continuous metrics do not depend on the decision rule.
    cont <- continuous_metrics(est, truth, !in_A)
    # Leakage sums energy over the active window before dividing (Section 2.2.4).
    cont[["leakage"]] <- sum(est[in_window & !active_ch]^2) / sum(est[in_window]^2)
    cont_250 <- continuous_metrics(est[topo], truth[topo], !active_ch[topo])
    names(cont_250) <- paste0(names(cont_250), "_250")
    main <- sens[sens$rule == "global" & sens$baseline_end %in% s$main_baseline_end &
                 sens$odds == s$main_odds, ]
    main <- cbind(main, t(cont), t(cont_250),
                  sesoi_global = sesoi$delta[sesoi$rule == "global" &
                                             sesoi$baseline_end %in% s$main_baseline_end][1])
    main$max_rhat_g <- max(conv_g$rhat, na.rm = TRUE)
    main$min_ess_bulk_g <- min(conv_g$ess_bulk, na.rm = TRUE)
    main$min_ess_tail_g <- min(conv_g$ess_tail, na.rm = TRUE)
    main$n_unstable_decisions <- unstable_decisions(
        draws, rep(main$sesoi_global, ncol(draws)), s$main_odds,
        posterior::nchains(fit)
    )
    write.csv(main, file.path(out_dir, "evaluation_main.csv"), row.names = FALSE)
    invisible(main)
}

if (!interactive() && sys.nframe() == 0L) {
    a <- commandArgs(trailingOnly = TRUE)
    if (length(a) != 4L) {
        stop("Usage: Rscript code/04_evaluate/04_evaluate_bgam.R <fit.rds> <sim.rds> <gt.rds> <out_dir>")
    }
    main <- evaluate_bgam_fit(a[[1]], a[[2]], a[[3]], a[[4]])
    print(main[, c("F1", "IoU", "FPR", "onset_error_ms", "offset_error_ms",
                   "recovery", "leakage", "map_cor", "sesoi_global")])
}
