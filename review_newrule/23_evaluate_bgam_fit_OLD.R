###########################################################
# 23_evaluate_bgam_fit.R
#
# Evaluate one BGAM fit with the thesis definitions (Sections 2.3, 3.1, 3.3).
#
# Run:
#   Rscript code/23_evaluate_bgam_fit.R <fit.rds> <sim.rds> <gt.rds> <out_dir>
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

project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = getwd()))

# The population truth function lives in ground_truth.R (GitHub) or
# paired_radius_superpopulation_truth.R (LIDO).
truth_file <- file.path(project_root, "code", c(
    "ground_truth.R", "paired_radius_superpopulation_truth.R"
))
source(truth_file[file.exists(truth_file)][1])

# Fixed evaluation settings from the thesis.
evaluation_settings <- function() {
    list(
        true_onset = 0.160,
        true_offset = 0.340,
        main_odds = 10,
        odds_grid = c(1, 3, 6, 10, 20, 50, 100),
        baseline_ends = c(0.10, 0.15),
        main_baseline_end = 0.10,
        topo_time = 0.25,
        min_cluster_size = 2L,
        min_onset_channels = 3L,
        min_onset_run = 2L
    )
}

# ---------------------------------------------------------------
# Onset/offset extraction (Section 3.2.3), same rule as 04_compare_methods.R
# ---------------------------------------------------------------

find_components <- function(mask, adj) {
    n_times <- nrow(mask)
    idx <- which(mask, arr.ind = TRUE)
    if (nrow(idx) == 0L) return(NULL)
    pos <- integer(length(mask))
    pos[(idx[, 2] - 1L) * n_times + idx[, 1]] <- seq_len(nrow(idx))
    from <- which(idx[, 1] < n_times)
    to <- pos[(idx[from, 2] - 1L) * n_times + idx[from, 1] + 1L]
    e_from <- from[to > 0L]
    e_to <- to[to > 0L]
    pairs <- which(upper.tri(adj) & adj == 1L, arr.ind = TRUE)
    for (p in seq_len(nrow(pairs))) {
        n1 <- pos[(pairs[p, 1] - 1L) * n_times + seq_len(n_times)]
        n2 <- pos[(pairs[p, 2] - 1L) * n_times + seq_len(n_times)]
        both <- n1 > 0L & n2 > 0L
        e_from <- c(e_from, n1[both])
        e_to <- c(e_to, n2[both])
    }
    g <- make_empty_graph(n = nrow(idx), directed = FALSE)
    if (length(e_from)) g <- add_edges(g, rbind(e_from, e_to))
    list(idx = idx, cl = components(g))
}

extract_onset_offset <- function(mask, times, adj, s) {
    comp <- find_components(mask, adj)
    if (is.null(comp)) return(c(onset = NA_real_, offset = NA_real_))
    keep <- comp$cl$csize[comp$cl$membership] >= s$min_cluster_size
    if (!any(keep)) return(c(onset = NA_real_, offset = NA_real_))
    ids <- unique(comp$cl$membership[keep])
    largest <- ids[which.max(comp$cl$csize[ids])]
    n_detected <- tabulate(comp$idx[comp$cl$membership == largest, 1],
                           nbins = length(times))
    runs <- rle(n_detected >= s$min_onset_channels)
    good <- which(runs$values & runs$lengths >= s$min_onset_run)
    if (!length(good)) return(c(onset = NA_real_, offset = NA_real_))
    starts <- cumsum(c(1L, head(runs$lengths, -1L)))
    ends <- cumsum(runs$lengths)
    c(onset = times[starts[good[1]]], offset = times[ends[tail(good, 1)]])
}

# ---------------------------------------------------------------
# Metrics (Section 3.3)
# ---------------------------------------------------------------

classification <- function(pred, truth) {
    tp <- sum(pred & truth); fp <- sum(pred & !truth)
    fn <- sum(!pred & truth); tn <- sum(!pred & !truth)
    c(TP = tp, FP = fp, FN = fn, TN = tn,
      precision = if (tp + fp) tp / (tp + fp) else NA_real_,
      recall = if (tp + fn) tp / (tp + fn) else NA_real_,
      F1 = if (2 * tp + fp + fn) 2 * tp / (2 * tp + fp + fn) else NA_real_,
      IoU = if (tp + fp + fn) tp / (tp + fp + fn) else NA_real_,
      FPR = if (fp + tn) fp / (fp + tn) else NA_real_)
}

continuous_metrics <- function(est, truth, outside) {
    c(map_cor = cor(est, truth),
      RMSE = sqrt(mean((est - truth)^2)),
      recovery = sum(est * truth) / sum(truth^2),
      leakage = sum(est[outside]^2) / sum(est^2))
}

# ---------------------------------------------------------------
# Convergence of the condition effect g (Section 2.4.4)
# ---------------------------------------------------------------

# R-hat and bulk/tail ESS of g at each cell; draws are stored chain by chain.
effect_convergence <- function(draws, n_chains) {
    n_iter <- nrow(draws) / n_chains
    stopifnot(n_iter == round(n_iter))
    per_cell <- function(j) {
        x <- matrix(draws[, j], nrow = n_iter, ncol = n_chains)
        c(rhat = posterior::rhat(x), ess_bulk = posterior::ess_bulk(x),
          ess_tail = posterior::ess_tail(x))
    }
    as.data.frame(t(vapply(seq_len(ncol(draws)), per_cell, numeric(3))))
}

# Number of cells whose detection could flip within 2 MCSE of the odds threshold.
unstable_decisions <- function(draws, delta, odds, n_chains) {
    p_cut <- odds / (1 + odds)
    n_iter <- nrow(draws) / n_chains
    exceed <- sweep(draws, 2, delta, ">")
    p <- colMeans(exceed)
    near <- which(abs(p - p_cut) < 0.05)
    mcse <- vapply(near, function(j) {
        posterior::mcse_mean(matrix(as.numeric(exceed[, j]), n_iter, n_chains))
    }, numeric(1))
    sum(abs(p[near] - p_cut) < 2 * mcse, na.rm = TRUE)
}

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

    # Convergence of the evaluated quantity g at every cell (Section 2.4.4).
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
    to_matrix <- function(x) matrix(x, nrow = length(times), ncol = length(channels))

    # SESOI from the posterior mean in the baseline window.
    sesoi <- do.call(rbind, lapply(s$baseline_ends, function(b) {
        base <- grid$time >= 0 & grid$time <= b + 1e-9
        el <- tapply(est[base], grid$channel[base], max)
        rbind(data.frame(rule = "global", baseline_end = b, channel = channels,
                         delta = max(0, max(est[base]))),
              data.frame(rule = "electrode", baseline_end = b, channel = channels,
                         delta = pmax(0, as.numeric(el[channels]))))
    }))
    sesoi <- rbind(data.frame(rule = "zero", baseline_end = NA, channel = channels,
                              delta = 0), sesoi)
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
            oo <- extract_onset_offset(to_matrix(det), times, adj, s)
            rows[[length(rows) + 1L]] <- data.frame(
                rule = rr$rule, baseline_end = rr$baseline_end, odds = o,
                t(classification(det, in_A)),
                t(setNames(classification(det[topo], active_ch[topo]),
                           paste0(c("TP", "FP", "FN", "TN", "precision", "recall",
                                    "F1", "IoU", "FPR"), "_250"))),
                FPR_window = mean(det[null_window]),
                onset = oo[["onset"]], offset = oo[["offset"]],
                onset_error_ms = 1000 * (oo[["onset"]] - s$true_onset),
                offset_error_ms = 1000 * (oo[["offset"]] - s$true_offset)
            )
        }
    }
    sens <- do.call(rbind, rows)
    # FPR: inactive electrodes x active window (Section 3.3.2); FPR_cells: whole grid.
    sens$FPR_cells <- sens$FPR
    sens$FPR <- sens$FPR_window
    sens$FPR_window <- NULL
    write.csv(sens, file.path(out_dir, "evaluation_sensitivity.csv"), row.names = FALSE)

    # Continuous metrics do not depend on the decision rule.
    cont <- continuous_metrics(est, truth, !in_A)
    # Leakage sums energy over the active window before dividing (Section 3.3).
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
        stop("Usage: Rscript code/23_evaluate_bgam_fit.R <fit.rds> <sim.rds> <gt.rds> <out_dir>")
    }
    main <- evaluate_bgam_fit(a[[1]], a[[2]], a[[3]], a[[4]])
    print(main[, c("F1", "IoU", "FPR", "onset_error_ms", "offset_error_ms",
                   "recovery", "leakage", "map_cor", "sesoi_global")])
}
