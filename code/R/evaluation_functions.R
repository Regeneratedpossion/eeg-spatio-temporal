# Functions shared by the evaluation of the BGAM (code/04_evaluate) and of the conventional methods (code/03_conventional):
# the fixed evaluation settings, the onset/offset extraction, the classification metrics, the electrode rings and the checks of g.
suppressPackageStartupMessages({
    library(dplyr)
    library(igraph)
})

evaluation_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = here::here()))

source(file.path(evaluation_root, "code", "R", "settings.R"))
source(file.path(evaluation_root, "code", "R", "superpopulation_truth.R"))

# Fixed evaluation settings from the thesis.
evaluation_settings <- function() {
    list(
        true_onset = 0.160,
        true_offset = 0.340,
        main_odds = MAIN_ODDS,
        odds_grid = c(1, 3, 6, 10, 20, 50, 100),
        main_baseline_end = 0.10,
        topo_time = 0.25
    )
}

# ---------------------------------------------------------------
# Onset/offset extraction (Section 2.2.4), used for the BGAM and the conventional methods
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
    if (is.null(comp)) {
        return(c(onset = NA_real_, offset = NA_real_, tied = 0))
    }
    largest <- which.max(comp$cl$csize)
    # 1 if another component has the same size as the largest one.
    tied <- as.numeric(sum(comp$cl$csize == max(comp$cl$csize)) > 1)
    time_idx <- comp$idx[comp$cl$membership == largest, 1]
    c(onset = times[min(time_idx)], offset = times[max(time_idx)], tied = tied)
}

# Ring of each electrode around the active set: 0 = active, 1 = adjacent to an
# active electrode, 2 = adjacent to ring 1, 3 = all others.
electrode_rings <- function(adj, channels, active_ch) {
    ring <- ifelse(active_ch, 0L, NA_integer_)
    for (r in 1:2) {
        prev <- which(ring == r - 1L)
        nb <- which(colSums(adj[prev, , drop = FALSE]) > 0 & is.na(ring))
        ring[nb] <- r
    }
    ring[is.na(ring)] <- 3L
    ring
}

# ---------------------------------------------------------------
# Metrics (Section 2.2.4)
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
# Convergence of the condition effect g (Section 2.2.5)
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
