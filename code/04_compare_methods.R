###########################################################
# 04_compare_methods.R
#
# Compare methods on the full time x channel grid.
# Each method creates a TRUE/FALSE mask; connected blobs give onset/offset.
###########################################################

source(file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()), "code/utils.R"))

library(reticulate)
library(igraph)

set.seed(42)

dir.create(file.path(project_root(), "results/comparisons"), showWarnings = FALSE, recursive = TRUE)

# Main comparison settings. Change these here to rerun the benchmark differently.
compare_method_settings <- function() {
    list(
        alpha = 0.1,
        n_perm = 1024L,
        perm_seed = 42L,
        tfce_start = 0.0,
        tfce_step = 0.2,
        total_draws = 8000L,
        min_cluster_size = 2L,
        min_onset_channels = 3L,
        min_onset_run = 2L
    )
}

settings <- compare_method_settings()

ALPHA <- settings$alpha
N_PERM <- settings$n_perm
PERM_SEED <- settings$perm_seed
TFCE_START <- settings$tfce_start
TFCE_STEP <- settings$tfce_step
TOTAL_DRAWS <- settings$total_draws
ODDS_FLOOR <- 1 / TOTAL_DRAWS
MIN_CLUSTER_SIZE <- settings$min_cluster_size
MIN_ONSET_CHANNELS <- settings$min_onset_channels
MIN_ONSET_RUN <- settings$min_onset_run

# MNE is optional; cluster and TFCE return NA if MNE is unavailable.
load_mne <- function() {
    tryCatch(import("mne"), error = function(e) NULL)
}

mne <- load_mne()

# ---------------------------------------------------------------
# Connected components on the time x channel mask
# ---------------------------------------------------------------

# Find connected TRUE blobs using time adjacency and electrode adjacency.
.find_st_components <- function(mask, adj_chan) {
    n_times    <- nrow(mask)
    n_channels <- ncol(mask)

    # TRUE cells become graph nodes.
    idx_true <- which(mask, arr.ind = TRUE)
    n_active <- nrow(idx_true)
    if (n_active == 0L) return(NULL)

    # Map each time-channel cell to its graph node id.
    lin <- (idx_true[, 2] - 1L) * n_times + idx_true[, 1]
    pos <- integer(n_times * n_channels)
    pos[lin] <- seq_len(n_active)

    # Time edges: same channel, next time point.
    time_from <- which(idx_true[, 1] < n_times)
    time_to <- pos[(idx_true[time_from, 2] - 1L) * n_times +
                   idx_true[time_from, 1] + 1L]
    keep_time <- time_to > 0L

    edges_from <- time_from[keep_time]
    edges_to <- time_to[keep_time]

    # Space edges: same time point, neighbouring channels.
    channel_pairs <- which(upper.tri(adj_chan) & adj_chan == 1L, arr.ind = TRUE)
    if (nrow(channel_pairs) > 0L) {
        for (p in seq_len(nrow(channel_pairs))) {
            c1 <- channel_pairs[p, 1]
            c2 <- channel_pairs[p, 2]
            node1 <- pos[(c1 - 1L) * n_times + seq_len(n_times)]
            node2 <- pos[(c2 - 1L) * n_times + seq_len(n_times)]
            both_active <- node1 > 0L & node2 > 0L

            edges_from <- c(edges_from, node1[both_active])
            edges_to <- c(edges_to, node2[both_active])
        }
    }

    # Connected components are the detected time-channel blobs.
    g <- make_empty_graph(n = n_active, directed = FALSE)
    if (length(edges_from) > 0L) {
        g <- add_edges(g, rbind(edges_from, edges_to))
    }
    list(idx_true = idx_true, cl = components(g))
}


# Keep only connected blobs with at least min_size cells.
.large_components <- function(mask, adj_chan, min_size = MIN_CLUSTER_SIZE) {
    # Find all blobs first, then mark which cells belong to large blobs.
    info <- .find_st_components(mask, adj_chan)
    if (is.null(info)) return(NULL)

    keep <- info$cl$csize[info$cl$membership] >= min_size
    list(info = info, keep = keep)
}

# Remove connected blobs smaller than min_size.
.filter_components <- function(mask, adj_chan, min_size = MIN_CLUSTER_SIZE) {
    # Return the original empty mask if no TRUE blob exists.
    comp <- .large_components(mask, adj_chan, min_size)
    if (is.null(comp)) return(mask)

    # Rebuild the mask using only cells from large blobs.
    out  <- matrix(FALSE, nrow(mask), ncol(mask))
    out[comp$info$idx_true[comp$keep, ]] <- TRUE
    out
}

# ---------------------------------------------------------------
# Onset / offset from the largest connected blob
# ---------------------------------------------------------------

.extract_onset_offset_2d <- function(mask, times, adj_chan,
                                      min_size = MIN_CLUSTER_SIZE,
                                      min_channels = MIN_ONSET_CHANNELS,
                                      min_run = MIN_ONSET_RUN) {
    # Work only with blobs large enough to count as a detection.
    comp <- .large_components(mask, adj_chan, min_size)
    if (is.null(comp) || !any(comp$keep)) {
        return(list(onset = NA_real_, offset = NA_real_))
    }

    # Use the largest detected blob.
    cl <- comp$info$cl
    large_ids <- unique(cl$membership[comp$keep])
    largest <- large_ids[which.max(cl$csize[large_ids])]
    time_idx <- comp$info$idx_true[cl$membership == largest, 1]
    n_detected <- tabulate(time_idx, nbins = length(times))

    # Onset/offset require enough active channels for enough consecutive times.
    runs <- rle(n_detected >= min_channels)
    good <- which(runs$values & runs$lengths >= min_run)
    if (length(good) == 0L) {
        return(list(onset = NA_real_, offset = NA_real_))
    }

    starts <- cumsum(c(1L, head(runs$lengths, -1L)))
    ends <- cumsum(runs$lengths)

    list(onset  = times[starts[good[1]]],
         offset = times[ends[tail(good, 1)]])
}

# ---------------------------------------------------------------
# Data conversion
# ---------------------------------------------------------------

# Convert long data to arr[participant, time, channel] = cond2 - cond1.
build_diff_array <- function(model_data) {
    chs   <- sort(unique(model_data$channel))
    ts    <- sort(unique(model_data$time))
    ppts  <- sort(unique(model_data$participant))

    diff_long <- model_data %>%
        select(participant, channel, time, condition, eeg_mean) %>%
        pivot_wider(names_from = condition, values_from = eeg_mean) %>%
        mutate(diff = cond2 - cond1) %>%
        mutate(
            ppt_idx = match(participant, ppts),
            time_idx = match(time, ts),
            channel_idx = match(channel, chs)
        )

    arr <- array(NA_real_,
                 dim = c(length(ppts), length(ts), length(chs)),
                 dimnames = list(ppts, as.character(ts), chs))

    # Fill arr[participant, time, channel] with cond2 - cond1.
    arr[cbind(diff_long$ppt_idx,
              diff_long$time_idx,
              diff_long$channel_idx)] <- diff_long$diff

    list(arr = arr, channels = chs, times = ts)
}

# ---------------------------------------------------------------
# Method masks
# ---------------------------------------------------------------

# Send arr to MNE cluster/TFCE and return a time x channel p-value matrix.
run_st_mne <- function(arr, adj_use, mne, tfce = FALSE,
                       n_perm = N_PERM,
                       perm_seed = PERM_SEED,
                       tfce_start = TFCE_START,
                       tfce_step = TFCE_STEP) {
    if (is.null(mne)) {
        return(matrix(NA_real_, dim(arr)[2], dim(arr)[3]))
    }
    sp <- import("scipy.sparse")
    adj_sparse <- sp$csr_matrix(matrix(as.double(adj_use),
                                        nrow = nrow(adj_use)))

    args <- list(arr,
                 adjacency      = adj_sparse,
                 n_permutations = n_perm,
                 n_jobs         = 1L,
                 seed           = perm_seed,
                 out_type       = "mask")
    if (tfce) args$threshold <- dict(start = tfce_start, step = tfce_step)
    out <- do.call(mne$stats$spatio_temporal_cluster_1samp_test, args)
    clusters  <- out[[2]]
    pvs       <- out[[3]]
    pmat <- matrix(1, dim(arr)[2], dim(arr)[3])
    for (k in seq_along(clusters)) {
        idx <- which(as.matrix(clusters[[k]]), arr.ind = TRUE)
        pmat[idx] <- pmin(pmat[idx], pvs[k])
    }
    pmat
}

# Convert BGAM posterior probabilities to a time x channel odds matrix.
.bgam_full_to_odds <- function(bgam_full, channels, times,
                                odds_floor = ODDS_FLOOR) {
    odds_long <- bgam_full %>%
        mutate(odds = pmax(p_sig, odds_floor) /
                      pmax(1 - p_sig, odds_floor))

    odds_mat <- odds_long %>%
        select(channel, time, odds) %>%
        pivot_wider(names_from = channel, values_from = odds) %>%
        arrange(time) %>%
        select(-time) %>%
        as.matrix()
    odds_mat[, channels, drop = FALSE]
}

# Safe t-test p-value for one time-channel cell.
ttest_p <- function(x) {
    tryCatch(t.test(x)$p.value, error = function(e) NA_real_)
}

# Convert p-values to a clean TRUE/FALSE mask.
p_mask <- function(p, alpha = ALPHA) {
    !is.na(p) & p < alpha
}

# Convert posterior odds to a clean TRUE/FALSE mask.
odds_mask <- function(odds, threshold = ODDS_THRESHOLD) {
    !is.na(odds) & odds > threshold
}

# Build one time x channel TRUE/FALSE mask for each comparison method.
build_method_masks <- function(arr, channels, times, adj_chan, bgam_full,
                               mne = NULL, alpha = ALPHA) {
    n_times <- length(times)
    n_chan <- length(channels)

    cell_pvals <- apply(arr, c(2, 3), ttest_p)
    cluster_p <- run_st_mne(arr, adj_chan, mne = mne, tfce = FALSE)
    tfce_p <- run_st_mne(arr, adj_chan, mne = mne, tfce = TRUE)
    odds_mat <- .bgam_full_to_odds(bgam_full, channels, times)

    list(
        raw_p = p_mask(cell_pvals, alpha),
        BH = p_mask(matrix(p.adjust(cell_pvals, "BH"), n_times, n_chan), alpha),
        BY = p_mask(matrix(p.adjust(cell_pvals, "BY"), n_times, n_chan), alpha),
        holm = p_mask(matrix(p.adjust(cell_pvals, "holm"), n_times, n_chan), alpha),
        cluster = p_mask(cluster_p, alpha),
        tfce = p_mask(tfce_p, alpha),
        BGAM_2D = .filter_components(odds_mask(odds_mat), adj_chan,
                                      min_size = MIN_CLUSTER_SIZE)
    )
}

# ===============================================================
# compare_methods_2d()
#
# Compare all methods for one single simulation on the full time x channel grid.
# Each method creates a detection mask; the largest connected blob gives
# onset and offset.
#
# Inputs:
#   sim_data     trial-level data
#   ground_truth true onset/offset
#   adj_chan     electrode adjacency matrix
#   bgam_full    BGAM posterior summary with channel, time, p_sig
#   mne          MNE Python module; NULL skips cluster/tfce
#   alpha        p-value threshold
#
# Returns:
#   per-simulation output, one row per method:
#     method, onset, offset, bias_onset_ms, bias_offset_ms
#
# Next aggregation step:
#   05_aggregate_results.R summarises Bias/MAE/RMSE/SD, F1_ST,
#   and TP/FP/FN diagnostics across simulations.
# ===============================================================

compare_methods_2d <- function(sim_data, ground_truth,
                                adj_chan,
                                bgam_full,
                                mne = NULL,
                                alpha = 0.1) {

    stopifnot(!is.null(bgam_full))

    sim_data_model <- aggregate_trials(sim_data)
    true_on        <- ground_truth$true_onset
    true_off       <- ground_truth$true_offset

    bundle   <- build_diff_array(sim_data_model)
    arr      <- bundle$arr
    channels <- bundle$channels
    times    <- bundle$times

    # Re-order adjacency rows/cols to match channel ordering.
    if (!identical(rownames(adj_chan), channels)) {
        ord <- match(channels, rownames(adj_chan))
        adj_chan <- adj_chan[ord, ord]
    }

    masks <- build_method_masks(arr, channels, times, adj_chan, bgam_full,
                                mne = mne, alpha = alpha)

    # Extract onset/offset from the largest component for every method.
    bind_rows(lapply(names(masks), function(nm) {
        oo <- .extract_onset_offset_2d(masks[[nm]], times, adj_chan,
                                        min_size = MIN_CLUSTER_SIZE)
        data.frame(
            method         = nm,
            onset          = oo$onset,
            offset         = oo$offset,
            bias_onset_ms  = (oo$onset  - true_on)  * 1000,
            bias_offset_ms = (oo$offset - true_off) * 1000
        )
    }))
}
