# Helpers of the conventional methods (used by 03_evaluate_conventional.R
# and 04_save_masks.R).

suppressPackageStartupMessages({
    library(dplyr)
    library(reticulate)
})

# Participant-level condition differences as arr[participant, time, channel].
difference_array <- function(sim, times, channels) {
    cell <- sim %>%
        group_by(participant, channel, time, condition) %>%
        summarise(eeg_mean = mean(eeg), .groups = "drop")
    wide <- tidyr::pivot_wider(cell, names_from = condition, values_from = eeg_mean)
    ppts <- sort(unique(wide$participant))
    arr <- array(NA_real_, dim = c(length(ppts), length(times), length(channels)))
    arr[cbind(match(wide$participant, ppts),
              match(round(wide$time, 6), round(times, 6)),
              match(wide$channel, channels))] <- wide$cond2 - wide$cond1
    stopifnot(!anyNA(arr))
    arr
}

# Cluster-mass or TFCE p-values per cell from MNE.
mne_cell_pvalues <- function(arr, adj, s, tfce) {
    mne <- import("mne")
    sp <- import("scipy.sparse")
    args <- list(arr, adjacency = sp$csr_matrix(matrix(as.double(adj), nrow = nrow(adj))),
                 n_permutations = s$n_perm, n_jobs = 1L, seed = s$perm_seed,
                 out_type = "mask")
    if (tfce) args$threshold <- dict(start = s$tfce_start, step = s$tfce_step)
    out <- do.call(mne$stats$spatio_temporal_cluster_1samp_test, args)
    pmat <- matrix(1, dim(arr)[2], dim(arr)[3])
    for (k in seq_along(out[[2]])) {
        idx <- which(as.matrix(out[[2]][[k]]), arr.ind = TRUE)
        pmat[idx] <- pmin(pmat[idx], out[[3]][k])
    }
    pmat
}
