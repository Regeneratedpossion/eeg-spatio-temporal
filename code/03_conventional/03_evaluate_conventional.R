###########################################################
# 03_evaluate_conventional.R
#
# Evaluate the conventional methods on one dataset with the same truth,
# extraction rule and metrics as the BGAM (code/04_evaluate/04_evaluate_bgam.R).
#
# Run:
#   Rscript code/03_conventional/03_evaluate_conventional.R <sim.rds> <gt.rds> <out_dir>
#
# Writes <out_dir>/evaluation_conventional.csv: one row per method and alpha.
###########################################################

suppressPackageStartupMessages({
    library(dplyr)
    library(reticulate)
})

project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = here::here()))
# Shared truth, extraction rule and metric functions.
source(file.path(project_root, "code", "R", "evaluation_functions.R"))

conventional_settings <- function() {
    list(alphas = c(0.05, 0.01, 0.005, 0.001), main_alpha = 0.05,
         n_perm = 1024L, perm_seed = 42L, tfce_start = 0.0, tfce_step = 0.2)
}

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

evaluate_conventional <- function(sim_file, gt_file, out_dir,
                                  s = evaluation_settings(),
                                  cs = conventional_settings(),
                                  spatial_adjacency = TRUE) {
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    sim <- readRDS(sim_file)
    gt <- readRDS(gt_file)
    channels <- sort(unique(sim$channel))
    times <- sort(unique(sim$time))
    adj <- readRDS(file.path(project_root, "resources", "adjacency",
                             "adjacency_xproj_yproj_61ch.rds"))
    adj <- adj[match(channels, rownames(adj)), match(channels, rownames(adj))]

    # Truth on the same time x channel layout as the masks.
    in_window <- times >= s$true_onset - 1e-9 & times <= s$true_offset + 1e-9
    active_ch <- channels %in% gt$active_electrodes
    ring_chan <- electrode_rings(adj, channels, active_ch)
    truth_A <- outer(in_window, active_ch, "&")
    topo <- which(abs(times - s$topo_time) < 1e-9)

    arr <- difference_array(sim, times, channels)

    # Cluster mass and TFCE use this adjacency; FALSE keeps only the time neighbours.
    mne_adj <- if (spatial_adjacency) adj else adj * 0

    cell_p <- apply(arr, c(2, 3), function(x) t.test(x)$p.value)
    p_list <- list(
        uncorrected = cell_p,
        BH = matrix(p.adjust(cell_p, "BH"), nrow(cell_p)),
        BY = matrix(p.adjust(cell_p, "BY"), nrow(cell_p)),
        Holm = matrix(p.adjust(cell_p, "holm"), nrow(cell_p)),
        cluster_mass = mne_cell_pvalues(arr, mne_adj, cs, tfce = FALSE),
        TFCE = mne_cell_pvalues(arr, mne_adj, cs, tfce = TRUE)
    )

    rows <- list()
    by_elec_rows <- list()
    for (m in names(p_list)) for (a in cs$alphas) {
        mask <- p_list[[m]] < a
        if (a == cs$main_alpha) {
            by_elec_rows[[length(by_elec_rows) + 1L]] <- data.frame(
                method = m, channel = channels, rate = colMeans(mask[in_window, , drop = FALSE]))
        }
        oo <- extract_onset_offset(mask, times, adj, s)
        cells <- classification(as.vector(mask), as.vector(truth_A))
        topo_m <- classification(mask[topo, ], active_ch)
        rows[[length(rows) + 1L]] <- data.frame(
            method = m, alpha = a, main = a == cs$main_alpha,
            t(cells[c("TP", "FP", "FN", "TN", "precision", "recall", "F1", "IoU")]),
            FPR = mean(mask[in_window, !active_ch]), FPR_250 = topo_m[["FPR"]],
            FPR_ring1 = mean(mask[in_window, ring_chan == 1L]),
            FPR_ring2 = mean(mask[in_window, ring_chan == 2L]),
            FPR_ring3 = mean(mask[in_window, ring_chan == 3L]),
            FPR_cells = cells[["FPR"]],
            F1_250 = topo_m[["F1"]], IoU_250 = topo_m[["IoU"]],
            onset = oo[["onset"]], offset = oo[["offset"]],
            onset_error_ms = 1000 * (oo[["onset"]] - s$true_onset),
            offset_error_ms = 1000 * (oo[["offset"]] - s$true_offset),
            largest_tied = oo[["tied"]]
        )
    }
    out <- do.call(rbind, rows)
    write.csv(do.call(rbind, by_elec_rows), file.path(out_dir, "detection_by_electrode_conventional.csv"),
              row.names = FALSE)
    write.csv(out, file.path(out_dir, "evaluation_conventional.csv"), row.names = FALSE)
    invisible(out)
}

if (!interactive() && sys.nframe() == 0L) {
    a <- commandArgs(trailingOnly = TRUE)
    if (length(a) != 3L) {
        stop("Usage: Rscript code/03_conventional/03_evaluate_conventional.R <sim.rds> <gt.rds> <out_dir>")
    }
    out <- evaluate_conventional(a[[1]], a[[2]], a[[3]])
    print(out[out$main, c("method", "F1", "IoU", "FPR", "onset_error_ms",
                          "offset_error_ms")], row.names = FALSE, digits = 3)
}
