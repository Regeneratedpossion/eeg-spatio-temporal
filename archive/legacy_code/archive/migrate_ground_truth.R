repair_temporal_truth <- function(ground_truth,
                                  true_onset = 0.160,
                                  true_offset = 0.342) {
    required <- c("continuous_effect", "active_electrodes", "active_cells")
    missing <- setdiff(required, names(ground_truth))
    if (length(missing)) {
        stop(
            "Ground truth is missing required fields: ",
            paste(missing, collapse = ", "),
            call. = FALSE
        )
    }
    if (!all(c("time", "channel", "effect") %in%
             names(ground_truth$continuous_effect))) {
        stop("continuous_effect is missing time, channel, or effect",
             call. = FALSE)
    }

    allowed_changes <- c(
        "true_onset", "true_offset", "temporal_threshold",
        "temporal_active_times", "active_cells"
    )
    stable_names <- setdiff(names(ground_truth), allowed_changes)
    stable_fields <- ground_truth[stable_names]

    analysis_times <- sort(unique(ground_truth$continuous_effect$time))
    temporal_active_times <- analysis_times[
        analysis_times >= true_onset & analysis_times <= true_offset
    ]
    if (!length(temporal_active_times)) {
        stop("The fixed temporal window contains no analysis time points",
             call. = FALSE)
    }

    ground_truth$true_onset <- true_onset
    ground_truth$true_offset <- true_offset
    ground_truth$temporal_threshold <- NULL
    ground_truth$temporal_active_times <- temporal_active_times
    ground_truth$active_cells <- expand.grid(
        channel = ground_truth$active_electrodes,
        time = temporal_active_times,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    )

    stopifnot(identical(ground_truth[stable_names], stable_fields))
    ground_truth
}

repair_ground_truth_file <- function(path, backup_dir) {
    if (!file.exists(path)) {
        stop("Ground-truth file does not exist: ", path, call. = FALSE)
    }
    if (!dir.exists(backup_dir)) {
        stop("Backup directory does not exist: ", backup_dir, call. = FALSE)
    }

    backup_path <- file.path(backup_dir, basename(dirname(path)), basename(path))
    dir.create(dirname(backup_path), recursive = TRUE, showWarnings = FALSE)
    if (!file.copy(path, backup_path, overwrite = FALSE, copy.date = TRUE)) {
        stop("Could not create backup: ", backup_path, call. = FALSE)
    }

    repaired <- repair_temporal_truth(readRDS(path))
    temporary_path <- tempfile(
        pattern = paste0(basename(path), "."),
        tmpdir = dirname(path)
    )
    on.exit(unlink(temporary_path), add = TRUE)
    saveRDS(repaired, temporary_path)
    stopifnot(identical(readRDS(temporary_path), repaired))
    if (!file.rename(temporary_path, path)) {
        stop("Could not atomically replace: ", path, call. = FALSE)
    }
    invisible(path)
}

if (!interactive() && sys.nframe() == 0L) {
    project_root <- Sys.getenv("PROJECT_ROOT", unset = getwd())
    backup_dir <- Sys.getenv("TEMPORAL_TRUTH_BACKUP_DIR", unset = "")
    if (!nzchar(backup_dir)) {
        stop("TEMPORAL_TRUTH_BACKUP_DIR must name an existing directory",
             call. = FALSE)
    }
    paths <- file.path(
        project_root,
        "data",
        c("montecarlo_r1sd", "montecarlo_rhalfmax", "montecarlo_r2sd"),
        "gt_001.rds"
    )
    vapply(paths, repair_ground_truth_file, character(1),
           backup_dir = backup_dir)
}
