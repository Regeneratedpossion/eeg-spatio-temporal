script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[[1]])
project_root <- normalizePath(file.path(dirname(script_path), ".."))
source(file.path(
    project_root, "code", "truth", "migrate_ground_truth.R"
))

analysis_times <- seq(0, 0.5, by = 0.01)
original <- list(
    true_onset = 0.20,
    true_offset = 0.30,
    temporal_threshold = 0.55,
    temporal_active_times = seq(0.20, 0.30, by = 0.01),
    continuous_effect = expand.grid(
        channel = c("CZ", "PZ"),
        time = analysis_times,
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    ),
    active_electrodes = c("CZ", "PZ"),
    active_cells = expand.grid(
        channel = c("CZ", "PZ"),
        time = seq(0.20, 0.30, by = 0.01),
        KEEP.OUT.ATTRS = FALSE,
        stringsAsFactors = FALSE
    ),
    marker = "must remain unchanged"
)
original$continuous_effect$effect <- seq_len(nrow(original$continuous_effect))

repaired <- repair_temporal_truth(original)
expected_times <- seq(0.16, 0.34, by = 0.01)
stopifnot(identical(repaired$true_onset, 0.160))
stopifnot(identical(repaired$true_offset, 0.342))
stopifnot(is.null(repaired$temporal_threshold))
stopifnot(isTRUE(all.equal(
    repaired$temporal_active_times,
    expected_times,
    tolerance = sqrt(.Machine$double.eps)
)))
stopifnot(identical(repaired$continuous_effect, original$continuous_effect))
stopifnot(identical(repaired$marker, original$marker))
stopifnot(nrow(repaired$active_cells) == 2L * length(expected_times))

work_dir <- tempfile("temporal-truth-test-")
backup_dir <- tempfile("temporal-truth-backup-")
dir.create(work_dir)
dir.create(backup_dir)
on.exit(unlink(c(work_dir, backup_dir), recursive = TRUE), add = TRUE)
truth_path <- file.path(work_dir, "gt_001.rds")
saveRDS(original, truth_path)
repair_ground_truth_file(truth_path, backup_dir)
stopifnot(identical(
    readRDS(file.path(backup_dir, basename(work_dir), "gt_001.rds")),
    original
))
stopifnot(identical(readRDS(truth_path), repaired))

cat("paired-radius temporal-truth repair tests passed\n")
