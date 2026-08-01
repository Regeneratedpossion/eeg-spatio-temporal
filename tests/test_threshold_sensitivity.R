script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[[1]])
project_root <- normalizePath(file.path(dirname(script_path), ".."))
Sys.setenv(PROJECT_ROOT = project_root)

source(file.path(project_root, "code", "26_posterior_odds_sensitivity.R"))

channels <- c("A", "B", "C")
times <- c(0, 1)
cells <- data.frame(
    channel = c("A", "B", "C", "A", "B", "C"),
    time = c(0, 0, 0, 1, 1, 1),
    odds = c(100, 50, 20, 6, 0, 0),
    practical_null = c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE)
)

truth <- cells[c("channel", "time")]
truth$true_active <- c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE)

adjacency <- matrix(
    c(0, 1, 0,
      1, 0, 1,
      0, 1, 0),
    nrow = 3,
    dimnames = list(channels, channels),
    byrow = TRUE
)

got <- compute_odds_sensitivity(
    cells = cells,
    truth_cells = truth,
    adjacency = adjacency,
    thresholds = c(1, 10, 50),
    min_cluster_size = 2L,
    min_onset_channels = 1L,
    min_onset_run = 1L
)

stopifnot(identical(got$odds_threshold, c(1, 10, 50)))
stopifnot(identical(got$n_pred, c(4L, 3L, 0L)))
stopifnot(all(diff(got$n_pred) <= 0))

at_10 <- got[got$odds_threshold == 10, ]
stopifnot(at_10$st_tp == 3L)
stopifnot(at_10$st_fp == 0L)
stopifnot(at_10$st_fn == 1L)
stopifnot(isTRUE(all.equal(at_10$st_iou, 0.75, tolerance = 1e-12)))
stopifnot(isTRUE(all.equal(at_10$st_F1, 6 / 7, tolerance = 1e-12)))
stopifnot(at_10$spatial_F1 == 1)
stopifnot(at_10$estimated_onset == 0)
stopifnot(at_10$estimated_offset == 0)

cat("threshold-sensitivity regression tests passed\n")
