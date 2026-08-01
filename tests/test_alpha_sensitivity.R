script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[[1]])
project_root <- normalizePath(file.path(dirname(script_path), ".."))
Sys.setenv(PROJECT_ROOT = project_root)

source(file.path(project_root, "code", "27_significance_alpha_sensitivity.R"))

channels <- c("A", "B", "C")
pvalue_cells <- data.frame(
    method = "raw_p",
    channel = c("A", "B", "C", "A", "B", "C"),
    time = c(0, 0, 0, 1, 1, 1),
    p_value = c(0.0005, 0.009, 0.049, 0.02, 0.8, 0.9),
    practical_null = c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE)
)
truth <- pvalue_cells[c("channel", "time")]
truth$true_active <- c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE)
adjacency <- matrix(
    c(0, 1, 0,
      1, 0, 1,
      0, 1, 0),
    nrow = 3,
    dimnames = list(channels, channels),
    byrow = TRUE
)

got <- compute_alpha_sensitivity(
    pvalue_cells = pvalue_cells,
    truth_cells = truth,
    adjacency = adjacency,
    alphas = c(0.05, 0.01, 0.005, 0.001),
    min_cluster_size = 2L,
    min_onset_channels = 1L,
    min_onset_run = 1L
)

stopifnot(identical(got$alpha, c(0.05, 0.01, 0.005, 0.001)))
stopifnot(identical(got$n_pred, c(4L, 2L, 0L, 0L)))
stopifnot(all(diff(got$n_pred) <= 0))

at_001 <- got[got$alpha == 0.01, ]
stopifnot(at_001$st_tp == 2L)
stopifnot(at_001$st_fp == 0L)
stopifnot(at_001$st_fn == 2L)
stopifnot(isTRUE(all.equal(at_001$st_iou, 0.5, tolerance = 1e-12)))
stopifnot(isTRUE(all.equal(at_001$st_F1, 2 / 3, tolerance = 1e-12)))
stopifnot(isTRUE(all.equal(at_001$spatial_F1, 0.8, tolerance = 1e-12)))

cat("alpha-sensitivity regression tests passed\n")
