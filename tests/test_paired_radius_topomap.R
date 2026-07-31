script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[[1]])
project_root <- normalizePath(file.path(dirname(script_path), ".."))
Sys.setenv(PROJECT_ROOT = project_root)

source(file.path(
    project_root, "code", "reporting", "07_plot_paired_radius_truth_pilot.R"
))

sensors <- data.frame(
    channel = c("CZ", "PZ"),
    xproj = c(0, 1),
    yproj = c(0, 0)
)
truth <- list(
    active_electrodes = "CZ",
    continuous_effect = data.frame(
        channel = rep(c("CZ", "PZ"), each = 2),
        time = rep(c(0.24, 0.25), 2),
        effect = c(1, 2, 0, 0)
    )
)

electrodes <- truth_topomap_electrodes(truth, sensors, 0.25)
stopifnot(nrow(electrodes) == 2L)
stopifnot(all(electrodes$time == 0.25))
stopifnot(electrodes$value[electrodes$channel == "CZ"] == 2)
stopifnot(electrodes$active[electrodes$channel == "CZ"])
stopifnot(!electrodes$active[electrodes$channel == "PZ"])

cat("paired radius topomap tests passed\n")
