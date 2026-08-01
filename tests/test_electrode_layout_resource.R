script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[[1]])
repository_root <- normalizePath(file.path(dirname(script_path), ".."))
Sys.setenv(PROJECT_ROOT = repository_root)

source(file.path(repository_root, "code", "01_simulate_eeg_data.R"))

# The generator must obtain coordinates from a formal resource and remain
# usable when no previously generated pilot dataset exists.
fixture_root <- tempfile("electrode-layout-resource-")
fixture_dir <- file.path(fixture_root, "resources", "electrodes")
dir.create(fixture_dir, recursive = TRUE)
fixture <- data.frame(
    channel = c("CZ", "PZ"),
    xproj = c(0, 0.25),
    yproj = c(0, -0.50)
)
write.csv(
    fixture,
    file.path(fixture_dir, "electrode_layout_61ch.csv"),
    row.names = FALSE
)

old_root <- Sys.getenv("PROJECT_ROOT")
on.exit({
    Sys.setenv(PROJECT_ROOT = old_root)
    unlink(fixture_root, recursive = TRUE)
}, add = TRUE)
Sys.setenv(PROJECT_ROOT = fixture_root)

layout <- load_electrode_layout()
stopifnot(identical(layout, fixture))

cat("electrode layout resource test passed\n")
