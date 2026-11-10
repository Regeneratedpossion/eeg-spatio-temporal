# Control analysis: cluster-mass and TFCE without spatial neighbours.
# Temporal adjacency and evaluation metrics remain unchanged.
#
# Usage:
# Rscript code/03_conventional/03_evaluate_conventional_no_space.R <sim.rds> <gt.rds> <out_dir>
#
# Use a separate output directory from the main analysis.

root <- normalizePath(
    Sys.getenv("PROJECT_ROOT", unset = here::here())
)

source(file.path(
    root, "code/03_conventional/03_evaluate_conventional.R"
))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3L) {
    stop(
        "Usage: Rscript code/03_conventional/",
        "03_evaluate_conventional_no_space.R ",
        "<sim.rds> <gt.rds> <out_dir>",
        call. = FALSE
    )
}

out <- evaluate_conventional(
    sim_file = args[[1]],
    gt_file = args[[2]],
    out_dir = args[[3]],
    spatial_adjacency = FALSE
)

print(
    out[out$main, c(
        "method", "F1", "onset_error_ms", "offset_error_ms"
    )],
    row.names = FALSE,
    digits = 3
)
