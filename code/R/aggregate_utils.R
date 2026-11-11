# Input/output helpers of the aggregation scripts.
# They find files and folders only; the statistics stay in
# aggregate_metrics.R.

# Results folder of one scenario.
scenario_results_dir <- function(res_dir, scenario) {
    file.path(
        res_dir,
        sprintf("scenario%s", scenario)
    )
}

# Fits that pass the convergence criterion (condition, seed).
read_passed_fits <- function() {
    read.table(
        here::here("results/passed_fits.txt"),
        col.names = c("cond", "seed")
    )
}

# Output folder of the summaries (AGG_OUT overrides it).
aggregate_out_dir <- function() {
    out_dir <- Sys.getenv(
        "AGG_OUT",
        unset = here::here("results", "fig_data")
    )

    dir.create(
        file.path(dirname(out_dir), "tables"),
        showWarnings = FALSE,
        recursive = TRUE
    )

    out_dir
}
