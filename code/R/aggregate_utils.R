# Input/output helpers of the aggregation scripts.
# They find files and folders only; the statistics stay in
# aggregate_metrics.R.

# Results folder of one scenario.
scenario_results_dir <- function(res_dir, scenario) {
    file.path(
        res_dir,
        sprintf("single_node_cond%s_ds_cr_seeds001_100", scenario)
    )
}

# Fits that pass the convergence criterion (condition, seed).
read_passed_fits <- function() {
    read.table(
        here::here("resources/passed_fits.txt"),
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
        file.path(out_dir, "tables"),
        showWarnings = FALSE,
        recursive = TRUE
    )

    out_dir
}
