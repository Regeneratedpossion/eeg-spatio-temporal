# Coverage and calibration of the posterior of g, by scenario and region (passed fits only).
# Usage: Rscript code/extra/05_agg_coverage.R
# Reads results/scenario*/seed_*/evaluation/coverage_z.csv (made by code/extra/04_coverage_z.R).
# Writes coverage_by_region.csv to results/fig_data (AGG_OUT overrides).
#   coverage: share of fits whose 90% interval contains the true g
#   z_mean:   mean of (posterior mean - true g) / posterior SD (bias)
#   z_sd:     SD of z over fits (above 1: intervals too narrow)
# Each statistic is computed per cell over the fits of a scenario, then averaged within the region.
suppressPackageStartupMessages(library(dplyr))
source(here::here("code/R/aggregate_utils.R"))

man <- read_passed_fits()
res_dir <- here::here("results")
out_dir <- aggregate_out_dir()
win <- c(0.160, 0.340)

# ring: 0 active, 1 adjacent, 2 and 3 distant
region_of <- function(ring, time) {
    in_win <- time >= win[1] - 1e-9 & time <= win[2] + 1e-9
    ifelse(!in_win, "outside window",
        ifelse(ring == 0L, "active",
            ifelse(ring == 1L, "adjacent", "distant")))
}

summaries <- lapply(sort(unique(man$cond)), function(sc) {
    d <- scenario_results_dir(res_dir, sc)
    files <- file.path(d, sprintf("seed_%03d", man$seed[man$cond == sc]), "evaluation", "coverage_z.csv")
    x <- bind_rows(lapply(files, read.csv))
    x$covered <- x$truth >= x$lo & x$truth <= x$hi
    x$z <- (x$mean - x$truth) / x$sd
    x$region <- region_of(x$ring, x$time)

    per_cell <- x |>
        group_by(region, channel, time) |>
        summarise(coverage = mean(covered), z_mean = mean(z), z_sd = sd(z),
                  width = mean(hi - lo), .groups = "drop")

    bind_rows(
        per_cell |> summarise(region = "all cells", n_cells = n(), across(c(coverage, z_mean, z_sd, width), mean)),
        per_cell |> group_by(region) |> summarise(n_cells = n(), across(c(coverage, z_mean, z_sd, width), mean))
    ) |> mutate(scenario = sc, n_fits = length(files), .before = 1)
})

res <- bind_rows(summaries)
write.csv(res, file.path(out_dir, "coverage_by_region.csv"), row.names = FALSE)
print(as.data.frame(mutate(res, across(where(is.numeric) & !c(scenario, n_fits, n_cells), ~ round(.x, 3)))), row.names = FALSE)
