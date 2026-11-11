# Run from the project root. Reads results/scenario*/seed_*/evaluation; writes to results/fig_data (AGG_OUT overrides).
source(here::here("code/R/aggregate_utils.R"))
res_dir <- here::here("results")
out_dir <- aggregate_out_dir()
out <- list()
for (c in c("1", "3", "2", "4")) {
  d <- scenario_results_dir(res_dir, c)
  b <- do.call(rbind, lapply(Sys.glob(file.path(d, "seed_*/evaluation/detection_by_electrode.csv")), read.csv)); b$method <- "BGAM"
  g <- do.call(rbind, lapply(Sys.glob(file.path(d, "seed_*/evaluation/detection_by_electrode_conventional.csv")), read.csv))
  x <- rbind(b[, c("method","channel","rate")], g[, c("method","channel","rate")])
  a <- aggregate(rate ~ method + channel, x, mean); a$n <- as.integer(aggregate(rate ~ method + channel, x, length)$rate); a$scenario <- c
  out[[c]] <- a
}
r <- do.call(rbind, out); write.csv(r, file.path(out_dir, "electrode_fpr.csv"), row.names = FALSE); cat("rows", nrow(r), "n range", range(r$n), "\n")
