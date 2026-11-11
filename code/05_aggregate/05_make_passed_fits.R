# Writes resources/passed_fits.txt: the fits that pass the convergence criterion of Section 2.2.5.
# A fit passes if, over all cells of g, max R-hat < 1.01 and min bulk ESS and min tail ESS > 400.
# Reads results/single_node_cond*_ds_cr_seeds001_100/seed_*/evaluation/evaluation_main.csv. Run from the project root.
source(here::here("code/R/aggregate_utils.R"))
root <- here::here()
rows <- do.call(rbind, lapply(c("1", "2", "3", "4"), function(cond) {
  files <- Sys.glob(file.path(scenario_results_dir(file.path(root, "results"), cond), "seed_*", "evaluation", "evaluation_main.csv"))
  do.call(rbind, lapply(files, function(f) {
    x <- read.csv(f)[1, c("max_rhat_g", "min_ess_bulk_g", "min_ess_tail_g")]
    data.frame(cond = cond, seed = as.integer(sub(".*seed_([0-9]+).*", "\\1", f)), x)
  }))
}))
rows <- rows[order(match(rows$cond, c("1", "3", "2", "4")), rows$seed), ]
passed <- rows[rows$max_rhat_g < 1.01 & rows$min_ess_bulk_g > 400 & rows$min_ess_tail_g > 400, ]
out <- Sys.getenv("PASSED_OUT", unset = file.path(root, "resources", "passed_fits.txt"))
write.table(passed[, c("cond", "seed")], out, quote = FALSE, row.names = FALSE, col.names = FALSE)
cat("fits per scenario:", table(factor(passed$cond, c("1", "3", "2", "4"))), "of", nrow(rows), "\n")
