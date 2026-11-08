# Recovered from a Claude Code session on 2026-10-16 (session 4068e1d2).
# Original command line: cat > truth_all.R
suppressMessages(library(dplyr))
root <- Sys.getenv("PROJECT_ROOT"); stopifnot(nzchar(root))   # Usage: PROJECT_ROOT=<project> Rscript make_truth_250.R
source(file.path(root, "code/R/superpopulation_truth.R"))
tmpl <- file.path(root, "resources/brms_meeg/code/erp_template.R")
sc <- c(I = 1, H = 2, F = 3, E = 4)
tr <- bind_rows(lapply(names(sc), function(c) {
  gt <- readRDS(file.path(root, sprintf("data/single_node_cond%s_seeds001_100/gt_001.rds", c)))
  t <- paired_radius_superpopulation_truth(gt, tmpl); t$scenario <- sc[[c]]; t
})) |> filter(abs(time - 0.25) < 1e-9)
write.csv(tr[, c("scenario", "channel", "time", "effect")], "truth_250.csv", row.names = FALSE)
