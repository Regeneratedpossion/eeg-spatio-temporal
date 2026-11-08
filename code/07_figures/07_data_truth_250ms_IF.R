# Recovered from a Claude Code session on 2026-10-09 (session 1840aaa0).
# Original: a scratch script (truth.R). Usage: PROJECT_ROOT=<project> Rscript make_truth_250ms_IF.R
suppressPackageStartupMessages(library(dplyr))
root <- Sys.getenv("PROJECT_ROOT"); stopifnot(nzchar(root))   # the simulation project (on the cluster or the local copy)
source(file.path(root, "code", "R", "superpopulation_truth.R"))
template <- file.path(root, "resources", "brms_meeg", "code", "erp_template.R")
adj <- readRDS(file.path(root, "resources", "adjacency", "adjacency_xproj_yproj_61ch.rds"))
out <- list()
for (c in c("I", "F")) {
  gt <- readRDS(sprintf("%s/data/single_node_cond%s_seeds001_100/gt_001.rds", root, c))
  tr <- paired_radius_superpopulation_truth(gt, template)
  ch <- rownames(adj); act <- ch %in% gt$active_electrodes
  ring <- ifelse(act, 0L, NA_integer_)
  for (r in 1:2) { prev <- which(ring == r - 1L); nb <- which(colSums(adj[prev, , drop = FALSE]) > 0 & is.na(ring)); ring[nb] <- r }
  ring[is.na(ring)] <- 3L
  t250 <- tr[abs(tr$time - 0.25) < 1e-9, c("channel", "effect")]
  t250$ring <- ring[match(t250$channel, ch)]
  t250$scenario <- c
  out[[c]] <- t250
}
write.csv(do.call(rbind, out), "truth_250ms_IF.csv", row.names = FALSE)
cat("ok\n")
