# Per fit: posterior mean, SD and 90% interval of g at every cell, with the true g and the electrode ring.
# g is the population-level condition effect used for the posterior odds (no participant terms, delta not subtracted).
# Usage: Rscript code/extra/04_coverage_z.R <condition> <seed>
# Writes results/scenario<N>/seed_<NNN>/evaluation/coverage_z.csv
suppressPackageStartupMessages(library(brms))
project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = here::here()))
source(file.path(project_root, "code", "04_evaluate", "04_evaluate_bgam.R"))

a <- commandArgs(trailingOnly = TRUE)
cond <- a[[1]]; seed <- as.integer(a[[2]]); tag <- sprintf("%03d", seed)
seed_dir <- file.path(project_root, "results", sprintf("scenario%s", cond), sprintf("seed_%s", tag))
data_dir <- file.path(project_root, "data", sprintf("scenario%s", cond))
out_file <- file.path(seed_dir, "evaluation", "coverage_z.csv")
if (file.exists(out_file)) quit(save = "no")

fit <- readRDS(file.path(seed_dir, "fit.rds"))
gt <- readRDS(file.path(data_dir, sprintf("gt_%s.rds", tag)))
sim <- readRDS(file.path(data_dir, sprintf("sim_%s.rds", tag)))
sensors <- unique(sim[, c("channel", "xproj", "yproj")]); sensors <- sensors[order(sensors$channel), ]
channels <- sensors$channel
times <- sort(unique(fit$data$time))
adj <- readRDS(file.path(project_root, "resources", "adjacency", "adjacency_xproj_yproj_61ch.rds"))
adj <- adj[match(channels, rownames(adj)), match(channels, rownames(adj))]

grid <- expand.grid(time = times, channel = channels, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
grid <- merge(grid, sensors, by = "channel", sort = FALSE)
grid <- grid[order(match(grid$channel, channels), grid$time), ]
grid$eeg_sem <- 1
levs <- levels(fit$data$condition)
nd1 <- transform(grid, condition = factor(levs[1], levels = levs))
nd2 <- transform(grid, condition = factor(levs[2], levels = levs))
draws <- posterior_epred(fit, newdata = nd2, re_formula = NA) - posterior_epred(fit, newdata = nd1, re_formula = NA)

template <- file.path(project_root, "resources", "brms_meeg", "code", "erp_template.R")
truth_df <- paired_radius_superpopulation_truth(gt, template)
truth <- truth_df$effect[match(paste(grid$channel, grid$time), paste(truth_df$channel, truth_df$time))]
ring <- electrode_rings(adj, channels, channels %in% gt$active_electrodes)[match(grid$channel, channels)]

q <- apply(draws, 2, quantile, probs = c(0.05, 0.95))
write.csv(data.frame(channel = grid$channel, time = grid$time, ring = ring, truth = truth,
                     mean = colMeans(draws), sd = apply(draws, 2, sd), lo = q[1, ], hi = q[2, ]),
          out_file, row.names = FALSE)
