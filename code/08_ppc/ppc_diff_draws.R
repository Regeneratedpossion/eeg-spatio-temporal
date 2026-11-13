# Per-draw difference waveforms of the PPC for one dataset (same 200 draws as code/08_ppc/ppc_one_dataset.R).
# Usage: Rscript ppc_diff_draws.R <scenario> <seed> <out_dir>
suppressPackageStartupMessages({library(brms); library(posterior)})

args <- commandArgs(trailingOnly = TRUE)
cond <- as.integer(args[[1]]); seed <- as.integer(args[[2]]); out_dir <- args[[3]]
B <- 200L
root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = "."))
tag <- sprintf("%03d", seed)
ppc_dir <- file.path(root, "results", "ppc_checks", sprintf("scenario%d_seed_%s", cond, tag))

fit <- readRDS(file.path(root, "results", paste0("scenario", cond), paste0("seed_", tag), "fit.rds"))
sim <- readRDS(file.path(root, "data", paste0("scenario", cond), paste0("sim_", tag, ".rds")))
gt <- readRDS(file.path(root, "data", paste0("scenario", cond), paste0("gt_", tag, ".rds")))
d <- fit$data

# Electrode labels are not in fit$data; recover them from the coordinates.
smap <- unique(sim[, c("channel", "xproj", "yproj")])
key <- function(x, y) sprintf("%.8f|%.8f", as.numeric(x), as.numeric(y))
d$channel <- as.character(smap$channel)[match(key(d$xproj, d$yproj), key(smap$xproj, smap$yproj))]

nd <- d[as.character(d$channel) %in% as.character(gt$active_electrodes), , drop = FALSE]
nd <- nd[order(as.character(nd$participant), as.character(nd$condition), as.character(nd$channel), nd$time), , drop = FALSE]
stopifnot(nrow(nd) == 18360L)

# Same draws as the PPC script.
set.seed(20261110)
ids <- sort(sample(seq_len(posterior::ndraws(fit)), B, replace = FALSE))
info <- readLines(file.path(ppc_dir, "run_info.txt"))
saved_ids <- as.integer(strsplit(sub("posterior_draw_ids=", "", info[grepl("^posterior_draw_ids=", info)]), ",")[[1]])
stopifnot(identical(ids, saved_ids))

yrep <- brms::posterior_predict(fit, newdata = nd, re_formula = NULL, draw_ids = ids)

# W maps the rows to the average difference wave (condition 2 minus condition 1).
times <- sort(unique(nd$time))
n_avg <- nlevels(factor(nd$participant)) * length(unique(nd$channel))
W <- matrix(0, nrow(nd), length(times))
W[cbind(seq_len(nrow(nd)), match(nd$time, times))] <- ifelse(as.character(nd$condition) == "cond2", 1, -1) / n_avg
obs <- as.vector(nd$eeg_mean %*% W)
rep <- yrep %*% W

# Check against the saved summary of the PPC script.
saved <- read.csv(file.path(ppc_dir, "wave_diff_active.csv"))
cat("max |obs - saved obs|      :", max(abs(obs - saved$obs)), "\n")
cat("max |rep median - saved q50|:", max(abs(apply(rep, 2, median) - saved$rep_q50)), "\n")
win <- times >= 0.160 - 1e-9 & times <= 0.340 + 1e-9
t_rep <- rowMeans(rep[, win, drop = FALSE]); t_obs <- mean(obs[win])
sdf <- read.csv(file.path(ppc_dir, "ppc_summary_diff.csv"))
sdf <- sdf[sdf$statistic == "mean_160_340_uV", ]
cat("T_mean obs", t_obs, "saved", sdf$observed, "| rep median", median(t_rep), "saved", sdf$rep_median, "\n")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(data.frame(time = times, obs = obs), file.path(out_dir, sprintf("ppc_obs_scenario%d_seed%s.csv", cond, tag)), row.names = FALSE)
write.csv(data.frame(draw = rep(seq_len(B), times = length(times)), time = rep(times, each = B), value = as.vector(rep)),
          file.path(out_dir, sprintf("ppc_rep_scenario%d_seed%s.csv", cond, tag)), row.names = FALSE)
cat("saved to", out_dir, "\n")
