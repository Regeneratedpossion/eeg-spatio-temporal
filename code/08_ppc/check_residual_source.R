# Usage: Rscript --vanilla code/08_ppc/check_residual_source.R <scenario:1-4> <seed> [n_draws=200]
suppressPackageStartupMessages({library(brms); library(dplyr)})

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 2L)
cond <- as.integer(args[1])
seed <- as.integer(args[2])
B <- if (length(args) >= 3L) as.integer(args[3]) else 200L
stopifnot(cond %in% 1:4, seed %in% 1:100, B >= 1)

root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = "."))
scen <- paste0("scenario", cond)
paths <- list(truth_file = file.path(root, "data", scen, sprintf("gt_%03d.rds", seed)), sim_file = file.path(root, "data", scen, sprintf("sim_%03d.rds", seed)))
# The fit is read from results/scenario<N>/. Set PPC_FIT_ROOT to read fits from another results folder.
fit_root <- Sys.getenv("PPC_FIT_ROOT", unset = file.path(root, "results"))
fit_path <- file.path(fit_root, scen, sprintf("seed_%03d", seed), "fit.rds")

fit <- readRDS(fit_path)
gt <- readRDS(paths$truth_file)
sim <- readRDS(paths$sim_file)

# 1. Exact ERP template used by the DGP
template <- new.env(parent = globalenv())
sys.source(file.path(root, "resources/brms_meeg/code/erp_template.R"), envir = template)

signal_at_onset <- function(onset_ms) {
  k <- match(onset_ms, template$Xf)
  if (is.na(k)) stop("Invalid participant onset")
  c(rep(0, k - 2L), template$erp, rep(0, template$Nf - k - length(template$erp) + 2L))
}

# 2. Recover electrode labels omitted from brms data
nd <- fit$data
sensors <- unique(sim[c("channel", "xproj", "yproj")])
key <- function(x, y) sprintf("%.8f|%.8f", x, y)
stopifnot(!anyDuplicated(key(sensors$xproj, sensors$yproj)))
nd$channel <- sensors$channel[match(key(nd$xproj, nd$yproj), key(sensors$xproj, sensors$yproj))]
stopifnot(!anyNA(nd$channel))

nd <- nd[nd$channel %in% gt$active_electrodes, , drop = FALSE]
nd <- nd[order(nd$condition, nd$participant, nd$channel, nd$time), , drop = FALSE]
stopifnot(nrow(nd) == 20L * 2L * 9L * 51L)

# 3. Participant-specific true mean: amplitude x spatial weight x temporal signal (cond2 only)
pp <- gt$participant_parameters
pidx <- match(as.character(nd$participant), pp$participant)
stopifnot(!anyNA(pidx))
spatial <- unique(gt$continuous_effect[c("channel", "spatial_weight")])
stopifnot(!anyDuplicated(spatial$channel))
widx <- match(nd$channel, spatial$channel)
stopifnot(!anyNA(widx))

onsets <- unique(pp$onset_ms)
signals <- setNames(lapply(onsets, signal_at_onset), as.character(onsets))
time_idx <- match(as.integer(round(nd$time * 1000)), template$Xf)
stopifnot(!anyNA(time_idx))
signal_value <- vapply(seq_len(nrow(nd)), function(i)
  signals[[as.character(pp$onset_ms[pidx[i]])]][time_idx[i]], numeric(1))

true_mean <- ifelse(as.character(nd$condition) == "cond2",
                    pp$amplitude[pidx] * spatial$spatial_weight[widx] * signal_value, 0)

# 4. True noise residual
y <- nd$eeg_mean
e_gt <- y - true_mean

# 5. BGAM predictions with participant effects
set.seed(20261110 + 1000L * cond + seed)
n_draws <- posterior::ndraws(fit)
stopifnot(B <= n_draws)
ids <- sample.int(n_draws, B)
mu <- brms::posterior_epred(fit, newdata = nd, re_formula = NULL, draw_ids = ids)
stopifnot(is.matrix(mu), nrow(mu) == B, ncol(mu) == nrow(nd))

# 6. Same series and lag statistic for both residuals
series_id <- interaction(nd$participant, nd$condition, nd$channel, drop = TRUE)
groups <- split(seq_len(nrow(nd)), series_id)
lag_mean <- function(e, index_list, lag) {
  mean(vapply(index_list, function(ii) {
    x <- e[ii]; n <- length(x)
    cor(x[seq_len(n - lag)], x[(lag + 1L):n])
  }, numeric(1)), na.rm = TRUE)
}

output <- list()
for (condition_name in c("cond1", "cond2")) {
  indices <- which(as.character(nd$condition) == condition_name)
  local_groups <- lapply(groups, function(ii) intersect(ii, indices))
  local_groups <- Filter(function(ii) length(ii) == 51L, local_groups)
  stopifnot(length(local_groups) == 20L * 9L)

  gt_lag1 <- lag_mean(e_gt, local_groups, 1L)
  gt_lag2 <- lag_mean(e_gt, local_groups, 2L)
  bg_lag1 <- bg_lag2 <- numeric(B)
  for (j in seq_len(B)) {
    e_bg <- y - mu[j, ]
    bg_lag1[j] <- lag_mean(e_bg, local_groups, 1L)
    bg_lag2[j] <- lag_mean(e_bg, local_groups, 2L)
  }
  bg_mean <- colMeans(mu)

  output[[condition_name]] <- data.frame(
    scenario = cond, seed = seed, condition = condition_name,
    gt_noise_lag1 = gt_lag1, bgam_lag1 = median(bg_lag1),
    gt_noise_lag2 = gt_lag2, bgam_lag2 = median(bg_lag2),
    lag1_difference = median(bg_lag1) - gt_lag1,
    lag2_difference = median(bg_lag2) - gt_lag2,
    mean_rmse = sqrt(mean((bg_mean[indices] - true_mean[indices])^2)))
}
result <- bind_rows(output)

out_dir <- file.path(Sys.getenv("PPC_OUT_ROOT", file.path(root, "results", "ppc_checks")), "residual_source")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(result, file.path(out_dir, sprintf("scenario%d_seed_%03d.csv", cond, seed)), row.names = FALSE)
print(result)
