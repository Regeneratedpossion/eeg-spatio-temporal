#!/usr/bin/env Rscript
# Model checks of one existing converged BGAM fit: posterior predictive checks and the residual-source diagnostic.
# The waveform is the difference wave (condition 2 minus condition 1), mean over 20 participants and the 9 active electrodes.
# Peak amplitude and latency are the maximum of the waveform between 160 and 340 ms.
# Usage: Rscript code/08_ppc/ppc_one_dataset.R <scenario:1-4> <seed:1-100> [n_draws=200]

suppressPackageStartupMessages({library(brms); library(posterior)})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("Usage: Rscript ppc_one_dataset.R <scenario:1-4> <seed:1-100> [n_draws=200]")
cond <- as.integer(args[[1]])
seed <- as.integer(args[[2]])
B <- if (length(args) >= 3L) as.integer(args[[3]]) else 200L
stopifnot(cond %in% 1:4, seed %in% 1:100, is.finite(B), B >= 20)

root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = getwd()), mustWork = TRUE)
tag <- sprintf("%03d", seed)
# The fit is read from results/scenario<N>/. Set PPC_FIT_ROOT to read fits from another results folder.
fit_root <- Sys.getenv("PPC_FIT_ROOT", unset = file.path(root, "results"))

fit_path <- file.path(fit_root, paste0("scenario", cond), paste0("seed_", tag), "fit.rds")
sim_path <- file.path(root, "data", paste0("scenario", cond), paste0("sim_", tag, ".rds"))
gt_path  <- file.path(root, "data", paste0("scenario", cond), paste0("gt_", tag, ".rds"))
paths <- c(fit_path, sim_path, gt_path)
if (!all(file.exists(paths))) stop("Missing input file(s):\n", paste(paths[!file.exists(paths)], collapse = "\n"))

out_dir <- file.path(Sys.getenv("PPC_OUT_ROOT", file.path(root, "results", "ppc_checks")), sprintf("scenario%d_seed_%s", cond, tag))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Load existing model and data ----
fit <- readRDS(fit_path)
sim <- readRDS(sim_path)
gt <- readRDS(gt_path)
d <- fit$data

required <- c("participant", "condition", "time", "xproj", "yproj", "eeg_sem", "eeg_mean")
if (!all(required %in% names(d))) stop("fit$data is missing columns: ", paste(setdiff(required, names(d)), collapse = ", "))

if (!"channel" %in% names(d)) {
  smap <- unique(sim[, c("channel", "xproj", "yproj")])
  key <- function(x, y) sprintf("%.8f|%.8f", as.numeric(x), as.numeric(y))
  smap_key <- key(smap$xproj, smap$yproj)
  if (anyDuplicated(smap_key)) stop("Spatial coordinates are not unique across electrode names")
  d$channel <- as.character(smap$channel)[match(key(d$xproj, d$yproj), smap_key)]
}
if (anyNA(d$channel)) stop("Could not align channels using xproj/yproj")

conditions <- c("cond1", "cond2")
if (!all(conditions %in% as.character(unique(d$condition)))) stop("Expected condition levels cond1 and cond2")
if (length(gt$active_electrodes) != 9) stop("Expected 9 active electrodes")

# ---- 2. Active electrodes, both conditions ----
keep <- as.character(d$channel) %in% as.character(gt$active_electrodes)
nd <- d[keep, , drop = FALSE]
nd <- nd[order(as.character(nd$participant), as.character(nd$condition), as.character(nd$channel), nd$time), , drop = FALSE]
if (nrow(nd) != 20L * 2L * 9L * 51L) stop("Expected 18360 active-electrode observations; got ", nrow(nd))

series_id <- interaction(nd$participant, nd$condition, nd$channel, drop = TRUE)
series <- split(seq_len(nrow(nd)), series_id)
if (length(series) != 20L * 2L * 9L || !all(lengths(series) == 51)) stop("Incomplete participant/condition/channel series")
if (!all(vapply(series, function(ix) all(diff(nd$time[ix]) > 0), logical(1)))) stop("Timepoints not increasing within series")

# ---- 3. Posterior predictive replications ----
set.seed(20261110)
if (B > posterior::ndraws(fit)) stop("Requested draws exceed total posterior draws")
ids <- sort(sample(seq_len(posterior::ndraws(fit)), B, replace = FALSE))
message("Generating ", B, " aligned posterior predictions for ", nrow(nd), " observations...")
mu <- brms::posterior_epred(fit, newdata = nd, re_formula = NULL, draw_ids = ids)
yrep <- brms::posterior_predict(fit, newdata = nd, re_formula = NULL, draw_ids = ids)
stopifnot(identical(dim(mu), dim(yrep)), nrow(mu) == B, ncol(mu) == nrow(nd),
          all(is.finite(mu)), all(is.finite(yrep)))

# ---- 4. PPC 1: waveform summary statistics ----
sel2 <- as.character(nd$condition) == "cond2"
times <- sort(unique(nd$time[sel2]))
if (length(times) != 51) stop("Expected exactly 51 timepoints")
win <- times >= 0.160 - 1e-9 & times <= 0.340 + 1e-9

add_stat <- function(n, observed, replicated) {
  q05 <- unname(quantile(replicated, 0.05)); q95 <- unname(quantile(replicated, 0.95))
  data.frame(statistic = n, observed = observed, rep_median = median(replicated),
             rep_q05 = q05, rep_q95 = q95, below_90 = observed < q05, above_90 = observed > q95)
}

# Difference wave: condition 2 minus condition 1.
# Pair each condition-2 row with the condition-1 row of the same participant, electrode and time point.
row_key <- paste(nd$participant, nd$channel, nd$time, sep = "|")
i2 <- which(as.character(nd$condition) == "cond2")
i1 <- match(row_key[i2], row_key[as.character(nd$condition) == "cond1"])
i1 <- which(as.character(nd$condition) == "cond1")[i1]
if (anyNA(i1)) stop("Could not pair condition-2 rows with condition-1 rows")
obs_diff <- nd$eeg_mean[i2] - nd$eeg_mean[i1]
rep_diff <- yrep[, i2, drop = FALSE] - yrep[, i1, drop = FALSE]
time_diff <- nd$time[i2]
wave_d <- do.call(rbind, lapply(times, function(tt) {
  ix <- which(abs(time_diff - tt) < 1e-9)
  if (length(ix) != 20L * 9L) stop("Expected 180 pairs per time point")
  reps <- rowMeans(rep_diff[, ix, drop = FALSE])
  data.frame(time = tt, obs = mean(obs_diff[ix]), rep_mean = mean(reps),
             rep_q05 = unname(quantile(reps, 0.05)), rep_q50 = median(reps), rep_q95 = unname(quantile(reps, 0.95)))
}))
write.csv(wave_d, file.path(out_dir, "wave_diff_active.csv"), row.names = FALSE)
rep_wave_d <- sapply(times, function(tt) rowMeans(rep_diff[, abs(time_diff - tt) < 1e-9, drop = FALSE]))
peak_ix_d <- which(win)[which.max(wave_d$obs[win])]
rep_peak_ix_d <- apply(rep_wave_d[, win, drop = FALSE], 1, which.max)
stat_table <- list(
  add_stat("peak_amplitude_uV", wave_d$obs[peak_ix_d], apply(rep_wave_d[, win, drop = FALSE], 1, max)),
  add_stat("peak_latency_ms", 1000 * times[peak_ix_d], 1000 * times[which(win)[rep_peak_ix_d]]),
  add_stat("mean_160_340_uV", mean(wave_d$obs[win]), rowMeans(rep_wave_d[, win, drop = FALSE])))

# ---- 5. PPC 2: residual SD and autocorrelation ----
resid_summary <- function(e, indices) {
  ss <- lapply(indices, function(ix) {
    v <- e[ix]
    c(sd = sd(v), lag1 = cor(v[-length(v)], v[-1]),
      lag2 = cor(v[-c(length(v) - 1, length(v))], v[-c(1, 2)]))
  })
  colMeans(do.call(rbind, ss), na.rm = TRUE)
}
noise_rows <- vector("list", 4L * B)
k <- 0L
sel_by_cond <- lapply(conditions, function(cn)
  series[vapply(series, function(ix) as.character(nd$condition[ix[1]]) == cn, logical(1))])
names(sel_by_cond) <- conditions
for (j in seq_len(B)) {
  eobs <- nd$eeg_mean - mu[j, ]
  erep <- yrep[j, ] - mu[j, ]
  for (condname in conditions) {
    selected <- sel_by_cond[[condname]]
    o <- resid_summary(eobs, selected)
    r <- resid_summary(erep, selected)
    k <- k + 1L
    noise_rows[[k]] <- data.frame(draw = j, condition = condname, type = "observed", sd = o["sd"], lag1 = o["lag1"], lag2 = o["lag2"])
    k <- k + 1L
    noise_rows[[k]] <- data.frame(draw = j, condition = condname, type = "replicated", sd = r["sd"], lag1 = r["lag1"], lag2 = r["lag2"])
  }
}
noise <- do.call(rbind, noise_rows)
write.csv(noise, file.path(out_dir, "noise_ppc.csv"), row.names = FALSE)

# ---- 6. Summary ----
for (condname in conditions) {
  for (st in c("sd", "lag1", "lag2")) {
    vals_obs <- noise[noise$condition == condname & noise$type == "observed", st]
    vals_rep <- noise[noise$condition == condname & noise$type == "replicated", st]
    stat_table[[length(stat_table) + 1L]] <- add_stat(paste0(condname, "_", st), median(vals_obs), vals_rep)
  }
}
summary <- do.call(rbind, stat_table)
write.csv(summary, file.path(out_dir, "ppc_summary.csv"), row.names = FALSE)

writeLines(c(paste0("scenario=", cond), paste0("seed=", seed), paste0("draws=", B),
             paste0("posterior_draw_ids=", paste(ids, collapse = ",")), paste0("fit=", fit_path)),
           file.path(out_dir, "run_info.txt"))

# ---- 7. Residual source: observed residuals against the true noise (simulation only) ----
# This part uses its own draws (seed 20261110 + 1000 * scenario + seed) and its own row order.
rm(mu, yrep); invisible(gc())
residual_source <- function(fit, gt, sim, root, cond, seed, B) {
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
  result <- do.call(rbind, output)

  result
}
result <- residual_source(fit, gt, sim, root, cond, seed, B)
rs_dir <- file.path(dirname(out_dir), "residual_source")
dir.create(rs_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(result, file.path(rs_dir, sprintf("scenario%d_seed_%03d.csv", cond, seed)), row.names = FALSE)

message("Saved PPC files in: ", out_dir)
print(summary, row.names = FALSE)
print(result)
