#!/usr/bin/env Rscript
# PPC on one existing converged BGAM fit.
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

# ---- 4. PPC 1: condition 2 average waveform ----
sel2 <- as.character(nd$condition) == "cond2"
times <- sort(unique(nd$time[sel2]))
if (length(times) != 51) stop("Expected exactly 51 timepoints")

wave <- do.call(rbind, lapply(times, function(tt) {
  ix <- which(sel2 & abs(nd$time - tt) < 1e-9)
  if (length(ix) != 20L * 9L) stop("Expected 180 observations per condition-2 timepoint")
  reps <- rowMeans(yrep[, ix, drop = FALSE])
  data.frame(time = tt, obs = mean(nd$eeg_mean[ix]), rep_mean = mean(reps),
             rep_q05 = unname(quantile(reps, 0.05)), rep_q50 = median(reps),
             rep_q95 = unname(quantile(reps, 0.95)))
}))
write.csv(wave, file.path(out_dir, "wave_cond2_active.csv"), row.names = FALSE)

# ---- 5. Waveform summary statistics ----
win <- times >= 0.160 - 1e-9 & times <= 0.340 + 1e-9
late <- times >= 0.270 - 1e-9 & times <= 0.340 + 1e-9
rep_wave <- sapply(times, function(tt) rowMeans(yrep[, sel2 & abs(nd$time - tt) < 1e-9, drop = FALSE]))
if (nrow(rep_wave) != B || ncol(rep_wave) != 51) stop("Bad wave shape")

stat_table <- list()
add_stat <- function(n, observed, replicated) {
  q05 <- unname(quantile(replicated, 0.05)); q95 <- unname(quantile(replicated, 0.95))
  data.frame(statistic = n, observed = observed, rep_median = median(replicated),
             rep_q05 = q05, rep_q95 = q95, below_90 = observed < q05, above_90 = observed > q95)
}
peak_ix <- which(win)[which.max(wave$obs[win])]
rep_peak_ix <- apply(rep_wave[, win, drop = FALSE], 1, which.max)
stat_table[[1]] <- add_stat("peak_amplitude_uV", wave$obs[peak_ix], apply(rep_wave[, win, drop = FALSE], 1, max))
stat_table[[2]] <- add_stat("peak_latency_ms", 1000 * times[peak_ix], 1000 * times[which(win)[rep_peak_ix]])
stat_table[[3]] <- add_stat("mean_160_340_uV", mean(wave$obs[win]), rowMeans(rep_wave[, win, drop = FALSE]))
stat_table[[4]] <- add_stat("mean_270_340_uV", mean(wave$obs[late]), rowMeans(rep_wave[, late, drop = FALSE]))

# ---- 6. PPC 2: residual SD and autocorrelation ----
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

# ---- 7. Summary ----
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

# ---- 8. Figure (pdf: compute nodes have no X11 display for png) ----
pdf(file.path(out_dir, "wave_cond2_active.pdf"), width = 9.5, height = 6)
plot(wave$time * 1000, wave$obs, type = "n", ylim = range(wave[, c("obs", "rep_q05", "rep_q95")]),
     xlab = "Time (ms)", ylab = "Mean EEG (uV)",
     main = sprintf("PPC: Condition 2, active electrodes | scenario %s seed %d", cond, seed))
polygon(c(wave$time, rev(wave$time)) * 1000, c(wave$rep_q05, rev(wave$rep_q95)),
        border = NA, col = adjustcolor("steelblue", alpha.f = 0.25))
lines(wave$time * 1000, wave$rep_mean, col = "steelblue", lwd = 2)
lines(wave$time * 1000, wave$obs, col = "darkorange", lwd = 2)
abline(v = c(160, 340), lty = 3)
legend("topright", legend = c("Observed", "Posterior replication mean", "90% replication interval"),
       col = c("darkorange", "steelblue", adjustcolor("steelblue", alpha.f = 0.25)),
       lwd = c(2, 2, 8), bty = "n")
dev.off()

message("Saved PPC files in: ", out_dir)
print(summary, row.names = FALSE)
