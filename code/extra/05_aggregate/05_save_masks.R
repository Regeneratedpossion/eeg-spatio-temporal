# Saves the main-setting detection masks of one seed: BGAM (odds 10) or the conventional methods (alpha 0.05).
# Usage: Rscript code/extra/05_aggregate/05_save_masks.R <condition> <seed> <bgam|conv>
suppressPackageStartupMessages(library(dplyr))
project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = here::here()))
source(file.path(project_root, "code", "03_conventional", "03_evaluate_conventional.R"))

a <- commandArgs(trailingOnly = TRUE)
cond <- a[[1]]; seed <- as.integer(a[[2]]); what <- a[[3]]
tag <- sprintf("%03d", seed)
seed_dir <- file.path(project_root, "results", sprintf("single_node_cond%s_ds_cr_seeds001_100", cond),
                      sprintf("seed_%s", tag))
data_dir <- file.path(project_root, "data", sprintf("single_node_cond%s_seeds001_100", cond))
out_file <- file.path(seed_dir, "evaluation", sprintf("masks_%s.rds", what))
if (file.exists(out_file)) quit(save = "no")

s <- evaluation_settings()
sim <- readRDS(file.path(data_dir, sprintf("sim_%s.rds", tag)))
gt <- readRDS(file.path(data_dir, sprintf("gt_%s.rds", tag)))
channels <- sort(unique(sim$channel)); times <- sort(unique(sim$time))
adj <- readRDS(file.path(project_root, "resources", "adjacency", "adjacency_xproj_yproj_61ch.rds"))
adj <- adj[match(channels, rownames(adj)), match(channels, rownames(adj))]
active_ch <- channels %in% gt$active_electrodes
truth <- outer(times >= s$true_onset - 1e-9 & times <= s$true_offset + 1e-9, active_ch, "&")

if (what == "bgam") {
  fit <- readRDS(file.path(seed_dir, "fit.rds"))
  sensors <- unique(sim[, c("channel", "xproj", "yproj")]); sensors <- sensors[order(sensors$channel), ]
  times <- sort(unique(fit$data$time))
  grid <- expand.grid(time = times, channel = channels, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  grid <- merge(grid, sensors, by = "channel", sort = FALSE)
  grid <- grid[order(match(grid$channel, channels), grid$time), ]
  grid$eeg_sem <- 1
  levs <- levels(fit$data$condition)
  nd1 <- transform(grid, condition = factor(levs[1], levels = levs))
  nd2 <- transform(grid, condition = factor(levs[2], levels = levs))
  draws <- brms::posterior_epred(fit, newdata = nd2, re_formula = NA) -
    brms::posterior_epred(fit, newdata = nd1, re_formula = NA)
  rm(fit); est <- colMeans(draws)
  base <- grid$time >= 0 & grid$time <= s$main_baseline_end + 1e-9
  delta <- max(0, max(est[base]))
  p <- colMeans(sweep(draws, 2, rep(delta, ncol(draws)), ">"))
  odds <- pmax(p, 1 / nrow(draws)) / pmax(1 - p, 1 / nrow(draws))
  masks <- list(BGAM = matrix(odds > s$main_odds, nrow = length(times), ncol = length(channels)))
} else {
  cs <- conventional_settings()
  arr <- difference_array(sim, times, channels)
  cell_p <- apply(arr, c(2, 3), function(x) t.test(x)$p.value)
  p_list <- list(
    uncorrected = cell_p,
    BH = matrix(p.adjust(cell_p, "BH"), nrow(cell_p)),
    BY = matrix(p.adjust(cell_p, "BY"), nrow(cell_p)),
    Holm = matrix(p.adjust(cell_p, "holm"), nrow(cell_p)),
    cluster_mass = mne_cell_pvalues(arr, adj, cs, tfce = FALSE),
    TFCE = mne_cell_pvalues(arr, adj, cs, tfce = TRUE)
  )
  masks <- lapply(p_list, function(p) p < cs$main_alpha)
}
saveRDS(list(masks = masks, truth = truth, times = times, channels = channels, active = active_ch,
             adj = adj), out_file)
cat("saved", out_file, "\n")
