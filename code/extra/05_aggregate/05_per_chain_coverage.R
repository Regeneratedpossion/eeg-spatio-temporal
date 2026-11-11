# Per fit: onset/offset of each chain (main rule) and the 95% interval of g at every cell.
# Usage: Rscript code/extra/05_aggregate/05_per_chain_coverage.R <condition> <seed>
suppressPackageStartupMessages(library(brms))
project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = here::here()))
source(file.path(project_root, "code", "04_evaluate", "04_evaluate_bgam.R"))

a <- commandArgs(trailingOnly = TRUE)
cond <- a[[1]]; seed <- as.integer(a[[2]]); tag <- sprintf("%03d", seed)
seed_dir <- file.path(project_root, "results", sprintf("single_node_cond%s_ds_cr_seeds001_100", cond), sprintf("seed_%s", tag))
data_dir <- file.path(project_root, "data", sprintf("single_node_cond%s_seeds001_100", cond))
out_dir <- file.path(seed_dir, "evaluation")
f_chain <- file.path(out_dir, "per_chain_onset.csv"); f_cov <- file.path(out_dir, "coverage_cells.csv")
if (file.exists(f_chain) && file.exists(f_cov)) quit(save = "no")

s <- evaluation_settings()
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

# Same SESOI as the main analysis.
sesoi <- read.csv(file.path(out_dir, "sesoi.csv"))
sesoi <- sesoi[sesoi$rule == "global", ]
delta <- sesoi$delta[match(grid$channel, sesoi$channel)]
n_chains <- posterior::nchains(fit)
chain <- rep(seq_len(n_chains), each = nrow(draws) / n_chains)
to_matrix <- function(x) matrix(x, nrow = length(times), ncol = length(channels))

rows <- lapply(seq_len(n_chains), function(k) {
    dk <- draws[chain == k, , drop = FALSE]
    p <- colMeans(sweep(dk, 2, delta, ">"))
    odds <- pmax(p, 1 / nrow(dk)) / pmax(1 - p, 1 / nrow(dk))
    oo <- extract_onset_offset(to_matrix(odds > s$main_odds), times, adj, s)
    data.frame(cond = cond, seed = seed, chain = k, onset = oo[["onset"]], offset = oo[["offset"]])
})
write.csv(do.call(rbind, rows), f_chain, row.names = FALSE)

template <- file.path(project_root, "resources", "brms_meeg", "code", "erp_template.R")
truth_df <- paired_radius_superpopulation_truth(gt, template)
truth <- truth_df$effect[match(paste(grid$channel, grid$time), paste(truth_df$channel, truth_df$time))]
lo <- apply(draws, 2, quantile, probs = 0.025); hi <- apply(draws, 2, quantile, probs = 0.975)
write.csv(data.frame(channel = grid$channel, time = grid$time, mean_diff = colMeans(draws), lo = lo, hi = hi,
                     truth = truth, covers = truth >= lo & truth <= hi), f_cov, row.names = FALSE)
