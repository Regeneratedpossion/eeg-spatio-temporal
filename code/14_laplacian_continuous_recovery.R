###########################################################
# 14_laplacian_continuous_recovery.R
#
# Continuous recovery diagnostics for Laplacian low-rank fits.
# This evaluates the reconstructed sensor x time posterior effect
# against the continuous Gaussian ground-truth field, instead of
# thresholding sensors into active/inactive labels.
###########################################################

suppressPackageStartupMessages({
    library(brms)
    library(dplyr)
    library(tidyr)
})

source(file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()), "code/utils.R"))

fit_subdir <- Sys.getenv("LAPLACE_FIT_SUBDIR", "laplacian_basis_brms_fixed")
fit_dir <- file.path(project_root(), "results", fit_subdir)
out_dir <- file.path(project_root(), "results", "diagnostics",
                     "laplacian_continuous")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

data_subdir <- Sys.getenv("LAPLACE_DATA_SUBDIR",
                          "montecarlo_p20_gaussian_hardlabel")
seed <- as.integer(Sys.getenv("LAPLACE_SEED", "1"))
adapt_tag <- Sys.getenv("LAPLACE_ADAPT_TAG", "0p999")
k_basis <- as.integer(Sys.getenv("LAPLACE_K", "20"))
time_k <- as.integer(Sys.getenv("LAPLACE_TIME_K", "10"))
iter <- as.integer(Sys.getenv("LAPLACE_ITER", "3000"))
warmup <- as.integer(Sys.getenv("LAPLACE_WARMUP", "1000"))
tag <- sprintf("laplacian_k%d_seed%03d_iter%d_warmup%d_adapt%s",
               k_basis, seed, iter, warmup, adapt_tag)
tag_suffix <- Sys.getenv("LAPLACE_TAG_SUFFIX", "")
if (nzchar(tag_suffix)) {
    tag <- paste0(tag, "_", tag_suffix)
}

summary_file <- file.path(
    fit_dir,
    sprintf("summary_all_iter%d_warmup%d_adapt%s.csv", iter, warmup, adapt_tag)
)
conv_summary <- read.csv(summary_file)
fit_files <- conv_summary$fit_file[order(conv_summary$basis_task)]
if (!all(file.exists(fit_files))) {
    stop("Some Laplacian fit files are missing.")
}

gt <- readRDS(data_path(data_subdir, sprintf("gt_%03d.rds", seed)))
sim <- readRDS(data_path(data_subdir, sprintf("sim_%03d.rds", seed)))
sim_model <- aggregate_trials(sim)
times <- sort(unique(sim_model$time))

cond_levels <- c("cond1", "cond2")

truth_temporal_bump <- function(time, ground_truth) {
    peak_t <- 0.250
    sigma_t <- 0.050
    z <- (time - peak_t) / sigma_t
    val <- exp(-z^2 / 2)
    val_trunc <- exp(-1.5^2 / 2)
    val[abs(z) > 1.5] <- val_trunc
    val_norm <- pmax((val - val_trunc) / (1 - val_trunc), 0)
    val_norm[time < ground_truth$true_onset | time > ground_truth$true_offset] <- 0
    ground_truth$bump_amplitude * val_norm
}

make_continuous_truth <- function(sim_model, ground_truth, times, channels) {
    sensors <- sim_model %>%
        distinct(channel, xproj, yproj) %>%
        filter(channel %in% channels)
    seed_coord <- sensors %>% filter(channel == ground_truth$seed_electrode)
    sensors <- sensors %>%
        mutate(
            spatial_distance = sqrt((xproj - seed_coord$xproj)^2 +
                                    (yproj - seed_coord$yproj)^2),
            spatial_weight_raw = exp(-(spatial_distance^2) /
                                     (2 * ground_truth$sigma_space^2))
        ) %>%
        arrange(match(channel, channels))

    temporal <- truth_temporal_bump(times, ground_truth)
    truth_mat <- outer(temporal, sensors$spatial_weight_raw)
    dimnames(truth_mat) <- list(as.character(times), sensors$channel)
    truth_mat
}

condition_score_draws <- function(obj, times) {
    fit <- obj$fit
    base <- data.frame(
        time = times,
        score = 0,
        score_sem = 1,
        participant = obj$score_data$participant[1]
    )
    nd1 <- base
    nd2 <- base
    nd1$condition <- factor("cond1", levels = cond_levels)
    nd2$condition <- factor("cond2", levels = cond_levels)

    posterior_epred(fit, newdata = nd2, re_formula = NA) -
        posterior_epred(fit, newdata = nd1, re_formula = NA)
}

first_obj <- readRDS(fit_files[1])
channels <- first_obj$channels
n_draws <- nrow(posterior_epred(
    first_obj$fit,
    newdata = data.frame(
        time = times[1],
        score = 0,
        score_sem = 1,
        participant = first_obj$score_data$participant[1],
        condition = factor("cond1", levels = cond_levels)
    ),
    re_formula = NA
))

effect_array <- array(
    0,
    dim = c(n_draws, length(times), length(channels)),
    dimnames = list(NULL, as.character(times), channels)
)

for (file in fit_files) {
    obj <- readRDS(file)
    score_diff <- condition_score_draws(obj, times)
    basis_vector <- obj$basis_vector[channels]
    for (j in seq_along(channels)) {
        effect_array[, , j] <- effect_array[, , j] +
            score_diff * basis_vector[j]
    }
}

estimate_mat <- apply(effect_array, c(2, 3), mean)
truth_mat <- make_continuous_truth(sim_model, gt, times, channels)
stopifnot(identical(dim(estimate_mat), dim(truth_mat)))

truth_vec <- as.vector(truth_mat)
estimate_vec <- as.vector(estimate_mat)
active_idx <- truth_vec > 0
truth_weight_idx <- truth_vec > 0
truth_top25_idx <- truth_vec >= 0.25 * max(truth_vec)
peak_time <- pick_time((gt$true_onset + gt$true_offset) / 2, times)
peak_idx <- match(peak_time, times)

safe_cor <- function(x, y) {
    if (sd(x) == 0 || sd(y) == 0) return(NA_real_)
    cor(x, y)
}

cosine_similarity <- function(x, y) {
    denom <- sqrt(sum(x^2)) * sqrt(sum(y^2))
    if (denom == 0) return(NA_real_)
    sum(x * y) / denom
}

soft_dice <- function(x, y) {
    x <- pmax(x, 0)
    y <- pmax(y, 0)
    if (max(x) > 0) x <- x / max(x)
    if (max(y) > 0) y <- y / max(y)
    denom <- sum(x) + sum(y)
    if (denom == 0) return(NA_real_)
    2 * sum(pmin(x, y)) / denom
}

truth_weighted_mean <- function(x, truth) {
    idx <- truth > 0
    sum(x[idx] * truth[idx]) / sum(truth[idx])
}

summary <- data.frame(
    fit_file = tag,
    dataset = data_subdir,
    k_basis = k_basis,
    iter = iter,
    warmup = warmup,
    adapt_tag = adapt_tag,
    spatiotemporal_cor_all = safe_cor(truth_vec, estimate_vec),
    spatiotemporal_cor_effect_window = safe_cor(truth_vec[active_idx],
                                                estimate_vec[active_idx]),
    spatiotemporal_rmse = sqrt(mean((estimate_vec - truth_vec)^2)),
    spatiotemporal_mae = mean(abs(estimate_vec - truth_vec)),
    spatiotemporal_cosine = cosine_similarity(truth_vec, estimate_vec),
    spatiotemporal_soft_dice = soft_dice(truth_vec, estimate_vec),
    peak_spatial_cor = safe_cor(truth_mat[peak_idx, ], estimate_mat[peak_idx, ]),
    peak_spatial_rmse = sqrt(mean((estimate_mat[peak_idx, ] -
                                   truth_mat[peak_idx, ])^2)),
    peak_spatial_mae = mean(abs(estimate_mat[peak_idx, ] - truth_mat[peak_idx, ])),
    peak_spatial_cosine = cosine_similarity(truth_mat[peak_idx, ],
                                            estimate_mat[peak_idx, ]),
    peak_spatial_soft_dice = soft_dice(truth_mat[peak_idx, ],
                                       estimate_mat[peak_idx, ]),
    truth_peak_uV = max(truth_vec),
    estimate_peak_uV = max(estimate_vec),
    peak_recovery_pct = 100 * max(estimate_vec) / max(truth_vec),
    truth_weighted_mean_truth_uV = truth_weighted_mean(truth_vec, truth_vec),
    truth_weighted_mean_estimate_uV = truth_weighted_mean(estimate_vec, truth_vec),
    truth_weighted_mean_recovery_pct =
        100 * truth_weighted_mean(estimate_vec, truth_vec) /
        truth_weighted_mean(truth_vec, truth_vec),
    truth_top25_mean_truth_uV = mean(truth_vec[truth_top25_idx]),
    truth_top25_mean_estimate_uV = mean(estimate_vec[truth_top25_idx]),
    truth_top25_mean_recovery_pct =
        100 * mean(estimate_vec[truth_top25_idx]) /
        mean(truth_vec[truth_top25_idx]),
    n_truth_weight_points = sum(truth_weight_idx),
    n_truth_top25_points = sum(truth_top25_idx)
)

full <- expand.grid(
    time = times,
    channel = channels,
    KEEP.OUT.ATTRS = FALSE
) %>%
    mutate(
        truth_effect = truth_vec,
        estimate_effect = estimate_vec,
        residual = estimate_effect - truth_effect
    )

summary_path <- file.path(out_dir, sprintf("summary_continuous_%s.csv", tag))
full_path <- file.path(out_dir, sprintf("full_continuous_%s.csv", tag))
write.csv(summary, summary_path, row.names = FALSE)
write.csv(full, full_path, row.names = FALSE)
cat("Saved:", summary_path, "\n")
cat("Saved:", full_path, "\n")
print(summary)
