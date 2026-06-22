###########################################################
# 12_laplacian_basis_brms_fixed.R
#
# Single Laplacian-basis workflow.
#
# Actions:
#   LAPLACE_ACTION=fit                 fit basis-score temporal GAMMs
#   LAPLACE_ACTION=diagnose_hard       threshold-based recovery summary
#   LAPLACE_ACTION=diagnose_continuous continuous Gaussian-field recovery
#
# The fit action:
# - runs basis tasks sequentially in one R process
# - compiles the brms/CmdStan model once
# - reuses the compiled model with update(..., recompile = FALSE)
###########################################################

source(file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()), "code/utils.R"))

setup_cmdstan_cache()

laplace_action <- Sys.getenv("LAPLACE_ACTION", "fit")

run_fit <- function() {
# K=20 + constant Laplacian basis model.
k_basis <- 20L
task_ids <- seq_len(k_basis)
data_subdir <- "montecarlo_p20_gaussian_hardlabel"
out_subdir <- "laplacian_basis_brms_p20_gaussian_hardlabel_constant_long/k20"
seed <- 1L
iter <- 3000L
warmup <- 1000L
time_k <- 10L
adapt_delta <- 0.999
max_treedepth <- 13L

model_formula <- bf(score | se(score_sem) ~
                        condition +
                        s(time, by = condition, k = time_k) +
                        (1 + condition | participant))

# Output folder.
out_dir <- file.path(project_root(), "results", out_subdir)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Graph-Laplacian spatial basis from the sensor adjacency matrix.
adjacency_file <- file.path(project_root(), "archive", "models_original_2026_05_04",
                            "adjacency_xproj_yproj_61ch.rds")
adj <- readRDS(adjacency_file)

# Clean adjacency to an undirected 0/1 graph with no self-connections.
adj <- as.matrix(adj)
adj <- ((adj + t(adj)) > 0) * 1
diag(adj) <- 0
channels <- rownames(adj)

# Graph Laplacian: L = D - A.
degree <- rowSums(adj)
laplacian <- diag(degree) - adj
eig <- eigen(laplacian, symmetric = TRUE)

# Low eigenvalues give smoother spatial basis functions.
ord <- order(eig$values)
eig_values <- eig$values[ord]
eig_vectors <- eig$vectors[, ord, drop = FALSE]
rownames(eig_vectors) <- channels

# Trial data are first aggregated to participant x condition x time x channel.
sim_data <- readRDS(data_path(data_subdir, sprintf("sim_%03d.rds", seed)))
sim_data_model <- aggregate_trials(sim_data)

# Put the 61 channels into columns so each row is one scalp map.
wide_by_channel <- function(value_col) {
    sim_data_model %>%
        select(participant, condition, time, channel, all_of(value_col)) %>%
        tidyr::pivot_wider(names_from = channel, values_from = all_of(value_col)) %>%
        arrange(participant, condition, time)
}

mean_wide <- wide_by_channel("eeg_mean")
sem_wide <- wide_by_channel("eeg_sem")

eeg_mean_by_channel <- as.matrix(mean_wide[, channels])
eeg_sem_by_channel <- as.matrix(sem_wide[, channels])

# Project one scalp map onto one Laplacian basis vector.
make_score_data <- function(task_id) {
    basis_index <- task_id
    basis_vector <- eig_vectors[, basis_index]
    score <- as.numeric(eeg_mean_by_channel %*% basis_vector)
    score_sem <- sqrt(as.numeric((eeg_sem_by_channel^2) %*% (basis_vector^2)))

    score_data <- mean_wide %>%
        select(participant, condition, time) %>%
        mutate(
            score = score,
            score_sem = score_sem,
            condition = factor(condition, levels = c("cond1", "cond2"))
        )
    contrasts(score_data$condition) <- c(-0.5, 0.5)

    list(
        task_id = task_id,
        basis_index = basis_index,
        basis_vector = basis_vector,
        basis_label = sprintf("lap%02d_lambda%.6f", task_id, eig_values[basis_index]),
        eigenvalue = eig_values[basis_index],
        participant_id = score_data$participant[1],
        score_data = score_data
    )
}

summarise_fit <- function(fit, meta, elapsed_sec, fit_file) {
    # One-row sampler diagnostic summary for this basis model.
    sampler_info <- nuts_params(fit)
    draw_summary <- posterior::summarise_draws(as_draws_df(fit))

    data.frame(
        basis_task = meta$task_id,
        basis_label = meta$basis_label,
        elapsed_sec = elapsed_sec,
        n_divergent = sum(sampler_info$Value[sampler_info$Parameter == "divergent__"]),
        n_treedepth = sum(sampler_info$Value[sampler_info$Parameter == "treedepth__"] >=
                              max_treedepth),
        max_rhat = max(draw_summary$rhat, na.rm = TRUE),
        min_ess_bulk = min(draw_summary$ess_bulk, na.rm = TRUE),
        fit_file = fit_file
    )
}

# Fit one brms temporal GAMM for each Laplacian basis score.
# The first basis is fitted with a regular brm() call.
# Later bases use update(..., recompile = FALSE) to reuse the compiled Stan model.
fit_model <- NULL
all_results <- list()

for (task_id in task_ids) {
    meta <- make_score_data(task_id)
    message(sprintf("Starting %s with adapt_delta=%.3f", meta$basis_label, adapt_delta))
    start_time <- Sys.time()

    if (is.null(fit_model)) {
        # Regular brms model fit for the first spatial basis score.
        fit <- brm(
            formula = model_formula,
            data = meta$score_data,
            family = gaussian(),
            chains = 4,
            cores = 4,
            iter = iter,
            warmup = warmup,
            seed = seed + meta$task_id,
            threads = threading(threads = 4),
            backend = "cmdstanr",
            refresh = 100,
            control = list(adapt_delta = adapt_delta,
                           max_treedepth = max_treedepth)
        )
    } else {
        # Same brms model, new basis-score data, no Stan recompilation.
        fit <- update(
            fit_model,
            newdata = meta$score_data,
            chains = 4,
            cores = 4,
            iter = iter,
            warmup = warmup,
            seed = seed + meta$task_id,
            refresh = 100,
            recompile = FALSE,
            control = list(adapt_delta = adapt_delta,
                           max_treedepth = max_treedepth)
        )
    }
    fit_model <- fit

    # Save the fitted basis model and record its path for reconstruction.
    fit_file <- file.path(out_dir, paste0("fit_", meta$basis_label, ".rds"))
    saveRDS(list(
        fit = fit,
        basis_task = meta$task_id,
        basis_vector = meta$basis_vector,
        channels = channels,
        participant_id = meta$participant_id
    ), fit_file)

    result <- summarise_fit(
        fit = fit,
        meta = meta,
        elapsed_sec = as.numeric(Sys.time() - start_time, units = "secs"),
        fit_file = fit_file
    )

    print(result)
    all_results[[length(all_results) + 1]] <- result
}

# This combined summary is used later to find all basis fit files for reconstruction.
combined_file <- file.path(out_dir, "summary_all.csv")
write.csv(dplyr::bind_rows(all_results), combined_file, row.names = FALSE)
message(sprintf("Wrote combined summary: %s", combined_file))
}

load_laplace_reconstruction <- function() {
    fit_dir <- file.path(project_root(), "results", out_subdir)
    adapt_tag <- gsub("\\.", "p", as.character(adapt_delta))
    tag <- sprintf("laplacian_k%d_seed%03d_iter%d_warmup%d_adapt%s",
                   k_basis, seed, iter, warmup, adapt_tag)

    summary_file <- file.path(fit_dir, "summary_all.csv")
    conv_summary <- read.csv(summary_file)
    fit_files <- conv_summary$fit_file[order(conv_summary$basis_task)]

    gt <- readRDS(data_path(data_subdir, sprintf("gt_%03d.rds", seed)))
    sim <- readRDS(data_path(data_subdir, sprintf("sim_%03d.rds", seed)))
    sim_model <- aggregate_trials(sim)
    times <- sort(unique(sim_model$time))
    cond_levels <- c("cond1", "cond2")

    # Posterior condition effect for one fitted basis-score model.
    condition_score_draws <- function(obj, times) {
        base <- data.frame(
            time = times,
            score = 0,
            score_sem = 1,
            participant = obj$participant_id
        )
        nd1 <- base
        nd2 <- base
        nd1$condition <- factor("cond1", levels = cond_levels)
        nd2$condition <- factor("cond2", levels = cond_levels)

        posterior_epred(obj$fit, newdata = nd2, re_formula = NA) -
            posterior_epred(obj$fit, newdata = nd1, re_formula = NA)
    }

    first_obj <- readRDS(fit_files[1])
    channels <- first_obj$channels
    effect_array <- NULL

    # Reconstruct sensor x time posterior effect:
    # score effect for each basis x spatial basis vector, summed over bases.
    for (file in fit_files) {
        obj <- readRDS(file)
        score_diff <- condition_score_draws(obj, times)
        basis_vector <- obj$basis_vector[channels]

        if (is.null(effect_array)) {
            effect_array <- array(
                0,
                dim = c(nrow(score_diff), length(times), length(channels)),
                dimnames = list(NULL, as.character(times), channels)
            )
        }

        for (j in seq_along(channels)) {
            effect_array[, , j] <- effect_array[, , j] +
                score_diff * basis_vector[j]
        }
    }

    list(
        tag = tag,
        data_subdir = data_subdir,
        seed = seed,
        gt = gt,
        sim_model = sim_model,
        conv_summary = conv_summary,
        effect_array = effect_array,
        times = times,
        channels = channels,
        n_draws = dim(effect_array)[1]
    )
}

f1_metrics <- function(TP, FP, FN) {
    precision <- TP / (TP + FP)
    recall <- TP / (TP + FN)
    F1 <- 2 * precision * recall / (precision + recall)

    c(TP = TP, FP = FP, FN = FN,
      precision = precision, recall = recall, F1 = F1)
}

run_hard_diagnostics <- function() {
    rec <- load_laplace_reconstruction()
    out_dir <- file.path(project_root(), "results", "diagnostics",
                         "laplacian_basis_fixed")
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

    effect_array <- rec$effect_array
    times <- rec$times
    channels <- rec$channels
    gt <- rec$gt
    conv_summary <- rec$conv_summary
    mean_mat <- apply(effect_array, c(2, 3), mean)
    lo_mat <- apply(effect_array, c(2, 3), quantile, probs = 0.025)
    hi_mat <- apply(effect_array, c(2, 3), quantile, probs = 0.975)
    p_sig_mat <- apply(effect_array > ROPE_UV, c(2, 3), mean)
    p_floor <- 1 / rec$n_draws
    odds_mat <- pmax(p_sig_mat, p_floor) / pmax(1 - p_sig_mat, p_floor)

    full <- expand.grid(time = times, channel = channels, KEEP.OUT.ATTRS = FALSE) %>%
        mutate(
            mean_diff = as.vector(mean_mat),
            lo = as.vector(lo_mat),
            hi = as.vector(hi_mat),
            p_sig = as.vector(p_sig_mat),
            odds = as.vector(odds_mat),
            pred_active = odds > ODDS_THRESHOLD,
            truth_active = channel %in% gt$active_electrodes &
                time >= gt$true_onset & time <= gt$true_offset
        )

    seed_ch <- gt$seed_electrode
    seed_curve <- full %>% filter(channel == seed_ch) %>% arrange(time)
    clusters <- find_clusters(seed_curve$time, seed_curve$odds, ODDS_THRESHOLD)
    onset_est <- clusters$cluster_onset[1]
    offset_est <- clusters$cluster_offset[1]
    temporal_pass <- !is.na(onset_est) &&
        abs((onset_est - gt$true_onset) * 1000) <= BIAS_TOLERANCE_MS &&
        abs((offset_est - gt$true_offset) * 1000) <= BIAS_TOLERANCE_MS

    peak_time <- pick_time((gt$true_onset + gt$true_offset) / 2, times)
    peak_space <- full %>% filter(time == peak_time) %>% arrange(desc(mean_diff))
    pred_active_space <- peak_space$channel[peak_space$lo > 0]
    spatial <- f1_metrics(
        TP = length(intersect(pred_active_space, gt$active_electrodes)),
        FP = length(setdiff(pred_active_space, gt$active_electrodes)),
        FN = length(setdiff(gt$active_electrodes, pred_active_space))
    )
    st <- f1_metrics(
        TP = sum(full$pred_active & full$truth_active, na.rm = TRUE),
        FP = sum(full$pred_active & !full$truth_active, na.rm = TRUE),
        FN = sum(!full$pred_active & full$truth_active, na.rm = TRUE)
    )

    seed_idx <- match(seed_ch, channels)
    peak_idx <- match(peak_time, times)
    seed_peak_draws <- effect_array[, peak_idx, seed_idx]
    window_idx <- which(times >= gt$true_onset & times <= gt$true_offset)
    window_draws <- rowMeans(effect_array[, window_idx, seed_idx, drop = FALSE][, , 1])
    truth_time <- seq(gt$true_onset, gt$true_offset, length.out = 5000)
    truth_mean <- mean(truth_temporal_bump(truth_time, gt))
    cz_neighbours <- c("CZ", "FCZ", "CPZ", "C1", "C2", "FC1", "FC2",
                       "CP1", "CP2", "Cz", "FCz", "CPz")

    summary <- data.frame(
        fit_file = rec$tag,
        dataset = rec$data_subdir,
        convergence_pass = max(conv_summary$n_divergent, na.rm = TRUE) == 0 &&
            max(conv_summary$n_treedepth, na.rm = TRUE) == 0 &&
            max(conv_summary$max_rhat, na.rm = TRUE) < 1.01 &&
            min(conv_summary$min_ess_bulk, na.rm = TRUE) > 400,
        temporal_pass = temporal_pass,
        spatial_top_pass = peak_space$channel[1] %in% cz_neighbours,
        peak_ci_covers_truth = unname(quantile(seed_peak_draws, 0.025)) <=
            gt$bump_amplitude &&
            gt$bump_amplitude <= unname(quantile(seed_peak_draws, 0.975)),
        mean_ci_covers_truth = unname(quantile(window_draws, 0.025)) <=
            truth_mean && truth_mean <= unname(quantile(window_draws, 0.975)),
        spatial_F1 = spatial["F1"],
        spatiotemporal_F1 = st["F1"],
        peak_recovery_pct = 100 * mean(seed_peak_draws) / gt$bump_amplitude,
        mean_recovery_pct = 100 * mean(window_draws) / truth_mean
    )

    summary_path <- file.path(out_dir, sprintf("summary_%s.csv", rec$tag))
    write.csv(summary, summary_path, row.names = FALSE)
    cat("Saved:", summary_path, "\n")
    print(summary)
}

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
    truth_mat <- outer(truth_temporal_bump(times, ground_truth),
                       sensors$spatial_weight_raw)
    dimnames(truth_mat) <- list(as.character(times), sensors$channel)
    truth_mat
}

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

run_continuous_diagnostics <- function() {
    rec <- load_laplace_reconstruction()
    out_dir <- file.path(project_root(), "results", "diagnostics",
                         "laplacian_continuous")
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

    estimate_mat <- apply(rec$effect_array, c(2, 3), mean)
    truth_mat <- make_continuous_truth(rec$sim_model, rec$gt,
                                       rec$times, rec$channels)
    truth_vec <- as.vector(truth_mat)
    estimate_vec <- as.vector(estimate_mat)
    active_idx <- truth_vec > 0
    truth_top25_idx <- truth_vec >= 0.25 * max(truth_vec)
    peak_time <- pick_time((rec$gt$true_onset + rec$gt$true_offset) / 2,
                           rec$times)
    peak_idx <- match(peak_time, rec$times)

    summary <- data.frame(
        fit_file = rec$tag,
        dataset = rec$data_subdir,
        spatiotemporal_cor_all = safe_cor(truth_vec, estimate_vec),
        spatiotemporal_cor_effect_window = safe_cor(truth_vec[active_idx],
                                                    estimate_vec[active_idx]),
        spatiotemporal_rmse = sqrt(mean((estimate_vec - truth_vec)^2)),
        spatiotemporal_mae = mean(abs(estimate_vec - truth_vec)),
        spatiotemporal_cosine = cosine_similarity(truth_vec, estimate_vec),
        spatiotemporal_soft_dice = soft_dice(truth_vec, estimate_vec),
        peak_spatial_cor = safe_cor(truth_mat[peak_idx, ],
                                    estimate_mat[peak_idx, ]),
        truth_peak_uV = max(truth_vec),
        estimate_peak_uV = max(estimate_vec),
        peak_recovery_pct = 100 * max(estimate_vec) / max(truth_vec),
        truth_weighted_mean_recovery_pct =
            100 * truth_weighted_mean(estimate_vec, truth_vec) /
            truth_weighted_mean(truth_vec, truth_vec),
        truth_top25_mean_recovery_pct =
            100 * mean(estimate_vec[truth_top25_idx]) /
            mean(truth_vec[truth_top25_idx])
    )

    summary_path <- file.path(out_dir, sprintf("summary_continuous_%s.csv", rec$tag))
    write.csv(summary, summary_path, row.names = FALSE)
    cat("Saved:", summary_path, "\n")
    print(summary)
}

switch(
    laplace_action,
    fit = run_fit(),
    diagnose_hard = run_hard_diagnostics(),
    diagnose_continuous = run_continuous_diagnostics(),
    stop("LAPLACE_ACTION must be one of: fit, diagnose_hard, diagnose_continuous.")
)
