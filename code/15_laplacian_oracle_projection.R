###########################################################
# 15_laplacian_oracle_projection.R
#
# Oracle projection diagnostic for Laplacian spatial bases.
# Projects the known continuous Gaussian truth field onto the
# selected spatial basis and reports the best possible recovery
# before fitting any Bayesian model.
###########################################################

suppressPackageStartupMessages({
    library(dplyr)
})

source(file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()), "code/utils.R"))

data_subdir <- Sys.getenv("LAPLACE_DATA_SUBDIR",
                          "montecarlo_p20_gaussian_hardlabel")
seed <- as.integer(Sys.getenv("LAPLACE_SEED", "1"))
k_values <- as.integer(strsplit(Sys.getenv("LAPLACE_K_VALUES", "10,15,20,30"),
                                ",")[[1]])

adj_file <- file.path(project_root(), "archive", "models_original_2026_05_04",
                      "adjacency_xproj_yproj_61ch.rds")
adj <- as.matrix(readRDS(adj_file))
adj <- ((adj + t(adj)) > 0) * 1
diag(adj) <- 0
channels <- rownames(adj)

laplacian <- diag(rowSums(adj)) - adj
eig <- eigen(laplacian, symmetric = TRUE)
ord <- order(eig$values)
eig_values <- eig$values[ord]
eig_vectors <- eig$vectors[, ord, drop = FALSE]
rownames(eig_vectors) <- channels

gt <- readRDS(data_path(data_subdir, sprintf("gt_%03d.rds", seed)))
sim <- readRDS(data_path(data_subdir, sprintf("sim_%03d.rds", seed)))
sim_model <- aggregate_trials(sim)
times <- sort(unique(sim_model$time))

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

sensors <- sim_model %>%
    distinct(channel, xproj, yproj) %>%
    filter(channel %in% channels)
seed_coord <- sensors %>% filter(channel == gt$seed_electrode)
sensors <- sensors %>%
    mutate(
        spatial_distance = sqrt((xproj - seed_coord$xproj)^2 +
                                (yproj - seed_coord$yproj)^2),
        spatial_weight_raw = exp(-(spatial_distance^2) /
                                 (2 * gt$sigma_space^2))
    ) %>%
    arrange(match(channel, channels))

truth_mat <- outer(truth_temporal_bump(times, gt), sensors$spatial_weight_raw)
truth_vec <- as.vector(truth_mat)

truth_weighted_mean <- function(x, truth) {
    idx <- truth > 0
    sum(x[idx] * truth[idx]) / sum(truth[idx])
}

safe_cor <- function(x, y) {
    if (sd(x) == 0 || sd(y) == 0) return(NA_real_)
    cor(x, y)
}

evaluate_basis <- function(k_basis, include_constant) {
    non_constant <- which(eig_values > 1e-8)
    selected <- if (include_constant) {
        seq_len(k_basis)
    } else {
        non_constant[seq_len(k_basis)]
    }
    basis <- eig_vectors[, selected, drop = FALSE]
    projection <- basis %*% t(basis)
    estimate_mat <- truth_mat %*% projection
    estimate_vec <- as.vector(estimate_mat)
    truth_top25_idx <- truth_vec >= 0.25 * max(truth_vec)

    data.frame(
        dataset = data_subdir,
        seed = seed,
        k_basis = k_basis,
        include_constant = include_constant,
        first_eigenvalue = min(eig_values[selected]),
        max_eigenvalue = max(eig_values[selected]),
        peak_recovery_pct = 100 * max(estimate_vec) / max(truth_vec),
        truth_weighted_mean_recovery_pct =
            100 * truth_weighted_mean(estimate_vec, truth_vec) /
            truth_weighted_mean(truth_vec, truth_vec),
        truth_top25_mean_recovery_pct =
            100 * mean(estimate_vec[truth_top25_idx]) /
            mean(truth_vec[truth_top25_idx]),
        spatiotemporal_cor = safe_cor(truth_vec, estimate_vec),
        rmse = sqrt(mean((estimate_vec - truth_vec)^2))
    )
}

summary <- do.call(rbind, lapply(k_values, function(k_basis) {
    rbind(
        evaluate_basis(k_basis, include_constant = FALSE),
        evaluate_basis(k_basis, include_constant = TRUE)
    )
}))

out_dir <- file.path(project_root(), "results", "diagnostics",
                     "laplacian_oracle_projection")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(out_dir, sprintf("oracle_projection_%s_seed%03d.csv",
                                       data_subdir, seed))
write.csv(summary, out_file, row.names = FALSE)
cat("Saved:", out_file, "\n")
print(summary)
