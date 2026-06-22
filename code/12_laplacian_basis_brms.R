###########################################################
# 12_laplacian_basis_brms.R
#
# Adjacency-informed low-rank spatial basis model.
# Projects 61-channel scalp maps onto graph-Laplacian eigenvectors,
# then fits one Bayesian temporal GAMM per spatial basis score.
###########################################################

source(file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()), "code/utils.R"))

setup_cmdstan_cache()

task_id <- as.integer(Sys.getenv(
    "SLURM_ARRAY_TASK_ID",
    unset = as.character(as.integer(c(commandArgs(trailingOnly = TRUE), 1)[1]))
))

k_basis <- as.integer(Sys.getenv("LAPLACE_K", "10"))
stopifnot(task_id >= 1, task_id <= k_basis)

data_subdir <- Sys.getenv("LAPLACE_DATA_SUBDIR",
                          "montecarlo_p5_gaussian_hardlabel")
seed <- as.integer(Sys.getenv("LAPLACE_SEED", "1"))
iter <- as.integer(Sys.getenv("LAPLACE_ITER", "3000"))
warmup <- as.integer(Sys.getenv("LAPLACE_WARMUP", "1000"))
time_k <- as.integer(Sys.getenv("LAPLACE_TIME_K", "10"))

adj_file <- file.path(project_root(), "archive", "models_original_2026_05_04",
                      "adjacency_xproj_yproj_61ch.rds")
out_dir <- file.path(project_root(), "results", "laplacian_basis_brms")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

adj <- readRDS(adj_file)
adj <- as.matrix(adj)
adj <- ((adj + t(adj)) > 0) * 1
diag(adj) <- 0
channels <- rownames(adj)
stopifnot(!is.null(channels), identical(channels, colnames(adj)))

degree <- rowSums(adj)
laplacian <- diag(degree) - adj
eig <- eigen(laplacian, symmetric = TRUE)
ord <- order(eig$values)
eig_values <- eig$values[ord]
eig_vectors <- eig$vectors[, ord, drop = FALSE]
rownames(eig_vectors) <- channels

non_constant <- which(eig_values > 1e-8)
if (length(non_constant) < k_basis) {
    stop("Not enough non-constant graph-Laplacian eigenvectors.")
}
basis_index <- non_constant[task_id]
basis_vector <- eig_vectors[, basis_index]
basis_label <- sprintf("lap%02d_lambda%.6f", task_id, eig_values[basis_index])

sim_data <- readRDS(data_path(data_subdir, sprintf("sim_%03d.rds", seed)))
sim_data_model <- aggregate_trials(sim_data)

missing_channels <- setdiff(channels, unique(sim_data_model$channel))
extra_channels <- setdiff(unique(sim_data_model$channel), channels)
if (length(missing_channels) || length(extra_channels)) {
    stop(sprintf("Channel mismatch. Missing: %s Extra: %s",
                 paste(missing_channels, collapse = ","),
                 paste(extra_channels, collapse = ",")))
}

mean_wide <- sim_data_model %>%
    select(participant, condition, time, channel, eeg_mean) %>%
    tidyr::pivot_wider(names_from = channel, values_from = eeg_mean) %>%
    arrange(participant, condition, time)

sem_wide <- sim_data_model %>%
    select(participant, condition, time, channel, eeg_sem) %>%
    tidyr::pivot_wider(names_from = channel, values_from = eeg_sem) %>%
    arrange(participant, condition, time)

Y <- as.matrix(mean_wide[, channels])
S <- as.matrix(sem_wide[, channels])

score <- as.numeric(Y %*% basis_vector)
score_sem <- sqrt(as.numeric((S^2) %*% (basis_vector^2)))
score_sem[!is.finite(score_sem) | score_sem <= 0] <- min(score_sem[score_sem > 0],
                                                         na.rm = TRUE)

score_data <- mean_wide %>%
    select(participant, condition, time) %>%
    mutate(
        score = score,
        score_sem = score_sem,
        condition = factor(condition, levels = c("cond1", "cond2"))
    )
contrasts(score_data$condition) <- c(-0.5, 0.5)

start_time <- Sys.time()

result <- tryCatch({
    fit <- brm(
        formula = bf(score | se(score_sem) ~
                         condition +
                         s(time, by = condition, k = time_k) +
                         (1 + condition | participant)),
        data = score_data,
        family = gaussian(),
        chains = 4,
        cores = 4,
        iter = iter,
        warmup = warmup,
        seed = seed + task_id,
        threads = threading(threads = 4),
        backend = "cmdstanr",
        refresh = 100,
        control = list(adapt_delta = 0.99, max_treedepth = 13)
    )

    fit_file <- file.path(out_dir, sprintf("fit_%s_iter%d_warmup%d.rds",
                                           basis_label, iter, warmup))
    saveRDS(list(
        fit = fit,
        basis_label = basis_label,
        basis_task = task_id,
        basis_index = basis_index,
        eigenvalue = eig_values[basis_index],
        basis_vector = basis_vector,
        channels = channels,
        k_basis = k_basis,
        score_data = score_data
    ), fit_file)

    elapsed_sec <- as.numeric(Sys.time() - start_time, units = "secs")
    np <- nuts_params(fit)
    n_div <- sum(np$Value[np$Parameter == "divergent__"])
    max_td <- fit$fit@stan_args[[1]]$control$max_treedepth
    n_tree <- sum(np$Value[np$Parameter == "treedepth__"] >= max_td)
    summ <- posterior::summarise_draws(as_draws_df(fit))

    data.frame(
        basis_task = task_id,
        basis_label = basis_label,
        basis_index = basis_index,
        eigenvalue = eig_values[basis_index],
        status = "ok",
        seed = seed,
        k_basis = k_basis,
        time_k = time_k,
        iter = iter,
        warmup = warmup,
        elapsed_sec = elapsed_sec,
        n_divergent = n_div,
        n_treedepth = n_tree,
        max_rhat = max(summ$rhat, na.rm = TRUE),
        n_rhat_gt_1_01 = sum(summ$rhat > 1.01, na.rm = TRUE),
        min_ess_bulk = min(summ$ess_bulk, na.rm = TRUE),
        mean_ess_bulk = mean(summ$ess_bulk, na.rm = TRUE),
        min_ess_tail = min(summ$ess_tail, na.rm = TRUE),
        fit_file = fit_file,
        error = NA_character_
    )
}, error = function(e) {
    data.frame(
        basis_task = task_id,
        basis_label = basis_label,
        basis_index = basis_index,
        eigenvalue = eig_values[basis_index],
        status = "error",
        seed = seed,
        k_basis = k_basis,
        time_k = time_k,
        iter = iter,
        warmup = warmup,
        elapsed_sec = as.numeric(Sys.time() - start_time, units = "secs"),
        n_divergent = NA_integer_,
        n_treedepth = NA_integer_,
        max_rhat = NA_real_,
        n_rhat_gt_1_01 = NA_integer_,
        min_ess_bulk = NA_real_,
        mean_ess_bulk = NA_real_,
        min_ess_tail = NA_real_,
        fit_file = NA_character_,
        error = conditionMessage(e)
    )
})

out_file <- file.path(out_dir, sprintf("summary_%s_iter%d_warmup%d.csv",
                                       basis_label, iter, warmup))
write.csv(result, out_file, row.names = FALSE)
print(result)
