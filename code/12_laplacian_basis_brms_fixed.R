###########################################################
# 12_laplacian_basis_brms_fixed.R
#
# Fixed Laplacian-basis run:
# - runs basis tasks sequentially in one R process
# - compiles the brms/CmdStan model once
# - reuses the compiled model with update(..., recompile = FALSE)
# - uses a higher adapt_delta to reduce divergent transitions
###########################################################

source(file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()), "code/utils.R"))

setup_cmdstan_cache()

k_basis <- as.integer(Sys.getenv("LAPLACE_K", "10"))
include_constant_basis <- tolower(Sys.getenv("LAPLACE_INCLUDE_CONSTANT", "false")) %in%
    c("true", "1", "yes")
task_spec <- Sys.getenv("LAPLACE_TASKS", unset = sprintf("1:%d", k_basis))
task_ids <- eval(parse(text = task_spec))
task_ids <- as.integer(task_ids)
stopifnot(all(task_ids >= 1), all(task_ids <= k_basis))

data_subdir <- Sys.getenv("LAPLACE_DATA_SUBDIR",
                          "montecarlo_p5_gaussian_hardlabel")
seed <- as.integer(Sys.getenv("LAPLACE_SEED", "1"))
iter <- as.integer(Sys.getenv("LAPLACE_ITER", "3000"))
warmup <- as.integer(Sys.getenv("LAPLACE_WARMUP", "1000"))
time_k <- as.integer(Sys.getenv("LAPLACE_TIME_K", "10"))
adapt_delta <- as.numeric(Sys.getenv("LAPLACE_ADAPT_DELTA", "0.999"))
max_treedepth <- as.integer(Sys.getenv("LAPLACE_MAX_TREEDEPTH", "13"))
random_effect_mode <- Sys.getenv("LAPLACE_RANDOM_EFFECT", "slope")
model_formula <- switch(
    random_effect_mode,
    slope = bf(score | se(score_sem) ~
                   condition +
                   s(time, by = condition, k = time_k) +
                   (1 + condition | participant)),
    intercept = bf(score | se(score_sem) ~
                       condition +
                       s(time, by = condition, k = time_k) +
                       (1 | participant)),
    none = bf(score | se(score_sem) ~
                  condition +
                  s(time, by = condition, k = time_k)),
    stop("LAPLACE_RANDOM_EFFECT must be one of: slope, intercept, none.")
)

adj_file <- file.path(project_root(), "archive", "models_original_2026_05_04",
                      "adjacency_xproj_yproj_61ch.rds")
out_subdir <- Sys.getenv("LAPLACE_OUT_SUBDIR", "laplacian_basis_brms_fixed")
out_dir <- file.path(project_root(), "results", out_subdir)
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
selected_basis_indices <- if (include_constant_basis) {
    seq_len(k_basis)
} else {
    non_constant[seq_len(k_basis)]
}
if (length(selected_basis_indices) < k_basis || any(is.na(selected_basis_indices))) {
    stop("Not enough graph-Laplacian eigenvectors.")
}

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

make_score_data <- function(task_id) {
    basis_index <- selected_basis_indices[task_id]
    basis_vector <- eig_vectors[, basis_index]
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

    list(
        task_id = task_id,
        basis_index = basis_index,
        basis_vector = basis_vector,
        basis_label = sprintf("lap%02d_lambda%.6f", task_id, eig_values[basis_index]),
        eigenvalue = eig_values[basis_index],
        score_data = score_data
    )
}

summarise_fit <- function(fit, meta, elapsed_sec, fit_file) {
    np <- nuts_params(fit)
    n_div <- sum(np$Value[np$Parameter == "divergent__"])
    n_tree <- sum(np$Value[np$Parameter == "treedepth__"] >= max_treedepth)
    summ <- posterior::summarise_draws(as_draws_df(fit))

    data.frame(
        basis_task = meta$task_id,
        basis_label = meta$basis_label,
        basis_index = meta$basis_index,
        eigenvalue = meta$eigenvalue,
        status = "ok",
        seed = seed,
        k_basis = k_basis,
        time_k = time_k,
        iter = iter,
        warmup = warmup,
        adapt_delta = adapt_delta,
        max_treedepth = max_treedepth,
        random_effect_mode = random_effect_mode,
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
}

fit_model <- NULL
all_results <- list()

for (task_id in task_ids) {
    meta <- make_score_data(task_id)
    message(sprintf("Starting %s with adapt_delta=%.3f", meta$basis_label, adapt_delta))
    start_time <- Sys.time()

    result <- tryCatch({
        if (is.null(fit_model)) {
            fit <- brm(
                formula = model_formula,
                data = meta$score_data,
                family = gaussian(),
                chains = 4,
                cores = 4,
                iter = iter,
                warmup = warmup,
                seed = seed + task_id,
                threads = threading(threads = 4),
                backend = "cmdstanr",
                refresh = 100,
                control = list(adapt_delta = adapt_delta,
                               max_treedepth = max_treedepth)
            )
        } else {
            fit <- update(
                fit_model,
                newdata = meta$score_data,
                chains = 4,
                cores = 4,
                iter = iter,
                warmup = warmup,
                seed = seed + task_id,
                refresh = 100,
                recompile = FALSE,
                control = list(adapt_delta = adapt_delta,
                               max_treedepth = max_treedepth)
            )
        }
        fit_model <- fit

        fit_file <- file.path(out_dir, sprintf("fit_%s_iter%d_warmup%d_adapt%s.rds",
                                               meta$basis_label, iter, warmup,
                                               gsub("\\.", "p", as.character(adapt_delta))))
        saveRDS(list(
            fit = fit,
            basis_label = meta$basis_label,
            basis_task = meta$task_id,
            basis_index = meta$basis_index,
            eigenvalue = meta$eigenvalue,
            basis_vector = meta$basis_vector,
            channels = channels,
            k_basis = k_basis,
            score_data = meta$score_data,
            adapt_delta = adapt_delta,
            max_treedepth = max_treedepth
        ), fit_file)

        summarise_fit(
            fit = fit,
            meta = meta,
            elapsed_sec = as.numeric(Sys.time() - start_time, units = "secs"),
            fit_file = fit_file
        )
    }, error = function(e) {
        data.frame(
            basis_task = meta$task_id,
            basis_label = meta$basis_label,
            basis_index = meta$basis_index,
            eigenvalue = meta$eigenvalue,
            status = "error",
            seed = seed,
            k_basis = k_basis,
            time_k = time_k,
            iter = iter,
            warmup = warmup,
            adapt_delta = adapt_delta,
            max_treedepth = max_treedepth,
            random_effect_mode = random_effect_mode,
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

    out_file <- file.path(out_dir, sprintf("summary_%s_iter%d_warmup%d_adapt%s.csv",
                                           meta$basis_label, iter, warmup,
                                           gsub("\\.", "p", as.character(adapt_delta))))
    write.csv(result, out_file, row.names = FALSE)
    print(result)
    all_results[[length(all_results) + 1]] <- result
}

combined_file <- file.path(out_dir, sprintf("summary_all_iter%d_warmup%d_adapt%s.csv",
                                            iter, warmup,
                                            gsub("\\.", "p", as.character(adapt_delta))))
write.csv(dplyr::bind_rows(all_results), combined_file, row.names = FALSE)
message(sprintf("Wrote combined summary: %s", combined_file))
