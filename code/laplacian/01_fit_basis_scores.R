###########################################################
# Fit temporal GAMMs to Laplacian basis scores.
#
# Each participant x condition x time scalp map is projected onto one
# spatial Laplacian basis vector. The resulting score is fitted with brms.
###########################################################

run_fit <- function() {
    out_dir <- file.path(project_root(), "results", out_subdir)
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

    # Build graph-Laplacian basis from the 61-channel adjacency matrix.
    adjacency_file <- file.path(project_root(), "archive", "models_original_2026_05_04",
                                "adjacency_xproj_yproj_61ch.rds")
    adj <- readRDS(adjacency_file)
    adj <- as.matrix(adj)
    adj <- ((adj + t(adj)) > 0) * 1
    diag(adj) <- 0
    channels <- rownames(adj)

    degree <- rowSums(adj)
    laplacian <- diag(degree) - adj
    eig <- eigen(laplacian, symmetric = TRUE)

    # Low eigenvalues correspond to smoother spatial basis functions.
    ord <- order(eig$values)
    eig_values <- eig$values[ord]
    eig_vectors <- eig$vectors[, ord, drop = FALSE]
    rownames(eig_vectors) <- channels

    # Trial-level data are aggregated before projection.
    sim_data <- readRDS(data_path(data_subdir, sprintf("sim_%03d.rds", seed)))
    sim_data_model <- aggregate_trials(sim_data)

    # Put channels into columns so each row is one scalp map.
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

    make_score_data <- function(task_id) {
        basis_vector <- eig_vectors[, task_id]
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
            basis_vector = basis_vector,
            basis_label = sprintf("lap%02d_lambda%.6f", task_id, eig_values[task_id]),
            participant_id = score_data$participant[1],
            score_data = score_data
        )
    }

    summarise_fit <- function(fit, meta, elapsed_sec, fit_file) {
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

    # First basis uses brm(); later bases reuse the compiled Stan model.
    fit_model <- NULL
    all_results <- list()

    for (task_id in task_ids) {
        meta <- make_score_data(task_id)
        message(sprintf("Starting %s with adapt_delta=%.3f", meta$basis_label, adapt_delta))
        start_time <- Sys.time()

        if (is.null(fit_model)) {
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

    combined_file <- file.path(out_dir, "summary_all.csv")
    write.csv(dplyr::bind_rows(all_results), combined_file, row.names = FALSE)
    message(sprintf("Wrote combined summary: %s", combined_file))
}
