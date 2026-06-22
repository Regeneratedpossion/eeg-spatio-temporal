###########################################################
# Fit temporal GAMMs to Laplacian basis scores.
#
# Each participant x condition x time scalp map is projected onto one
# spatial Laplacian basis vector. The resulting score is fitted with brms.
###########################################################

run_fit <- function() {
    # All basis fits and summaries are written to one fixed result folder.
    out_dir <- file.path(project_root(), "results", out_subdir)
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

    # Step 1: read the sensor adjacency matrix.
    # This matrix says which EEG sensors are neighbours on the scalp.
    adjacency_file <- file.path(project_root(), "archive", "models_original_2026_05_04",
                                "adjacency_xproj_yproj_61ch.rds")
    adj <- readRDS(adjacency_file)

    # Clean it into an undirected 0/1 graph.
    # After this, adj[i, j] = 1 means sensor i and sensor j are neighbours.
    adj <- as.matrix(adj)
    adj <- ((adj + t(adj)) > 0) * 1
    diag(adj) <- 0
    channels <- rownames(adj)

    # Step 2: graph Laplacian decomposition.
    # L = D - A, where D is the number of neighbours for each sensor.
    # Eigenvectors of L are spatial basis functions on the sensor graph.
    degree <- rowSums(adj)
    laplacian <- diag(degree) - adj
    eig <- eigen(laplacian, symmetric = TRUE)

    # Low eigenvalues correspond to smoother spatial patterns.
    # K=20 keeps the constant/global component plus the first smooth patterns.
    ord <- order(eig$values)
    eig_values <- eig$values[ord]
    eig_vectors <- eig$vectors[, ord, drop = FALSE]
    rownames(eig_vectors) <- channels

    # Step 3: aggregate trial-level data before spatial projection.
    # After aggregation, each row is participant x condition x time x channel.
    sim_data <- readRDS(data_path(data_subdir, sprintf("sim_%03d.rds", seed)))
    sim_data_model <- aggregate_trials(sim_data)

    # Step 4: put channels into columns.
    # Each row now contains one full 61-channel scalp map at one time point.
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
        # Step 5: project the scalp map onto one Laplacian basis vector.
        # score = low-dimensional representation of the 61-channel map.
        basis_vector <- eig_vectors[, task_id]
        score <- as.numeric(eeg_mean_by_channel %*% basis_vector)

        # Propagate SEM through the same linear projection.
        # Because SEM is uncertainty, variances are projected and then square-rooted.
        score_sem <- sqrt(as.numeric((eeg_sem_by_channel^2) %*% (basis_vector^2)))

        # brms fits one temporal GAMM to this basis score.
        # condition uses sum coding, so the condition coefficient is centred.
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
        # Minimal diagnostics needed to decide whether this basis fit is usable.
        sampler_info <- nuts_params(fit)
        draw_summary <- posterior::summarise_draws(as_draws_df(fit))
        divergent <- sampler_info$Value[sampler_info$Parameter == "divergent__"]
        treedepth <- sampler_info$Value[sampler_info$Parameter == "treedepth__"]

        data.frame(
            basis_task = meta$task_id,
            basis_label = meta$basis_label,
            elapsed_sec = elapsed_sec,
            n_divergent = sum(divergent),
            n_treedepth = sum(treedepth >= max_treedepth),
            max_rhat = max(draw_summary$rhat, na.rm = TRUE),
            min_ess_bulk = min(draw_summary$ess_bulk, na.rm = TRUE),
            fit_file = fit_file
        )
    }

    # Step 6: fit one temporal GAMM per spatial basis score.
    # The first basis uses brm() and compiles Stan.
    # Later bases use update(..., recompile = FALSE), which saves compile time.
    fit_model <- NULL
    all_results <- list()
    fit_control <- list(adapt_delta = adapt_delta, max_treedepth = max_treedepth)

    for (task_id in task_ids) {
        meta <- make_score_data(task_id)
        message(sprintf("Starting %s with adapt_delta=%.3f", meta$basis_label, adapt_delta))
        start_time <- Sys.time()

        if (is.null(fit_model)) {
            # First model: full brms call, including Stan compilation.
            fit <- brm(
                formula = model_formula, data = meta$score_data,
                family = gaussian(), chains = 4, cores = 4,
                iter = iter, warmup = warmup, seed = seed + meta$task_id,
                threads = threading(threads = 4),
                backend = "cmdstanr", refresh = 100,
                control = fit_control
            )
        } else {
            # Later models: same formula, different basis-score data.
            fit <- update(
                fit_model, newdata = meta$score_data,
                chains = 4, cores = 4,
                iter = iter, warmup = warmup, seed = seed + meta$task_id,
                refresh = 100, recompile = FALSE,
                control = fit_control
            )
        }
        fit_model <- fit

        # Save only what reconstruction needs:
        # fitted brms object, basis vector, channel order, and participant id.
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

    # This CSV is used by the reconstruction script to find all basis fits.
    combined_file <- file.path(out_dir, "summary_all.csv")
    write.csv(dplyr::bind_rows(all_results), combined_file, row.names = FALSE)
    message(sprintf("Wrote combined summary: %s", combined_file))
}
