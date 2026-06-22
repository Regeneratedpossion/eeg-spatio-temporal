###########################################################
# 17_laplacian_score_diagnostics.R
#
# Data-side diagnostics for projected Laplacian basis scores.
# Used to inspect problematic basis tasks before changing the
# Bayesian model structure.
###########################################################

suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
})

source(file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()), "code/utils.R"))

data_subdir <- Sys.getenv("LAPLACE_DATA_SUBDIR",
                          "montecarlo_p20_gaussian_hardlabel")
seed <- as.integer(Sys.getenv("LAPLACE_SEED", "1"))
k_basis <- as.integer(Sys.getenv("LAPLACE_K", "20"))
include_constant_basis <- tolower(Sys.getenv("LAPLACE_INCLUDE_CONSTANT", "true")) %in%
    c("true", "1", "yes")
task_spec <- Sys.getenv("LAPLACE_TASKS", "1:20")
task_ids <- as.integer(eval(parse(text = task_spec)))

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

non_constant <- which(eig_values > 1e-8)
selected_basis_indices <- if (include_constant_basis) {
    seq_len(k_basis)
} else {
    non_constant[seq_len(k_basis)]
}

sim <- readRDS(data_path(data_subdir, sprintf("sim_%03d.rds", seed)))
sim_model <- aggregate_trials(sim)

mean_wide <- sim_model %>%
    select(participant, condition, time, channel, eeg_mean) %>%
    pivot_wider(names_from = channel, values_from = eeg_mean) %>%
    arrange(participant, condition, time)

sem_wide <- sim_model %>%
    select(participant, condition, time, channel, eeg_sem) %>%
    pivot_wider(names_from = channel, values_from = eeg_sem) %>%
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
    mean_wide %>%
        select(participant, condition, time) %>%
        mutate(
            basis_task = task_id,
            basis_index = basis_index,
            basis_label = sprintf("lap%02d_lambda%.6f", task_id,
                                  eig_values[basis_index]),
            eigenvalue = eig_values[basis_index],
            score = score,
            score_sem = score_sem,
            z_score = score / score_sem
        )
}

score_data <- bind_rows(lapply(task_ids, make_score_data))

score_summary <- score_data %>%
    group_by(basis_task, basis_label, eigenvalue) %>%
    summarise(
        n = n(),
        score_mean = mean(score),
        score_sd = sd(score),
        score_min = min(score),
        score_max = max(score),
        sem_mean = mean(score_sem),
        sem_min = min(score_sem),
        sem_max = max(score_sem),
        z_abs_max = max(abs(z_score)),
        z_abs_p99 = unname(quantile(abs(z_score), 0.99)),
        .groups = "drop"
    )

condition_summary <- score_data %>%
    group_by(basis_task, basis_label, condition, time) %>%
    summarise(mean_score = mean(score), .groups = "drop_last") %>%
    pivot_wider(names_from = condition, values_from = mean_score) %>%
    ungroup() %>%
    mutate(cond_diff = cond2 - cond1) %>%
    group_by(basis_task, basis_label) %>%
    summarise(
        max_abs_condition_diff = max(abs(cond_diff)),
        peak_condition_diff = cond_diff[which.max(abs(cond_diff))],
        peak_time = time[which.max(abs(cond_diff))],
        condition_diff_sd = sd(cond_diff),
        .groups = "drop"
    )

participant_summary <- score_data %>%
    group_by(basis_task, basis_label, participant, condition) %>%
    summarise(mean_score = mean(score), .groups = "drop_last") %>%
    pivot_wider(names_from = condition, values_from = mean_score) %>%
    ungroup() %>%
    mutate(participant_diff = cond2 - cond1) %>%
    group_by(basis_task, basis_label) %>%
    summarise(
        participant_diff_sd = sd(participant_diff),
        participant_diff_min = min(participant_diff),
        participant_diff_max = max(participant_diff),
        participant_diff_abs_max = max(abs(participant_diff)),
        .groups = "drop"
    )

basis_vector_summary <- lapply(task_ids, function(task_id) {
    basis_index <- selected_basis_indices[task_id]
    basis_vector <- eig_vectors[, basis_index]
    data.frame(
        basis_task = task_id,
        basis_label = sprintf("lap%02d_lambda%.6f", task_id,
                              eig_values[basis_index]),
        basis_index = basis_index,
        eigenvalue = eig_values[basis_index],
        vector_sum = sum(basis_vector),
        vector_mean = mean(basis_vector),
        vector_sd = sd(basis_vector),
        vector_min = min(basis_vector),
        vector_max = max(basis_vector),
        positive_channels = sum(basis_vector > 0),
        negative_channels = sum(basis_vector < 0)
    )
}) %>% bind_rows()

summary <- score_summary %>%
    left_join(condition_summary, by = c("basis_task", "basis_label")) %>%
    left_join(participant_summary, by = c("basis_task", "basis_label")) %>%
    left_join(basis_vector_summary, by = c("basis_task", "basis_label",
                                           "eigenvalue"))

out_dir <- file.path(project_root(), "results", "diagnostics",
                     "laplacian_score")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
tag <- sprintf("%s_seed%03d_k%d_%s",
               data_subdir, seed, k_basis,
               ifelse(include_constant_basis, "constant", "nonconstant"))

summary_file <- file.path(out_dir, sprintf("score_summary_%s.csv", tag))
score_file <- file.path(out_dir, sprintf("score_data_%s.csv", tag))
write.csv(summary, summary_file, row.names = FALSE)
write.csv(score_data, score_file, row.names = FALSE)

cat("Saved:", summary_file, "\n")
cat("Saved:", score_file, "\n")
print(summary)
