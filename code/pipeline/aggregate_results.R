###########################################################
# 04_aggregate_results.R
#
# Final aggregation after Monte Carlo jobs finish.
#
# Inputs:
#   results/comparisons/<MC_COMPARE_SUBDIR>/sim_*.csv
#   data/<MC_DATA_SUBDIR>/sim_*.rds
#   data/<MC_DATA_SUBDIR>/gt_*.rds
#   results/<MC_SUMMARY_SUBDIR>/bgam_summary_*.rds
#
# Outputs:
#   <MC_AGG_PREFIX>_all.csv
#   <MC_AGG_PREFIX>_st_f1.csv
#   <MC_AGG_PREFIX>_paper_table.csv
###########################################################

source(file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()), "code/utils.R"))
source(file.path(project_root(), "code/pipeline/compare_methods.R"))

paths <- mc_paths()
dir.create(paths$out_dir, showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------
# Load per-simulation onset/offset results
# ---------------------------------------------------------------

load_clean_results <- function(paths) {
    files <- list.files(paths$compare_dir, pattern = "^sim_\\d+\\.csv$",
                        full.names = TRUE)
    stopifnot(length(files) > 0)

    all_results <- bind_rows(lapply(files, read.csv))
    bad_sims <- all_results %>%
        filter(n_div > 0) %>%
        pull(sim_id) %>%
        unique()

    list(
        all = all_results,
        clean = all_results %>% filter(!(sim_id %in% bad_sims))
    )
}

# ---------------------------------------------------------------
# Spatio-temporal F1 and FP/FN diagnosis
# ---------------------------------------------------------------

mask_metrics <- function(pred_mask, truth_mask) {
    TP <- sum(pred_mask & truth_mask, na.rm = TRUE)
    FP <- sum(pred_mask & !truth_mask, na.rm = TRUE)
    FN <- sum(!pred_mask & truth_mask, na.rm = TRUE)
    TN <- sum(!pred_mask & !truth_mask, na.rm = TRUE)

    data.frame(
        TP = TP,
        FP = FP,
        FN = FN,
        TN = TN,
        precision = TP / (TP + FP),
        recall = TP / (TP + FN),
        F1_ST = 2 * TP / (2 * TP + FP + FN)
    )
}

compute_st_metrics <- function(seed, adj_mat) {
    files <- montecarlo_files(seed, paths)
    sim <- readRDS(files$sim)
    gt <- readRDS(files$gt)
    bg_summary <- readRDS(files$summary)

    bundle <- build_diff_array(aggregate_trials(sim))
    channels <- bundle$channels
    times <- bundle$times
    adj_use <- adj_mat[match(channels, rownames(adj_mat)),
                       match(channels, rownames(adj_mat))]

    truth_channels <- channels %in% gt$active_electrodes
    truth_mask <- outer(times, truth_channels, FUN = function(tt, active)
        (tt >= gt$true_onset & tt <= gt$true_offset) & active)

    masks <- build_method_masks(
        arr = bundle$arr,
        channels = channels,
        times = times,
        adj_chan = adj_use,
        bgam_full = bg_summary$bgam_full,
        mne = mne,
        alpha = ALPHA
    )

    bind_rows(lapply(names(masks), function(method) {
        data.frame(seed = seed, method = method,
                   mask_metrics(masks[[method]], truth_mask))
    }))
}

summarise_st_metrics <- function(st_metrics) {
    st_metrics %>%
        group_by(method) %>%
        summarise(
            mean_TP = mean(TP, na.rm = TRUE),
            mean_FP = mean(FP, na.rm = TRUE),
            mean_FN = mean(FN, na.rm = TRUE),
            FP_rate = sum(FP, na.rm = TRUE) / sum(FP + TN, na.rm = TRUE),
            FN_rate = sum(FN, na.rm = TRUE) / sum(FN + TP, na.rm = TRUE),
            mean_precision = mean(precision, na.rm = TRUE),
            mean_recall = mean(recall, na.rm = TRUE),
            mean_F1_ST = mean(F1_ST, na.rm = TRUE),
            sd_F1_ST = sd(F1_ST, na.rm = TRUE),
            .groups = "drop"
        )
}

# ---------------------------------------------------------------
# Final method table
# ---------------------------------------------------------------

make_final_table <- function(clean_results, st_summary) {
    clean_results %>%
        group_by(method) %>%
        summarise(
            n_sims = n(),
            bias_onset_s = mean(bias_onset_ms, na.rm = TRUE) / 1000,
            mae_onset_s = mean(abs(bias_onset_ms), na.rm = TRUE) / 1000,
            rmse_onset_s = sqrt(mean(bias_onset_ms^2, na.rm = TRUE)) / 1000,
            sd_onset_s = sd(bias_onset_ms, na.rm = TRUE) / 1000,
            bias_offset_s = mean(bias_offset_ms, na.rm = TRUE) / 1000,
            mae_offset_s = mean(abs(bias_offset_ms), na.rm = TRUE) / 1000,
            rmse_offset_s = sqrt(mean(bias_offset_ms^2, na.rm = TRUE)) / 1000,
            sd_offset_s = sd(bias_offset_ms, na.rm = TRUE) / 1000,
            .groups = "drop"
        ) %>%
        left_join(
            st_summary %>%
                mutate(F1_ST = sprintf("%.3f +/- %.3f", mean_F1_ST, sd_F1_ST)) %>%
                select(method, F1_ST, mean_FP, mean_FN, FP_rate, FN_rate),
            by = "method"
        ) %>%
        arrange(mae_onset_s + mae_offset_s)
}

# ---------------------------------------------------------------
# Main aggregation
# ---------------------------------------------------------------

results <- load_clean_results(paths)
all_results <- results$all
clean_results <- results$clean

write.csv(all_results, file.path(paths$out_dir, paste0(paths$prefix, "_all.csv")),
          row.names = FALSE)

available_seeds <- sort(unique(clean_results$sim_id))
adj_mat <- readRDS(adjacency_path("adjacency_xproj_yproj_61ch.rds"))

st_metrics <- bind_rows(lapply(available_seeds, compute_st_metrics,
                               adj_mat = adj_mat))
st_summary <- summarise_st_metrics(st_metrics)

write.csv(st_metrics, file.path(paths$out_dir, paste0(paths$prefix, "_st_f1.csv")),
          row.names = FALSE)

paper_table <- make_final_table(clean_results, st_summary)

print(paper_table, n = Inf, width = Inf)

write.csv(paper_table,
          file.path(paths$out_dir, paste0(paths$prefix, "_paper_table.csv")),
          row.names = FALSE)
