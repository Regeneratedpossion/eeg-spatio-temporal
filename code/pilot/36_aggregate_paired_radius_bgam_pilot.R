suppressPackageStartupMessages(library(dplyr))

source(file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()), "code/utils.R"))

rules <- c("1sd", "half-max", "2sd")
root <- file.path(project_root(), "results", "paired_radius_bgam_pilot")

comparison <- bind_rows(lapply(rules, function(rule) {
    read.csv(file.path(root, rule, "diagnostic_summary.csv"))
})) %>%
    transmute(
        radius = radius_rule,
        jump = jump_weight,
        active_n = n_active,
        EDF = sprintf("space %.2f/%.0f; time %.2f/%.0f",
                      spatial_edf, spatial_edf_upper,
                      temporal_edf, temporal_edf_upper),
        peak_ratio,
        FP,
        FPR,
        recovery,
        leakage,
        spatial_edf_fraction,
        temporal_edf_fraction,
        spatial_edf_near_upper,
        temporal_edf_near_upper
    )
write.csv(comparison, file.path(root, "radius_comparison.csv"),
          row.names = FALSE)

method_metrics <- bind_rows(lapply(rules, function(rule) {
    read.csv(file.path(root, rule, "method_metrics_250ms.csv"))
})) %>%
    group_by(radius_rule) %>%
    mutate(
        FPR_rank = rank(FPR, ties.method = "min", na.last = "keep"),
        F1_rank = rank(-F1, ties.method = "min", na.last = "keep")
    ) %>%
    ungroup()
write.csv(method_metrics, file.path(root, "method_ranking.csv"),
          row.names = FALSE)

orders <- method_metrics %>%
    arrange(radius_rule, FPR_rank, F1_rank, method) %>%
    group_by(radius_rule) %>%
    summarise(
        FPR_order = paste(method[order(FPR_rank, F1_rank)], collapse = " > "),
        F1_order = paste(method[order(F1_rank, FPR_rank)], collapse = " > "),
        .groups = "drop"
    )
consistency <- data.frame(
    metric = c("FPR", "F1"),
    ranking_consistent = c(
        length(unique(orders$FPR_order)) == 1L,
        length(unique(orders$F1_order)) == 1L
    ),
    orders_by_radius = c(
        paste(paste(orders$radius_rule, orders$FPR_order, sep = ": "),
              collapse = " | "),
        paste(paste(orders$radius_rule, orders$F1_order, sep = ": "),
              collapse = " | ")
    )
)
write.csv(consistency, file.path(root, "ranking_consistency.csv"),
          row.names = FALSE)

print(comparison)
print(method_metrics)
print(consistency)
