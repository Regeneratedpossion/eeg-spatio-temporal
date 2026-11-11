# Summarise performance across statistical thresholds.
# Run from the project root. Reads results/single_node_cond*/seed_*/evaluation.
# Writes to results/fig_data (AGG_OUT overrides).

source(here::here("code/R/aggregate_utils.R"))

res_dir <- here::here("results")
out_dir <- aggregate_out_dir()


# 1. Helper functions ----------------------------------------

read_results <- function(file) {
    tryCatch({
        x <- read.csv(file)

        columns <- c(
            "rule", "baseline_end", "odds",
            "method", "alpha", "F1", "FPR",
            "onset_error_ms", "offset_error_ms"
        )

        x[, intersect(names(x), columns)]
    }, error = function(e) NULL)
}

mae <- function(errors) {
    mean(abs(errors[is.finite(errors)]))
}

n_valid <- function(errors) {
    sum(is.finite(errors))
}

summarise_threshold <- function(scenario, method, threshold, data) {
    data.frame(
        scenario = scenario,
        method = method,
        thr = threshold,
        n = nrow(data),
        F1 = mean(data$F1),
        FPR = mean(data$FPR),
        on_MAE = mae(data$onset_error_ms),
        off_MAE = mae(data$offset_error_ms),
        N_valid = n_valid(data$onset_error_ms)
    )
}


# 2. Summarise each scenario ---------------------------------

out <- list()

for (scenario in c("1", "3", "2", "4")) {

    scenario_dir <- scenario_results_dir(res_dir, scenario)

    # BGAM: posterior odds thresholds.
    bgam_files <- Sys.glob(file.path(
        scenario_dir,
        "seed_*/evaluation/evaluation_sensitivity.csv"
    ))

    bgam <- do.call(
        rbind,
        lapply(bgam_files, read_results)
    )

    bgam <- bgam[
        bgam$rule == "global" &
            abs(bgam$baseline_end - 0.10) < 1e-9,
    ]

    for (odds in sort(unique(bgam$odds))) {
        selected <- bgam[bgam$odds == odds, ]

        out[[length(out) + 1]] <- summarise_threshold(
            scenario, "BGAM", odds, selected
        )
    }

    # Conventional methods: alpha thresholds.
    conventional_files <- Sys.glob(file.path(
        scenario_dir,
        "seed_*/evaluation/evaluation_conventional.csv"
    ))

    conventional <- do.call(
        rbind,
        lapply(conventional_files, read_results)
    )

    for (method in unique(conventional$method)) {
        for (alpha in sort(unique(conventional$alpha))) {

            selected <- conventional[
                conventional$method == method &
                    conventional$alpha == alpha,
            ]

            out[[length(out) + 1]] <- summarise_threshold(
                scenario, method, alpha, selected
            )
        }
    }
}


# 3. Save the summary ----------------------------------------

result <- do.call(rbind, out)

write.csv(
    result,
    file.path(out_dir, "threshold_sweep_20261009.csv"),
    row.names = FALSE
)

cat(
    nrow(result),
    "rows; seeds per cell:",
    range(result$n),
    "\n"
)
