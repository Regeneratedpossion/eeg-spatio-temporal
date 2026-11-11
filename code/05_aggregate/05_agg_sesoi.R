# Compare the two SESOI rules (delta = 0 and the global baseline rule).
# Usage: Rscript code/05_aggregate/05_agg_sesoi.R [all|passed]
# all (default): every fit. passed: only the fits in resources/passed_fits.txt.
# Run from the project root. Reads results/single_node_cond*/seed_*/evaluation.
# Writes results/fig_data/tables/sesoi_rule_20261009.csv (all) or
# sesoi_rule_passed.csv (passed); AGG_OUT overrides the output folder.

mode <- commandArgs(TRUE)[1]
if (is.na(mode)) mode <- "all"
stopifnot(mode %in% c("all", "passed"))

res_dir <- here::here("results")
out_dir <- Sys.getenv(
    "AGG_OUT",
    unset = here::here("results", "fig_data")
)

dir.create(
    file.path(out_dir, "tables"),
    showWarnings = FALSE,
    recursive = TRUE
)

if (mode == "passed") {
    man <- read.table(
        here::here("resources/passed_fits.txt"),
        col.names = c("cond", "seed")
    )
}


# 1. Read the sensitivity results of one file ----------------

read_results <- function(file) {
    tryCatch({
        x <- read.csv(file)

        if (!("post_offset" %in% names(x))) return(NULL)

        x[x$odds == 10, c(
            "rule", "baseline_end", "F1", "FPR", "recall", "post_offset"
        )]
    }, error = function(e) NULL)
}


# 2. Summarise each scenario ---------------------------------

out <- list()

for (scenario in c("1", "3", "2", "4")) {

    scenario_dir <- file.path(
        res_dir,
        sprintf(
            "single_node_cond%s_ds_cr_seeds001_100",
            scenario
        )
    )

    if (mode == "passed") {
        files <- file.path(
            scenario_dir,
            sprintf("seed_%03d", man$seed[man$cond == scenario]),
            "evaluation/evaluation_sensitivity.csv"
        )
    } else {
        files <- Sys.glob(file.path(
            scenario_dir,
            "seed_*/evaluation/evaluation_sensitivity.csv"
        ))
    }

    results <- do.call(rbind, lapply(files, read_results))

    for (rule in c("zero", "global")) {
        selected <- results[
            results$rule == rule &
                (rule == "zero" |
                    abs(results$baseline_end - 0.10) < 1e-9),
        ]

        out[[length(out) + 1]] <- data.frame(
            scenario = scenario,
            rule = rule,
            n = nrow(selected),
            F1 = mean(selected$F1),
            FPR = mean(selected$FPR),
            post_offset = mean(selected$post_offset),
            recall = mean(selected$recall)
        )
    }
}


# 3. Save the summary ----------------------------------------

result <- do.call(rbind, out)

write.csv(
    result,
    file.path(
        out_dir, "tables",
        if (mode == "passed") "sesoi_rule_passed.csv"
        else "sesoi_rule_20261009.csv"
    ),
    row.names = FALSE
)

print(result)
