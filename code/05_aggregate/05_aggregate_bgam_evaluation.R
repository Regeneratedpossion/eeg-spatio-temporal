###########################################################
# 05_aggregate_bgam_evaluation.R
#
# Summarise per-seed BGAM evaluations (code/04_evaluate/04_evaluate_bgam.R) for one
# results directory, with Monte Carlo standard errors (Section 2.2.5) and the
# two inclusion policies of Section 2.2.5.
#
# Run:
#   Rscript code/05_aggregate/05_aggregate_bgam_evaluation.R <results_dir>
#
# <results_dir> contains seed_XXX/{evaluation/evaluation_main.csv,
# convergence.csv, nuts_summary.csv}. Writes <results_dir>/evaluation_summary.csv
# and <results_dir>/evaluation_per_seed.csv. If seed_XXX/evaluation/
# evaluation_conventional.csv exists (code/03_conventional/03_evaluate_conventional.R), it also
# writes evaluation_conventional_summary.csv and evaluation_comparison.csv.
###########################################################

suppressPackageStartupMessages(library(dplyr))

# Convergence criteria from Section 2.2.5, applied to the evaluated effect g.
converged <- function(x) {
    x$max_rhat_g < 1.01 & x$min_ess_bulk_g > 400 & x$min_ess_tail_g > 400
}

# Mean-type summary with MCSE = SD / sqrt(N).
mean_row <- function(v, measure) {
    v <- v[is.finite(v)]
    n <- length(v)
    data.frame(measure = measure, n_valid = n, estimate = mean(v),
               sd = if (n > 1) sd(v) else NA_real_,
               mcse = if (n > 1) sd(v) / sqrt(n) else NA_real_)
}

# Bias, MAE, RMSE and SD of a timing error (Section 2.2.4), with MCSEs.
timing_rows <- function(e, label) {
    e <- e[is.finite(e)]
    n <- length(e)
    rmse <- sqrt(mean(e^2))
    rbind(
        data.frame(measure = paste0(label, "_bias"), n_valid = n, estimate = mean(e),
                   sd = sd(e), mcse = sd(e) / sqrt(n)),
        data.frame(measure = paste0(label, "_MAE"), n_valid = n, estimate = mean(abs(e)),
                   sd = sd(abs(e)), mcse = sd(abs(e)) / sqrt(n)),
        # Delta-method MCSE for the RMSE.
        data.frame(measure = paste0(label, "_RMSE"), n_valid = n, estimate = rmse,
                   sd = NA_real_, mcse = sd(e^2) / (2 * rmse * sqrt(n))),
        data.frame(measure = paste0(label, "_SD"), n_valid = n, estimate = sd(e),
                   sd = NA_real_, mcse = sd(e) / sqrt(2 * (n - 1)))
    )
}

summarise_policy <- function(d, policy) {
    metric_cols <- c("F1", "IoU", "FPR", "precision", "recall",
                     "F1_250", "IoU_250", "FPR_250", "FPR_cells",
                     "map_cor", "RMSE", "recovery", "leakage",
                     "recovery_250", "leakage_250", "sesoi_global")
    metric_cols <- intersect(metric_cols, names(d))
    out <- rbind(
        do.call(rbind, lapply(metric_cols, function(m) mean_row(d[[m]], m))),
        timing_rows(d$onset_error_ms, "onset_error_ms"),
        timing_rows(d$offset_error_ms, "offset_error_ms"),
        data.frame(measure = "timing_missing_rate", n_valid = nrow(d),
                   estimate = mean(is.na(d$onset_error_ms)), sd = NA_real_,
                   mcse = sqrt(mean(is.na(d$onset_error_ms)) *
                               (1 - mean(is.na(d$onset_error_ms))) / nrow(d)))
    )
    cbind(policy = policy, n_fits = nrow(d), out)
}

aggregate_bgam_evaluation <- function(results_dir) {
    seeds <- list.files(results_dir, pattern = "^seed_[0-9]+$")
    per_seed <- do.call(rbind, lapply(seeds, function(s) {
        f <- file.path(results_dir, s, "evaluation", "evaluation_main.csv")
        cv <- file.path(results_dir, s, "convergence.csv")
        ns <- file.path(results_dir, s, "nuts_summary.csv")
        if (!all(file.exists(c(f, cv, ns)))) return(NULL)
        cbind(seed = as.integer(sub("seed_", "", s)), read.csv(f),
              read.csv(cv)[, c("max_rhat", "min_ess_bulk", "min_ess_tail")],
              read.csv(ns)[, c("divergences", "treedepth_hits")])
    }))
    per_seed$converged <- converged(per_seed)
    write.csv(per_seed, file.path(results_dir, "evaluation_per_seed.csv"), row.names = FALSE)

    summary <- rbind(
        summarise_policy(per_seed, "all_fits"),
        summarise_policy(per_seed[per_seed$converged, ], "converged_only"),
        data.frame(policy = "all_fits", n_fits = nrow(per_seed), measure = "failure_rate",
                   n_valid = nrow(per_seed), estimate = mean(!per_seed$converged),
                   sd = NA_real_,
                   mcse = sqrt(mean(!per_seed$converged) * mean(per_seed$converged) /
                               nrow(per_seed)))
    )
    write.csv(summary, file.path(results_dir, "evaluation_summary.csv"), row.names = FALSE)

    conv <- read_conventional(results_dir, seeds)
    if (!is.null(conv)) compare_with_conventional(results_dir, per_seed, conv)
    invisible(summary)
}

# Conventional-method results at the main alpha, one row per seed and method.
read_conventional <- function(results_dir, seeds) {
    files <- file.path(results_dir, seeds, "evaluation", "evaluation_conventional.csv")
    have <- file.exists(files)
    if (!any(have)) return(NULL)
    do.call(rbind, lapply(which(have), function(i) {
        x <- read.csv(files[i])
        cbind(seed = as.integer(sub("seed_", "", seeds[i])), x[x$main, ])
    }))
}

# Per-method summaries and paired BGAM-minus-method differences (Section 2.2.5).
compare_with_conventional <- function(results_dir, bgam, conv) {
    by_method <- do.call(rbind, lapply(split(conv, conv$method), function(d) {
        cbind(method = d$method[1], summarise_policy(d, "all_datasets"))
    }))
    write.csv(by_method, file.path(results_dir, "evaluation_conventional_summary.csv"),
              row.names = FALSE)

    measures <- c(F1 = "F1", IoU = "IoU", FPR = "FPR",
                  abs_onset_error_ms = "onset_error_ms",
                  abs_offset_error_ms = "offset_error_ms")
    paired <- do.call(rbind, lapply(split(conv, conv$method), function(d) {
        m <- merge(bgam, d, by = "seed", suffixes = c("_bgam", "_method"))
        do.call(rbind, lapply(names(measures), function(lbl) {
            b <- m[[paste0(measures[[lbl]], "_bgam")]]
            o <- m[[paste0(measures[[lbl]], "_method")]]
            if (grepl("^abs_", lbl)) { b <- abs(b); o <- abs(o) }
            diff <- b - o
            ok <- is.finite(diff)
            data.frame(method = d$method[1], measure = lbl, n_paired = sum(ok),
                       bgam = mean(b[ok]), method_value = mean(o[ok]),
                       difference = mean(diff[ok]),
                       mcse_difference = sd(diff[ok]) / sqrt(sum(ok)),
                       mcse_if_unpaired = sqrt(var(b[ok]) / sum(ok) + var(o[ok]) / sum(ok)))
        }))
    }))
    write.csv(paired, file.path(results_dir, "evaluation_comparison.csv"), row.names = FALSE)
    invisible(paired)
}

if (!interactive() && sys.nframe() == 0L) {
    a <- commandArgs(trailingOnly = TRUE)
    if (length(a) != 1L) stop("Usage: Rscript code/05_aggregate/05_aggregate_bgam_evaluation.R <results_dir>")
    s <- aggregate_bgam_evaluation(a[[1]])
    print(s, digits = 3, row.names = FALSE)
}
