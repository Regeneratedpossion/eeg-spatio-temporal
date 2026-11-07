###########################################################
# 00_check_environment.R
#
# Check the R/Python packages and resource files needed by
# the hard-threshold + SEM BGAM pipeline.
#
# LIDO:
#   module load R/4.4.2-gcc141-base
#   Rscript code/00_check_environment.R
###########################################################

# A 41-character "=" divider line, printed before/after each report section.
RULE <- strrep("=", 41)
project_root <- Sys.getenv("PROJECT_ROOT", getwd())

# Print OK/MISSING per item and return the ones that failed.
check_items <- function(items, ok_fn) {
    width <- max(nchar(items))
    ok <- vapply(items, ok_fn, logical(1))
    cat(sprintf("  %-*s %s\n", width, items,
                ifelse(ok, "OK", "MISSING")), sep = "")
    items[!ok]
}

cat(RULE, "\n Environment check\n", RULE, "\n\n", sep = "")

# ---- R packages (checked first; other categories may depend on them) ----

required_r <- c(
    "brms", "cmdstanr", "dplyr", "tidyr", "tidybayes", "posterior",
    "ggplot2", "patchwork", "reticulate", "igraph"
)
cat("--- R packages ---\n")
missing_r <- check_items(
    required_r,
    function(package) requireNamespace(package, quietly = TRUE)
)

required_files <- c(
    "code/utils.R",
    "code/paper_config.R",
    "code/01_generate_data.R",
    "code/03_fit_bgam.R",
    "code/04_compare_methods.R",
    "code/07_fit_paired_radius_bgam.R",
    "code/ground_truth.R",
    "resources/adjacency/adjacency_xproj_yproj_61ch.rds",
    "resources/electrodes/electrode_layout_61ch.csv",
    "resources/brms_meeg/code/eeg_noise.R",
    "resources/brms_meeg/code/erp_template.R",
    "resources/brms_meeg/code/meanpower.txt"
)

# ---- Everything else: one table, one loop ----
# `requires` gates a category on an R package checked above (skipped if
# that package is missing, since its own check would error).

resource_categories <- list(
    list(
        label = "CmdStan", requires = "cmdstanr",
        items = "CmdStan",
        ok_fn = function(item) {
            # cmdstan_path() errors if CmdStan isn't installed.
            path <- tryCatch(cmdstanr::cmdstan_path(), error = function(e) "")
            nzchar(path)
        }
    ),
    list(
        label = "Python modules", requires = "reticulate",
        items = c("mne", "scipy", "numpy"),
        ok_fn = reticulate::py_module_available
    ),
    list(
        label = "project files", requires = NULL,
        items = required_files,
        ok_fn = function(file) file.exists(file.path(project_root, file))
    )
)

missing_groups <- list("R packages" = missing_r)
for (rc in resource_categories) {
    if (!is.null(rc$requires) && rc$requires %in% missing_r) next
    cat(sprintf("\n--- %s ---\n", rc$label))
    missing_groups[[rc$label]] <- check_items(rc$items, rc$ok_fn)
}

# ---- Summary ----
# Item-level detail was already printed above; just say which categories
# still have problems.

cat("\n", RULE, "\n", sep = "")
problem_categories <- names(missing_groups)[lengths(missing_groups) > 0]
if (length(problem_categories)) {
    cat(" Problems found in:", paste(problem_categories, collapse = ", "), "\n")
} else {
    cat(" Environment looks ready.\n")
}
cat(RULE, "\n", sep = "")

if (length(problem_categories)) quit(status = 1)
