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

RULE <- strrep("=", 41)
project_root <- Sys.getenv("PROJECT_ROOT", unset = getwd())

# Print one OK/MISSING line per item; return the items that failed.
check_items <- function(items, ok_fn) {
    width <- max(nchar(items))
    ok <- vapply(items, function(item) {
        isTRUE(tryCatch(ok_fn(item), error = function(e) FALSE))
    }, logical(1))
    cat(sprintf("  %-*s %s\n", width, items,
                ifelse(ok, "OK", "MISSING")), sep = "")
    items[!ok]
}

cat(RULE, "\n Environment check\n", RULE, "\n\n", sep = "")

# ---- R packages ----
required_r <- c(
    "brms", "cmdstanr", "dplyr", "tidyr", "tidybayes", "posterior",
    "ggplot2", "patchwork", "reticulate", "igraph"
)

cat("--- R packages ---\n")
missing_r <- check_items(
    required_r,
    function(package) requireNamespace(package, quietly = TRUE)
)

# ---- CmdStan ----
missing_cmdstan <- character(0)
if (!("cmdstanr" %in% missing_r)) {
    cat("\n--- CmdStan ---\n")
    cmdstan_path <- tryCatch(
        cmdstanr::cmdstan_path(),
        error = function(e) ""
    )
    cmdstan_ok <- length(cmdstan_path) == 1L && nzchar(cmdstan_path)
    cat("  cmdstan:",
        if (cmdstan_ok) cmdstan_path else "NOT FOUND", "\n")
    if (!cmdstan_ok) missing_cmdstan <- "CmdStan"
}

# ---- Python modules ----
missing_python <- character(0)
if (!("reticulate" %in% missing_r)) {
    cat("\n--- Python modules ---\n")
    cat("  Python:",
        tryCatch(
            reticulate::py_config()$python,
            error = function(e) "NONE"
        ),
        "\n")
    missing_python <- check_items(
        c("mne", "scipy", "numpy"),
        reticulate::py_module_available
    )
}

# ---- Project resources ----
required_files <- c(
    "code/utils.R",
    "code/paper_config.R",
    "code/01_simulate_eeg_data.R",
    "code/02_generate_montecarlo.R",
    "code/03_fit_bgam.R",
    "code/04_compare_methods.R",
    "code/07_fit_paired_radius_bgam.R",
    "code/ground_truth.R",
    "data/sim_data_pilot.rds",
    "resources/adjacency/adjacency_xproj_yproj_61ch.rds",
    "resources/brms_meeg/code/eeg_noise.R",
    "resources/brms_meeg/code/erp_template.R",
    "resources/brms_meeg/code/meanpower.txt"
)
cat("\n--- Project resources ---\n")
missing_files <- check_items(
    required_files,
    function(file) file.exists(file.path(project_root, file))
)

# ---- Summary ----
cat("\n", RULE, "\n", sep = "")
if (length(missing_r)) {
    cat(" Missing R packages:\n",
        paste0("   - ", missing_r, "\n"), sep = "")
}
if (length(missing_cmdstan)) {
    cat(" Missing CmdStan installation.\n")
}
if (length(missing_python)) {
    cat(" Missing Python modules:\n",
        paste0("   - ", missing_python, "\n"), sep = "")
}
if (length(missing_files)) {
    cat(" Missing project files:\n",
        paste0("   - ", missing_files, "\n"), sep = "")
}

ready <- !length(c(
    missing_r,
    missing_cmdstan,
    missing_python,
    missing_files
))
if (ready) cat(" Environment looks ready.\n")
cat(RULE, "\n", sep = "")

if (!ready) quit(status = 1)
