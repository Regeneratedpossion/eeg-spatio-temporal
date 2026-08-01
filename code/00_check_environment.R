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

cat("=========================================\n")
cat(" Environment check\n")
cat("=========================================\n\n")

project_root <- Sys.getenv("PROJECT_ROOT", unset = getwd())

# ---- R packages ----
required_r <- c(
    "brms", "cmdstanr", "dplyr", "tidyr", "tidybayes", "posterior",
    "ggplot2", "patchwork", "reticulate", "igraph"
)

cat("--- R packages ---\n")
missing_r <- character(0)
for (pkg in required_r) {
    ok <- requireNamespace(pkg, quietly = TRUE)
    cat(sprintf("  %-12s %s\n", pkg, if (ok) "OK" else "MISSING"))
    if (!ok) missing_r <- c(missing_r, pkg)
}

# ---- CmdStan ----
cat("\n--- CmdStan ---\n")
if (requireNamespace("cmdstanr", quietly = TRUE)) {
    cmdstan_path <- tryCatch(cmdstanr::cmdstan_path(), error = function(e) NA)
    cat(sprintf("  cmdstan path: %s\n",
                if (is.na(cmdstan_path)) "NOT FOUND" else cmdstan_path))
}

# ---- Python modules ----
cat("\n--- Python modules ---\n")
if (requireNamespace("reticulate", quietly = TRUE)) {
    suppressPackageStartupMessages(library(reticulate))
    cat(sprintf("  Python: %s\n",
                tryCatch(py_config()$python, error = function(e) "NONE")))
    for (mod in c("mne", "scipy", "numpy")) {
        ok <- py_module_available(mod)
        cat(sprintf("  %-12s %s\n", mod, if (ok) "OK" else "MISSING"))
    }
}

# ---- Project resources ----
cat("\n--- Project resources ---\n")
required_files <- c(
    "code/utils.R",
    "code/01_simulate_eeg_data.R",
    "code/03_fit_bgam.R",
    "code/04_compare_methods.R",
    "resources/adjacency/adjacency_xproj_yproj_61ch.rds",
    "resources/brms_meeg/code/eeg_noise.R",
    "resources/brms_meeg/code/erp_template.R",
    "resources/brms_meeg/code/meanpower.txt"
)

missing_files <- character(0)
for (file in required_files) {
    ok <- file.exists(file.path(project_root, file))
    cat(sprintf("  %-55s %s\n", file, if (ok) "OK" else "MISSING"))
    if (!ok) missing_files <- c(missing_files, file)
}

cat("\n=========================================\n")
if (length(missing_r) == 0 && length(missing_files) == 0) {
    cat(" Environment looks ready.\n")
} else {
    if (length(missing_r) > 0) {
        cat(" Missing R packages:\n")
        cat(paste0("   - ", missing_r, "\n"))
    }
    if (length(missing_files) > 0) {
        cat(" Missing project files:\n")
        cat(paste0("   - ", missing_files, "\n"))
    }
}
cat("=========================================\n")
