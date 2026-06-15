###########################################################
# check_env.R — verify all dependencies are available
#
# Run this on the cluster before launching any analysis script.
#
# Usage on LIDO:
#   module load R/4.4.2-gcc141-base
#   cd bgam_eeg_project
#   Rscript code/check_env.R
###########################################################

cat("=========================================\n")
cat(" Environment check\n")
cat("=========================================\n\n")

# ---- R packages ----
pkgs <- c(
    # core stats / Bayesian
    "brms", "tidybayes", "posterior", "MASS",
    # data wrangling
    "dplyr", "tidyr", "tidyverse", "purrr",
    # plotting
    "ggplot2", "patchwork", "scales",
    # EEG-specific
    "eegkit", "neurogam",
    # Python bridge
    "reticulate"
)

cat("--- R packages ---\n")
missing_r <- character(0)
for (p in pkgs) {
    ok <- requireNamespace(p, quietly = TRUE)
    cat(sprintf("  %-15s %s\n", p, if (ok) "OK" else "MISSING"))
    if (!ok) missing_r <- c(missing_r, p)
}

# ---- cmdstanr (often installed separately) ----
cat("\n--- cmdstanr backend ---\n")
ok_cmdstanr <- requireNamespace("cmdstanr", quietly = TRUE)
cat(sprintf("  %-15s %s\n", "cmdstanr",
            if (ok_cmdstanr) "OK" else "MISSING"))
if (!ok_cmdstanr) missing_r <- c(missing_r, "cmdstanr")
if (ok_cmdstanr) {
    cs_path <- tryCatch(cmdstanr::cmdstan_path(), error = function(e) NA)
    cat(sprintf("  cmdstan path  : %s\n",
                if (is.na(cs_path)) "NOT FOUND" else cs_path))
}

# ---- Python modules (via reticulate) ----
cat("\n--- Python modules ---\n")
if (!"reticulate" %in% rownames(installed.packages())) {
    cat("  (reticulate not installed; skipping Python check)\n")
} else {
    suppressPackageStartupMessages(library(reticulate))
    cat(sprintf("  Python found  : %s\n",
                tryCatch(py_config()$python, error = function(e) "NONE")))
    py_modules <- c("mne", "scipy", "numpy")
    for (m in py_modules) {
        ok <- py_module_available(m)
        cat(sprintf("  %-15s %s\n", m, if (ok) "OK" else "MISSING"))
    }
}

# ---- Project files ----
cat("\n--- Project files ---\n")
project_root <- Sys.getenv("PROJECT_ROOT", unset = getwd())
required_files <- c(
    "data/sim_data_pilot.rds",
    "data/ground_truth_pilot.rds",
    "results/fits/bgam_fit_v4_k12.rds",
    "code/utils.R"
)
for (f in required_files) {
    p  <- file.path(project_root, f)
    ok <- file.exists(p)
    cat(sprintf("  %-44s %s\n", f, if (ok) "OK" else "MISSING"))
}

# ---- External: brms_meeg repo (for 00_generate_data.R) ----
cat("\n--- External brms_meeg repo (eeg_noise.R, erp_template.R) ---\n")
brms_meeg <- Sys.getenv("BRMS_MEEG", unset = "")
if (!nzchar(brms_meeg)) {
    # try common locations
    cands <- c(
        file.path(project_root, "brms_meeg/code"),
        file.path(dirname(project_root), "brms_meeg/code"),
        "/work/smjifan1/brms_meeg/code"
    )
    for (c in cands) if (dir.exists(c)) { brms_meeg <- c; break }
}
if (nzchar(brms_meeg) && dir.exists(brms_meeg)) {
    cat(sprintf("  Path           : %s\n", brms_meeg))
    for (f in c("eeg_noise.R", "erp_template.R", "utils.R", "meanpower.txt")) {
        ok <- file.exists(file.path(brms_meeg, f))
        cat(sprintf("  %-15s %s\n", f, if (ok) "OK" else "MISSING"))
    }
} else {
    cat("  brms_meeg path NOT FOUND\n")
    cat("  → set BRMS_MEEG env var, or place brms_meeg/ next to project\n")
}

# ---- Summary ----
cat("\n=========================================\n")
if (length(missing_r) == 0) {
    cat(" All R packages OK\n")
} else {
    cat(" Missing R packages:\n")
    cat(paste0("   - ", missing_r, "\n"))
    cat("\n Install with (run inside R):\n")
    cat(sprintf('   install.packages(c(%s),\n',
                paste0('"', missing_r, '"', collapse = ", ")))
    cat('                    repos = "https://cloud.r-project.org")\n')
}
cat("=========================================\n")
