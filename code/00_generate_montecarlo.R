###########################################################
# 00_generate_montecarlo.R
#
# Generate Monte Carlo datasets for later BGAM fitting and method comparison.
#
# Run:
#   Rscript code/00_generate_montecarlo.R
#
# Optional setting:
#   MC_DATA_SUBDIR=...   output folder under data/
###########################################################

source(file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()), "code/utils.R"))

# Load make_dataset().
source(file.path(project_root(), "code/00_generate_data.R"))

# Generate and save one simulation and its ground truth.
generate_dataset_for_seed <- function(seed, out_dir) {
    sim_path <- file.path(out_dir, sprintf("sim_%03d.rds", seed))
    gt_path <- file.path(out_dir, sprintf("gt_%03d.rds", seed))

    dataset <- make_dataset(seed = seed, n_ppt = 20L)
    saveRDS(dataset$sim_data, sim_path)
    saveRDS(dataset$ground_truth, gt_path)
}

main <- function() {
    n_sim <- 30L
    out_dir <- mc_paths()$data_dir
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

    invisible(lapply(seq_len(n_sim), generate_dataset_for_seed, out_dir = out_dir))
}

main()
