# Generate one EEG dataset for scenario COND (1, 2, 3 or 4).
# Usage: COND=4 Rscript code/01_simulate/01_generate_dataset.R 1
# Seed comes from argument 1 or SLURM_ARRAY_TASK_ID.
# Writes data/single_node_cond<COND>_seeds001_100/sim_<seed>.rds and gt_<seed>.rds.
# DATA_OUT overrides the default output directory.

source(here::here("code/01_simulate/01_generate_data.R"))

# 1. Read scenario and seed ----------------------------------

cond <- Sys.getenv("COND")

if (!cond %in% c("1", "2", "3", "4")) {
    stop("COND must be 1, 2, 3 or 4", call. = FALSE)
}

args <- commandArgs(trailingOnly = TRUE)
seed_input <- if (length(args)) {
    args[1]
} else {
    Sys.getenv("SLURM_ARRAY_TASK_ID")
}

seed <- suppressWarnings(as.numeric(seed_input))

if (length(seed) != 1L || !is.finite(seed) ||
    seed != floor(seed) || !seed %in% 1:100) {
    stop("seed must be an integer from 1 to 100", call. = FALSE)
}

seed <- as.integer(seed)


# 2. Prepare output files ------------------------------------

out_dir <- Sys.getenv(
    "DATA_OUT",
    unset = data_path(
        sprintf("single_node_cond%s_seeds001_100", cond)
    )
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sim_file <- file.path(out_dir, sprintf("sim_%03d.rds", seed))
gt_file <- file.path(out_dir, sprintf("gt_%03d.rds", seed))

if (file.exists(sim_file) || file.exists(gt_file)) {
    stop("Output already exists for seed ", seed, call. = FALSE)
}


# 3. Generate and save ---------------------------------------

d <- make_dataset(
    seed = seed,
    scenario = paste0("cond_", cond)
)

saveRDS(d$sim_data, sim_file)
saveRDS(d$ground_truth, gt_file)

cat("generated cond", cond, "seed", seed, "\n")
