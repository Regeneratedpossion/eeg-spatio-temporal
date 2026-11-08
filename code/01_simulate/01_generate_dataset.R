# Generates one dataset of the thesis: scenario COND (I, H, F or E), seed from SLURM_ARRAY_TASK_ID or argument 1.
# Usage: COND=E Rscript code/01_simulate/01_generate_dataset.R 1
# Writes data/single_node_cond<COND>_seeds001_100/sim_<seed>.rds and gt_<seed>.rds (DATA_OUT overrides the folder).
suppressPackageStartupMessages(library(dplyr))
source(here::here("code/01_simulate/01_generate_data.R"))

cond <- Sys.getenv("COND")
if (!cond %in% c("I", "H", "F", "E")) stop("COND must be I, H, F or E", call. = FALSE)
seed <- as.integer(c(commandArgs(trailingOnly = TRUE), Sys.getenv("SLURM_ARRAY_TASK_ID"))[1])
if (is.na(seed) || !seed %in% 1:100) stop("seed must be an integer from 1 to 100", call. = FALSE)

out_dir <- Sys.getenv("DATA_OUT", unset = data_path(sprintf("single_node_cond%s_seeds001_100", cond)))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
sim_file <- file.path(out_dir, sprintf("sim_%03d.rds", seed))
if (file.exists(sim_file) || file.exists(file.path(out_dir, sprintf("gt_%03d.rds", seed)))) stop("Output already exists for seed ", seed, call. = FALSE)

d <- make_dataset(seed = seed, scenario = paste0("cond_", cond))
saveRDS(d$sim_data, sim_file)
saveRDS(d$ground_truth, file.path(out_dir, sprintf("gt_%03d.rds", seed)))
cat("generated cond", cond, "seed", seed, "\n")
