script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_path <- sub("^--file=", "", script_arg[[1]])
project_root <- normalizePath(file.path(dirname(script_path), ".."))
Sys.setenv(PROJECT_ROOT = project_root)

source(file.path(project_root, "code", "pipeline", "generate_montecarlo.R"))

params <- simulation_parameters(n_ppt = 20L, scenario = "high")
specs <- radius_specifications(params$sigma_space)
manifest <- paired_scenario_manifest(params, specs[1, ], seeds = 1:5)

stopifnot(identical(manifest$scenario, "high"))
stopifnot(identical(manifest$n_ppt, 20L))
stopifnot(identical(manifest$n_trials, 20))
stopifnot(identical(manifest$signal_mode, "hard_signal"))
stopifnot(identical(manifest$radius_rule, "1sd"))
stopifnot(identical(manifest$replicate_seeds, 1:5))

out_root <- tempfile("paired-radius-test-")
dir.create(out_root)
on.exit(unlink(out_root, recursive = TRUE), add = TRUE)
output_dirs <- setNames(
    file.path(out_root, specs$data_subdir),
    specs$radius_rule
)
for (path in output_dirs) dir.create(path)

fixture <- lapply(specs$radius_rule, function(rule) {
    list(
        sim_data = data.frame(rule = rule),
        ground_truth = list(radius_rule = rule)
    )
})
names(fixture) <- specs$radius_rule
save_paired_replicate(3L, fixture, output_dirs)

for (rule in specs$radius_rule) {
    sim <- readRDS(file.path(output_dirs[[rule]], "sim_003.rds"))
    gt <- readRDS(file.path(output_dirs[[rule]], "gt_003.rds"))
    stopifnot(identical(sim$rule, rule))
    stopifnot(identical(gt$radius_rule, rule))
}

cat("paired radius orchestration tests passed\n")
