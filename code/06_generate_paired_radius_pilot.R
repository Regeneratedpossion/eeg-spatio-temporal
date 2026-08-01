source(file.path(
    Sys.getenv("PROJECT_ROOT", unset = getwd()),
    "code/02_generate_montecarlo.R"
))

run_paired_radius_simulation(
    seeds = 1:5,
    n_ppt = 20L,
    scenario = "high"
)
