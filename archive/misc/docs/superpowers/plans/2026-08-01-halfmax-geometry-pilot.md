# Half-Max BGAM Geometry Pilot Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build, verify, upload, and submit a three-arm single-chain half-max BGAM geometry pilot that isolates the effects of lower `adapt_delta` and a better-conditioned tensor basis.

**Architecture:** A focused R experiment module defines the immutable arm table, builds one of two brms formulae, fits arms A/B/C sequentially, retains raw sampler artifacts, and writes per-arm plus combined diagnostics. A single Slurm script gives all arms the same node and process environment. A source-level regression test exercises configuration, formula construction, diagnostic summarization, and overwrite protection without running MCMC.

**Tech Stack:** R 4.4, brms, cmdstanr, posterior, dplyr, base R tests, Bash, Slurm, TU Dortmund HPC.

## Global Constraints

- Input is `data/montecarlo_rhalfmax/sim_001.rds` with truth `data/montecarlo_rhalfmax/gt_001.rds`.
- Every arm uses scenario `high`, radius rule `half-max`, 20 participants, seed 1, `k_xy = 16`, and `k_time = 12`.
- Every arm uses `metric = "diag_e"`, one chain, four threads, 1250 iterations, 1000 warmup iterations, and `max_treedepth = 10`.
- Arm A is `t2(ps, ps, full = FALSE)` with `adapt_delta = 0.80`.
- Arm B is `t2(ps, ps, full = FALSE)` with `adapt_delta = 0.75`.
- Arm C is `te(ps, ps, np = FALSE)` with `adapt_delta = 0.80`.
- Arms run sequentially in one Slurm allocation; existing paired-radius outputs and the running `dense_halfmax` job are read-only.
- The pilot is descriptive. It does not automatically submit a combined arm or choose a final production configuration.

---

### Task 1: Arm specification, formula, and diagnostic helpers

**Files:**
- Create: `experiments/geometry_pilot_halfmax/run_geometry_pilot_halfmax.R`
- Create: `tests/test_geometry_pilot_halfmax.R`

**Interfaces:**
- Consumes: `aggregate_trials()` from `code/utils.R`; half-max simulation and truth RDS files.
- Produces: `geometry_pilot_arms() -> data.frame`, `geometry_pilot_formula(arm_id) -> brmsformula`, `summarize_geometry_nuts(nuts, elapsed) -> data.frame`, and `geometry_pilot_paths(project_dir) -> list`.

- [ ] **Step 1: Write the failing configuration and formula test**

```r
project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = getwd()))
Sys.setenv(GEOMETRY_PILOT_TESTING = "true")
source(file.path(
    project_root,
    "experiments/geometry_pilot_halfmax/run_geometry_pilot_halfmax.R"
))

arms <- geometry_pilot_arms()
stopifnot(
    identical(arms$arm, c("A", "B", "C")),
    identical(arms$smooth, c("t2", "t2", "te")),
    identical(arms$adapt_delta, c(0.80, 0.75, 0.80)),
    identical(arms$metric, rep("diag_e", 3)),
    identical(arms$iter, rep(1250L, 3)),
    identical(arms$warmup, rep(1000L, 3))
)

formula_a <- paste(deparse(geometry_pilot_formula("A")$formula), collapse = " ")
formula_b <- paste(deparse(geometry_pilot_formula("B")$formula), collapse = " ")
formula_c <- paste(deparse(geometry_pilot_formula("C")$formula), collapse = " ")
stopifnot(
    grepl("t2\\(", formula_a),
    identical(formula_a, formula_b),
    grepl("full = FALSE", formula_a, fixed = TRUE),
    grepl("te\\(", formula_c),
    grepl("np = FALSE", formula_c, fixed = TRUE)
)
```

- [ ] **Step 2: Run the test and verify the missing experiment module causes failure**

Run:

```bash
env PROJECT_ROOT="$PWD" Rscript --vanilla tests/test_geometry_pilot_halfmax.R
```

Expected: FAIL because `run_geometry_pilot_halfmax.R` or `geometry_pilot_arms()` does not exist.

- [ ] **Step 3: Implement immutable paths, arms, and formula construction**

The module must source `code/utils.R` and expose:

```r
geometry_pilot_paths <- function(project_dir = project_root()) {
    experiment_dir <- file.path(project_dir, "experiments", "geometry_pilot_halfmax")
    list(
        project_dir = project_dir,
        experiment_dir = experiment_dir,
        results_dir = file.path(experiment_dir, "results"),
        input_file = file.path(project_dir, "data", "montecarlo_rhalfmax", "sim_001.rds"),
        truth_file = file.path(project_dir, "data", "montecarlo_rhalfmax", "gt_001.rds")
    )
}

geometry_pilot_arms <- function() {
    data.frame(
        arm = c("A", "B", "C"),
        smooth = c("t2", "t2", "te"),
        adapt_delta = c(0.80, 0.75, 0.80),
        metric = rep("diag_e", 3),
        iter = rep(1250L, 3),
        warmup = rep(1000L, 3),
        stringsAsFactors = FALSE
    )
}
```

`geometry_pilot_formula()` must use the common likelihood and random-effects structure:

```r
eeg_mean | se(eeg_sem, sigma = TRUE) ~
    condition +
    (1 + condition | participant) +
    SMOOTH
```

For A/B, `SMOOTH` is the approved `t2(..., full = FALSE)` term. For C it is the approved `te(..., np = FALSE)` term. Both use `by = condition`, `bs = c("ps", "ps")`, `k = c(16, 12)`, and `d = c(2, 1)`.

- [ ] **Step 4: Add a failing diagnostic-summary fixture**

Append a synthetic `nuts_params` fixture with two divergences, known leapfrog quantiles, and one saturated treedepth. Assert that `summarize_geometry_nuts()` reports those exact values and uses the supplied warmup/sample elapsed times.

- [ ] **Step 5: Implement `summarize_geometry_nuts()` minimally**

The function must subset `divergent__`, `treedepth__`, `stepsize__`, `n_leapfrog__`, and `accept_stat__`, then return one row with exact fields:

```r
c(
    divergences,
    treedepth_saturated,
    treedepth_saturated_rate,
    stepsize,
    mean_leapfrog,
    median_leapfrog,
    q90_leapfrog,
    max_leapfrog,
    mean_accept_stat,
    warmup_seconds,
    sampling_seconds,
    sampling_seconds_per_iteration
)
```

- [ ] **Step 6: Run the focused test**

Run:

```bash
env PROJECT_ROOT="$PWD" Rscript --vanilla tests/test_geometry_pilot_halfmax.R
```

Expected: PASS with `half-max geometry pilot helper tests passed`.

- [ ] **Step 7: Commit the tested helpers**

```bash
git add experiments/geometry_pilot_halfmax/run_geometry_pilot_halfmax.R tests/test_geometry_pilot_halfmax.R
git commit -m "feat: define halfmax geometry pilot arms"
```

---

### Task 2: Sequential fitting, artifacts, and comparison output

**Files:**
- Modify: `experiments/geometry_pilot_halfmax/run_geometry_pilot_halfmax.R`
- Modify: `tests/test_geometry_pilot_halfmax.R`
- Create: `experiments/geometry_pilot_halfmax/submit_geometry_pilot_halfmax.sh`

**Interfaces:**
- Consumes: helpers from Task 1 and brms/cmdstanr model objects.
- Produces: `run_geometry_arm()`, `copy_sampler_artifacts()`, `summarize_geometry_predictions()`, `run_geometry_pilot()`, and a Slurm entry point.

- [ ] **Step 1: Add failing orchestration assertions**

Add source-level assertions that the experiment module contains:

```r
for (required in c(
    "chains = 1", "iter = arm$iter", "warmup = arm$warmup",
    "threads = brms::threading(threads = 4)",
    "metric = arm$metric", "save_warmup = TRUE", "refresh = 25",
    "max_treedepth = 10", "seed = 1"
)) stopifnot(grepl(required, source_text, fixed = TRUE))
```

Also assert that arm output paths are distinct and live only under `experiments/geometry_pilot_halfmax/results/`.

- [ ] **Step 2: Run the focused test and verify it fails on missing orchestration**

Run the Task 1 test command. Expected: FAIL on the first missing required fitting argument.

- [ ] **Step 3: Implement input validation and overwrite protection**

Before compiling, read both RDS files and require:

```r
stopifnot(
    identical(truth$scenario, "high"),
    identical(truth$radius_rule, "half-max"),
    length(unique(sim_data$participant)) == 20L
)
```

Each arm receives `results/<arm>/`. If `fit.rds` or `complete.flag` exists, stop rather than overwrite. Temporary and compiled model files use experiment-specific directories derived from `SLURM_JOB_ID`.

- [ ] **Step 4: Implement sequential arm execution**

For each arm, aggregate trials once, apply the existing `c(-0.5, 0.5)` condition contrast, build the arm formula, compile before starting the sampling timer, and fit using the common arguments. Use `backend = "cmdstanr"`, `metric = "diag_e"`, `chains = 1`, `cores = 1`, and four within-chain threads.

Retain the brms default prior families and save `prior_summary(fit)` for every arm so the report can verify the same prior family is applied to common parameter classes.

- [ ] **Step 5: Preserve raw sampler artifacts**

`copy_sampler_artifacts()` must support both a runtime `CmdStanMCMC` object with `$output_files()` and a brms-wrapped `stanfit` with `@stan_args`. Copy CSV, metric JSON, configuration JSON, and profile files into `results/<arm>/sampler/` before the R session exits. Stop if no sampler CSV is found.

- [ ] **Step 6: Generate per-arm diagnostics and fixed-grid predictions**

Write:

```text
results/<arm>/fit.rds
results/<arm>/config.rds
results/<arm>/prior_summary.csv
results/<arm>/nuts_summary.csv
results/<arm>/posterior_mean_grid.rds
results/<arm>/complete.flag
```

The prediction grid contains every channel and time for both conditions with population-level prediction (`re_formula = NA`). Store the cond2-minus-cond1 posterior mean and report CZ at 250 ms.

- [ ] **Step 7: Generate the combined comparison only after A/B/C complete**

Write `results/geometry_comparison.csv` with the NUTS fields plus CZ@250 ms. Join each arm's grid to Arm A and add posterior-mean-grid RMSE and maximum absolute difference. Do not submit or generate Arm D.

- [ ] **Step 8: Implement the Slurm wrapper**

Use one task, four CPUs, 64 GB memory, the `ultralong` partition, and a seven-day limit. Load `R/4.4.2-gcc141-base`, set project-specific `TMPDIR`, `CMDSTAN_CACHE_DIR`, and `R_LIBS_USER`, then call the R module once. Enable `END,FAIL` email notification using the project's existing configured address.

- [ ] **Step 9: Run local static and regression verification**

Run:

```bash
env PROJECT_ROOT="$PWD" Rscript --vanilla tests/test_geometry_pilot_halfmax.R
Rscript -e 'parse(file="experiments/geometry_pilot_halfmax/run_geometry_pilot_halfmax.R"); cat("R parse passed\n")'
bash -n experiments/geometry_pilot_halfmax/submit_geometry_pilot_halfmax.sh
git diff --check -- experiments/geometry_pilot_halfmax tests/test_geometry_pilot_halfmax.R
```

Expected: all commands exit 0.

- [ ] **Step 10: Run the full local test suite**

Run every `tests/test_*.R` file under the project with `PROJECT_ROOT="$PWD"`. Expected: all exit 0.

- [ ] **Step 11: Commit orchestration**

```bash
git add experiments/geometry_pilot_halfmax tests/test_geometry_pilot_halfmax.R
git commit -m "feat: add halfmax geometry pilot runner"
```

---

### Task 3: Upload, remote verification, and Slurm submission

**Files:**
- Upload: `experiments/geometry_pilot_halfmax/run_geometry_pilot_halfmax.R`
- Upload: `experiments/geometry_pilot_halfmax/submit_geometry_pilot_halfmax.sh`
- Upload: `tests/test_geometry_pilot_halfmax.R`

**Interfaces:**
- Consumes: verified local experiment files and the existing TU Dortmund project.
- Produces: one submitted Slurm job ID and a verified initial queue/log state.

- [ ] **Step 1: Verify the existing TU Dortmund job is untouched**

Run `squeue -u "$USER"` through `ssh gw01` and record the current `dense_halfmax` job ID and state. Do not cancel, requeue, or modify it.

- [ ] **Step 2: Upload only the pilot files**

Use `rsync -avR` from the local project root to `/work/smjifan1/bgam_eeg_project/`. Do not synchronize the whole dirty worktree.

- [ ] **Step 3: Run remote syntax and input checks**

On `gw01`, load the same R module and library path as the Slurm script, then run the focused test, parse the R module, run `bash -n`, and verify both half-max RDS inputs exist. Expected: all exit 0.

- [ ] **Step 4: Submit exactly one sequential pilot job**

From the remote project root run:

```bash
sbatch experiments/geometry_pilot_halfmax/submit_geometry_pilot_halfmax.sh
```

Capture the returned job ID.

- [ ] **Step 5: Verify queue and startup metadata**

Run `squeue` and `scontrol show job` for the new ID. Confirm one node, one task, four CPUs, 64 GB memory, seven-day time limit, and the expected command/work directory. Check that the existing `dense_halfmax` job remains present or has ended naturally rather than being modified by this work.

- [ ] **Step 6: Report submission state**

Report the new job ID, partition/state/reason, output log path, the unchanged status of `dense_halfmax`, and the exact commands for subsequent progress checks. Do not claim sampling has begun unless the log contains `Start sampling`.

