# DS Threading Throughput Seeds 2--10 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Prepare, test, upload, and preflight a paired full-fit comparison of no within-chain threading versus four within-chain threads for half-max seeds 2--10, then stop for user approval before submitting any Slurm job.

**Architecture:** Generate nine isolated half-max datasets once. Run two Slurm arrays over the same seeds: Arm `off` uses four parallel chains and four CPUs with no `threads` argument; Arm `threads4` uses four parallel chains, four threads per chain, and sixteen CPUs. Store each arm-seed independently, aggregate only completed pairs, and normalize performance by both wall time and CPU-hours.

**Tech Stack:** R 4.4.2, brms, cmdstanr/CmdStan 2.38.0, posterior, Slurm arrays, Bash.

## Global Constraints

- Seeds are exactly integers 2 through 10; seed 1 is not regenerated or refitted.
- Scenario is `high`, participant count is 20, radius is `half-max`, and predictors remain raw/unstandardized.
- Formula is `ds x ps`, `k = c(16, 12)`, `d = c(2, 1)`, `m = list(c(1, 0.5), c(2, 2))`, and `full = FALSE`.
- Both arms use four chains, 3,000 iterations, 1,000 warmup, brms default priors, `adapt_delta = 0.8`, `max_treedepth = 10`, and `diag_e`.
- Arm `off` omits `threads` and requests four CPUs; Arm `threads4` uses `threading(threads = 4)` and requests sixteen CPUs.
- Both arrays use `ultralong` and CPU feature `xeon_e52640v4`.
- Existing data, fits, manifests, logs, and the dirty worktree are read-only; implementation creates only new experiment and test files.
- No `sbatch` command may run until the user has reviewed the implemented files, resource budget, and exact submission commands and explicitly approves submission.

---

## File Structure

- Create `experiments/threading_compare_ds_full/threading_compare_helpers.R`: immutable configuration, formula, arm mapping, paths, input validation, and diagnostic summaries.
- Create `experiments/threading_compare_ds_full/generate_inputs.R`: generate and validate isolated half-max seeds 2--10.
- Create `experiments/threading_compare_ds_full/fit_seed.R`: run one full arm-seed fit selected by environment variables.
- Create `experiments/threading_compare_ds_full/aggregate_results.R`: join completed arm pairs and compute wall-time, CPU-hour, throughput, and ESS-efficiency ratios.
- Create `experiments/threading_compare_ds_full/submit_generate.sh`: one-node input generation job.
- Create `experiments/threading_compare_ds_full/submit_off.sh`: seeds 2--10 array, four CPUs per task, no within-chain threading.
- Create `experiments/threading_compare_ds_full/submit_threads4.sh`: seeds 2--10 array, sixteen CPUs per task, four threads per chain.
- Create `experiments/threading_compare_ds_full/submit_aggregate.sh`: one-CPU aggregation job dependent on both arrays.
- Create `tests/test_threading_compare_ds_full.R`: source-level regression tests for scientific invariants, seed mapping, resource mapping, paths, and overwrite guards.

---

### Task 1: Immutable experiment configuration and formula

**Files:**
- Create: `tests/test_threading_compare_ds_full.R`
- Create: `experiments/threading_compare_ds_full/threading_compare_helpers.R`

**Interfaces:**
- Produces: `threading_compare_config(project_root) -> list`
- Produces: `threading_compare_arm(arm) -> list`
- Produces: `threading_compare_formula() -> brmsformula`
- Produces: `threading_compare_paths(config, arm, seed) -> list`
- Produces: `validate_threading_seed(seed) -> integer`

- [ ] **Step 1: Write the failing configuration test**

Create the test bootstrap and assertions:

```r
script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
root <- normalizePath(file.path(
    dirname(sub("^--file=", "", script_arg[[1]])), ".."
))
Sys.setenv(PROJECT_ROOT = root, THREADING_COMPARE_TESTING = "true")
helper <- file.path(
    root, "experiments/threading_compare_ds_full/threading_compare_helpers.R"
)
stopifnot(file.exists(helper))
source(helper)

cfg <- threading_compare_config(root)
stopifnot(
    identical(cfg$seeds, 2:10),
    identical(cfg$iter, 3000L),
    identical(cfg$warmup, 1000L),
    identical(cfg$chains, 4L),
    identical(cfg$cores, 4L),
    identical(cfg$adapt_delta, 0.8),
    identical(cfg$max_treedepth, 10L),
    identical(cfg$metric, "diag_e")
)

off <- threading_compare_arm("off")
threads4 <- threading_compare_arm("threads4")
stopifnot(
    identical(off$threads, 0L), identical(off$cpus, 4L),
    identical(threads4$threads, 4L), identical(threads4$cpus, 16L)
)
stopifnot(inherits(try(threading_compare_arm("bad"), silent = TRUE),
                     "try-error"))
stopifnot(identical(validate_threading_seed("2"), 2L))
stopifnot(inherits(try(validate_threading_seed("1"), silent = TRUE),
                     "try-error"))
```

- [ ] **Step 2: Run the test and verify RED**

Run:

```bash
Rscript --vanilla tests/test_threading_compare_ds_full.R
```

Expected: failure because `threading_compare_helpers.R` does not exist.

- [ ] **Step 3: Implement the immutable helpers**

Implement the exact arm and model configuration:

```r
threading_compare_config <- function(project_root) {
    list(
        project_root = normalizePath(project_root),
        seeds = 2:10,
        data_dir = file.path(
            project_root, "data/threading_compare_halfmax_high_p20"
        ),
        results_dir = file.path(
            project_root, "results/threading_compare_ds_full"
        ),
        chains = 4L, cores = 4L, iter = 3000L, warmup = 1000L,
        adapt_delta = 0.8, max_treedepth = 10L, metric = "diag_e"
    )
}

threading_compare_arm <- function(arm) {
    arms <- list(
        off = list(name = "off", threads = 0L, cpus = 4L),
        threads4 = list(name = "threads4", threads = 4L, cpus = 16L)
    )
    if (length(arm) != 1L || is.null(arms[[arm]])) {
        stop("arm must be off or threads4", call. = FALSE)
    }
    arms[[arm]]
}

validate_threading_seed <- function(seed) {
    seed <- suppressWarnings(as.integer(seed))
    if (length(seed) != 1L || is.na(seed) || !(seed %in% 2:10)) {
        stop("seed must be one integer in 2:10", call. = FALSE)
    }
    seed
}

threading_compare_formula <- function() {
    brms::bf(
        eeg_mean | se(eeg_sem, sigma = TRUE) ~
            condition + (1 + condition | participant) +
            t2(
                xproj, yproj, time, by = condition,
                bs = c("ds", "ps"), k = c(16, 12), d = c(2, 1),
                m = list(c(1, 0.5), c(2, 2)), full = FALSE
            )
    )
}
```

`threading_compare_paths()` must format seeds with three digits and return
arm-specific `fit.rds`, `config.rds`, `runtime.csv`, `nuts_summary.csv`,
`convergence.csv`, `sampler/`, `cmdstan_cache/`, and `complete.flag` paths.

- [ ] **Step 4: Add exact formula and path assertions**

Append assertions that normalize whitespace in `deparse()` and require these
fragments:

```r
formula_text <- gsub(
    "[[:space:]]+", " ",
    paste(deparse(threading_compare_formula()), collapse = " ")
)
required <- c(
    "eeg_mean | se(eeg_sem, sigma = TRUE)",
    "condition + (1 + condition | participant)",
    'bs = c("ds", "ps")', "k = c(16, 12)", "d = c(2, 1)",
    "m = list(c(1, 0.5), c(2, 2))", "full = FALSE"
)
stopifnot(all(vapply(required, grepl, logical(1),
                     x = formula_text, fixed = TRUE)))

off_002 <- threading_compare_paths(cfg, "off", 2L)
t4_002 <- threading_compare_paths(cfg, "threads4", 2L)
stopifnot(off_002$fit_file != t4_002$fit_file)
stopifnot(grepl("off/seed_002/fit.rds$", off_002$fit_file))
stopifnot(grepl("threads4/seed_002/fit.rds$", t4_002$fit_file))
```

- [ ] **Step 5: Run the test and verify GREEN**

Run:

```bash
Rscript --vanilla tests/test_threading_compare_ds_full.R
```

Expected: `threading comparison helper tests passed` with exit code 0.

- [ ] **Step 6: Commit Task 1**

```bash
git add tests/test_threading_compare_ds_full.R \
  experiments/threading_compare_ds_full/threading_compare_helpers.R
git commit -m "test: define ds threading comparison invariants"
```

---

### Task 2: Isolated generation and validation for seeds 2--10

**Files:**
- Modify: `tests/test_threading_compare_ds_full.R`
- Create: `experiments/threading_compare_ds_full/generate_inputs.R`
- Create: `experiments/threading_compare_ds_full/submit_generate.sh`

**Interfaces:**
- Consumes: `threading_compare_config()` and existing `make_paired_datasets()`
- Produces: `generate_threading_inputs(config) -> invisible(validation_table)`
- Produces: nine `sim_002.rds`--`sim_010.rds`, nine matching truth files, `manifest.rds`, `validation.csv`, and `complete.flag`

- [ ] **Step 1: Add failing generation-source tests**

Assert that the generation script exists, contains `make_paired_datasets`,
selects `datasets[["half-max"]]`, iterates over `config$seeds`, refuses a
non-empty output directory, and writes `complete.flag` last. Assert that the
Slurm file contains:

```text
#SBATCH --partition=ultralong
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
Rscript --vanilla generate_inputs.R
```

- [ ] **Step 2: Run the test and verify RED**

Run the test and expect failure because the generation files are absent.

- [ ] **Step 3: Implement generation with an all-or-stop manifest guard**

Source `code/00_generate_montecarlo.R`, create the output directory only when
absent or empty, and use:

```r
for (seed in config$seeds) {
    datasets <- make_paired_datasets(
        seed = seed, n_ppt = 20L, scenario = "high"
    )
    halfmax <- datasets[["half-max"]]
    saveRDS(halfmax$sim_data,
            file.path(config$data_dir, sprintf("sim_%03d.rds", seed)))
    saveRDS(halfmax$ground_truth,
            file.path(config$data_dir, sprintf("gt_%03d.rds", seed)))
    rm(datasets, halfmax)
    invisible(gc())
}
```

Validate each pair with exact checks for scenario `high`, radius rule
`half-max`, 20 participants, 61 channels, 51 time points after aggregation,
nine active electrodes, finite `eeg`, `xproj`, `yproj`, and `time`, and zero
continuous truth outside active electrodes. Write the manifest and validation
table before writing `complete.flag`.

- [ ] **Step 4: Implement the one-CPU generation wrapper**

The wrapper must set `PROJECT_ROOT`, use a generation-specific temporary
directory, load R 4.4.2, change to the experiment directory, and execute the
generation script. It must not submit downstream jobs.

- [ ] **Step 5: Verify GREEN and parse/syntax checks**

Run:

```bash
Rscript --vanilla tests/test_threading_compare_ds_full.R
Rscript --vanilla -e 'parse(file="experiments/threading_compare_ds_full/generate_inputs.R")'
bash -n experiments/threading_compare_ds_full/submit_generate.sh
```

Expected: all exit 0.

- [ ] **Step 6: Commit Task 2**

```bash
git add tests/test_threading_compare_ds_full.R \
  experiments/threading_compare_ds_full/generate_inputs.R \
  experiments/threading_compare_ds_full/submit_generate.sh
git commit -m "feat: generate isolated halfmax threading inputs"
```

---

### Task 3: Paired full-fit runner and Slurm arrays

**Files:**
- Modify: `tests/test_threading_compare_ds_full.R`
- Create: `experiments/threading_compare_ds_full/fit_seed.R`
- Create: `experiments/threading_compare_ds_full/submit_off.sh`
- Create: `experiments/threading_compare_ds_full/submit_threads4.sh`

**Interfaces:**
- Consumes: `THREADING_ARM=off|threads4`, `SLURM_ARRAY_TASK_ID=2..10`
- Produces: one guarded full fit and diagnostic bundle per arm-seed

- [ ] **Step 1: Add failing model-call and resource tests**

The source test must require `chains = 4`, `cores = 4`, `iter = 3000`,
`warmup = 1000`, `metric = "diag_e"`, `save_warmup = TRUE`, and conditionally
add `brms::threading(threads = 4L)` only for Arm `threads4`.

Require both shell scripts to contain `#SBATCH --array=2-10`,
`#SBATCH --constraint=xeon_e52640v4`, `#SBATCH --time=4-00:00:00`, and
`#SBATCH --mem=32G`. Require `--cpus-per-task=4` with no
`STAN_NUM_THREADS` in `submit_off.sh`, and `--cpus-per-task=16` plus
`STAN_NUM_THREADS=4` in `submit_threads4.sh`.

- [ ] **Step 2: Run the test and verify RED**

Run the test and expect failure because the runner and array files are absent.

- [ ] **Step 3: Implement one guarded arm-seed fit**

Build `common_arguments` exactly as follows:

```r
common_arguments <- list(
    formula = threading_compare_formula(),
    data = model_data,
    family = gaussian(),
    chains = config$chains,
    cores = config$cores,
    iter = config$iter,
    warmup = config$warmup,
    seed = seed,
    backend = "cmdstanr",
    control = list(
        adapt_delta = config$adapt_delta,
        max_treedepth = config$max_treedepth
    ),
    metric = config$metric,
    save_warmup = TRUE,
    output_dir = paths$sampler_dir,
    refresh = 100L
)
if (arm$threads == 4L) {
    common_arguments$threads <- brms::threading(threads = 4L)
}
fit <- do.call(brms::brm, common_arguments)
```

Before fitting, validate the isolated input checksums, metadata, contrasts,
dimensions, and finite values. Refuse any non-empty result directory. Use a
per-arm/per-seed cache and temporary directory. Save the fit before derived
diagnostics so a diagnostic failure cannot lose completed sampling.

- [ ] **Step 4: Record comparable diagnostics**

Write runtime fields `arm`, `seed`, `hostname`, `cpu_feature`, `allocated_cpus`,
`wall_seconds`, `wall_hours`, and `cpu_hours`. From `nuts_params()` write
divergences, treedepth hits, mean/median leapfrog, mean step size, and sampling
seconds per leapfrog where CmdStan elapsed metadata is available. From
`posterior::summarise_draws()` write maximum R-hat and minimum bulk/tail ESS,
then calculate ESS per wall hour and ESS per CPU-hour. Write `complete.flag`
last.

- [ ] **Step 5: Implement both array wrappers**

Both wrappers set `PROJECT_ROOT`, `THREADING_ARM`, seed-specific `TMPDIR`, and
arm-seed `CMDSTAN_CACHE_DIR`; load R 4.4.2; and run `fit_seed.R`. The `off`
wrapper must explicitly `unset STAN_NUM_THREADS`. The `threads4` wrapper must
export `STAN_NUM_THREADS=4`.

- [ ] **Step 6: Verify GREEN and static checks**

Run:

```bash
Rscript --vanilla tests/test_threading_compare_ds_full.R
Rscript --vanilla -e 'parse(file="experiments/threading_compare_ds_full/fit_seed.R")'
bash -n experiments/threading_compare_ds_full/submit_off.sh
bash -n experiments/threading_compare_ds_full/submit_threads4.sh
```

Expected: all exit 0.

- [ ] **Step 7: Commit Task 3**

```bash
git add tests/test_threading_compare_ds_full.R \
  experiments/threading_compare_ds_full/fit_seed.R \
  experiments/threading_compare_ds_full/submit_off.sh \
  experiments/threading_compare_ds_full/submit_threads4.sh
git commit -m "feat: add paired ds threading fit arrays"
```

---

### Task 4: Paired aggregation and throughput report

**Files:**
- Modify: `tests/test_threading_compare_ds_full.R`
- Create: `experiments/threading_compare_ds_full/aggregate_results.R`
- Create: `experiments/threading_compare_ds_full/submit_aggregate.sh`

**Interfaces:**
- Consumes: completed runtime, NUTS, and convergence CSVs from both arms
- Produces: `paired_seed_metrics.csv`, `aggregate_summary.csv`,
  `incomplete_pairs.csv`, and aggregate `complete.flag`

- [ ] **Step 1: Add a failing synthetic aggregation test**

Create two-arm synthetic rows for seeds 2 and 3 and assert paired ratios:

```r
fixture <- data.frame(
    arm = rep(c("off", "threads4"), each = 2),
    seed = rep(2:3, times = 2),
    wall_hours = c(28, 30, 14, 15),
    allocated_cpus = c(4, 4, 16, 16),
    min_ess_bulk = c(1000, 1100, 1000, 1100)
)
paired <- pair_threading_metrics(fixture)
stopifnot(
    all(paired$wall_ratio_threads4_over_off == 0.5),
    all(paired$cpu_hour_ratio_threads4_over_off == 2),
    all(paired$throughput_ratio_off_over_threads4 == 2)
)
```

- [ ] **Step 2: Run the test and verify RED**

Expected: failure because `pair_threading_metrics()` is absent.

- [ ] **Step 3: Implement complete-pair aggregation**

Join by seed, exclude any seed missing either completion flag, calculate:

```text
wall_ratio_threads4_over_off = threads4_wall / off_wall
cpu_hours                    = wall_hours * allocated_cpus
cpu_hour_ratio_threads4_over_off = threads4_cpu_hours / off_cpu_hours
throughput_ratio_off_over_threads4 = cpu_hour_ratio_threads4_over_off
ess_cpu_hour_ratio_off_over_threads4 = off_ess_per_cpu_hour / threads4_ess_per_cpu_hour
```

Report paired medians, means, ranges, and all nine seed-level values. Keep
incomplete seeds in a separate audit table rather than silently dropping them.

- [ ] **Step 4: Implement the one-CPU aggregate wrapper**

The wrapper requests one node, one task, one CPU, 16 GB, and one hour. It runs
only after both arrays have completed successfully and never modifies fits.

- [ ] **Step 5: Verify GREEN and static checks**

Run the full source test, parse all R files, and run `bash -n` on all four shell
files. Expected: all exit 0.

- [ ] **Step 6: Commit Task 4**

```bash
git add tests/test_threading_compare_ds_full.R \
  experiments/threading_compare_ds_full/aggregate_results.R \
  experiments/threading_compare_ds_full/submit_aggregate.sh
git commit -m "feat: aggregate paired threading throughput"
```

---

### Task 5: Remote preflight and mandatory user audit

**Files:**
- Read-only verification of all files created in Tasks 1--4
- Remote new directory: `/work/smjifan1/bgam_eeg_project/experiments/threading_compare_ds_full/`

**Interfaces:**
- Produces: verified remote scripts and a submission brief
- Does not produce: Slurm jobs

- [ ] **Step 1: Run fresh local verification**

```bash
Rscript --vanilla tests/test_threading_compare_ds_full.R
Rscript --vanilla -e 'for (f in Sys.glob("experiments/threading_compare_ds_full/*.R")) { parse(file=f); cat("parse_ok", f, "\n") }'
for f in experiments/threading_compare_ds_full/*.sh; do bash -n "$f"; done
git diff --check
```

Expected: every command exits 0.

- [ ] **Step 2: Upload only the new experiment files**

Create the new remote experiment directory and use `rsync -av` with an
explicit file list for the eight experiment files plus
`tests/test_threading_compare_ds_full.R`. Upload the test to the existing
remote `tests/` directory. Do not synchronize the dirty repository or existing
data, results, logs, or fits.

- [ ] **Step 3: Run remote preflight without generation or sampling**

On `gw01`, load R 4.4.2 and verify:

```bash
Rscript --vanilla tests/test_threading_compare_ds_full.R
Rscript --vanilla -e 'stopifnot(requireNamespace("brms"), requireNamespace("cmdstanr"), requireNamespace("posterior"))'
for f in experiments/threading_compare_ds_full/*.sh; do bash -n "$f"; done
test -x /work/smjifan1/.cmdstan/cmdstan-2.38.0/bin/stanc
sinfo -N -o "%f" | grep -q xeon_e52640v4
```

Also verify the isolated data and result directories do not yet exist or are
empty. Print all resolved `#SBATCH` directives.

- [ ] **Step 4: Calculate and present the maximum resource budget**

The audit must state:

```text
Generation: 1 job x 1 CPU
Arm off: 9 array tasks x 4 CPU = 36 concurrent CPU maximum
Arm threads4: 9 array tasks x 16 CPU = 144 concurrent CPU maximum
Both fit arrays together: 180 concurrent CPU maximum
Maximum requested fit CPU-hours at the 4-day limit:
  off      = 9 x 4 x 96 = 3,456 CPU-hours
  threads4 = 9 x 16 x 96 = 13,824 CPU-hours
```

Also provide realistic estimates separately from hard limits; do not present
estimates as guarantees.

- [ ] **Step 5: Present the exact commands and stop for approval**

Show, but do not run:

```bash
THREAD_COMPARE_GEN_JOB=$(sbatch --parsable \
  experiments/threading_compare_ds_full/submit_generate.sh)

THREAD_COMPARE_OFF_JOB=$(sbatch --parsable \
  --dependency=afterok:${THREAD_COMPARE_GEN_JOB} \
  experiments/threading_compare_ds_full/submit_off.sh)

THREAD_COMPARE_T4_JOB=$(sbatch --parsable \
  --dependency=afterok:${THREAD_COMPARE_GEN_JOB} \
  experiments/threading_compare_ds_full/submit_threads4.sh)

THREAD_COMPARE_AGG_JOB=$(sbatch --parsable \
  --dependency=afterok:${THREAD_COMPARE_OFF_JOB}:${THREAD_COMPARE_T4_JOB} \
  experiments/threading_compare_ds_full/submit_aggregate.sh)
```

Stop and request explicit approval. Do not execute any command in this step.

---

### Task 6: Submit only after explicit approval

**Files:**
- No file changes

**Interfaces:**
- Consumes: explicit user approval of the Task 5 audit
- Produces: four Slurm job IDs and initial queue verification

- [ ] **Step 1: Execute the exact approved dependency chain**

Run the four commands shown in Task 5 without altering resources or paths.

- [ ] **Step 2: Verify queue state and dependencies**

Use `squeue` and `scontrol show job` to confirm generation is eligible/running,
both arrays depend on generation, aggregation depends on both arrays, seeds are
2--10, and requested CPUs are 1, 4, 16, and 1 respectively.

- [ ] **Step 3: Report job IDs and monitoring paths**

Report all IDs, current states, dependency relationships, log patterns, output
paths, and the fact that seed 1 was untouched.
