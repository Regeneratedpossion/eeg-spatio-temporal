# Block-Orthogonalization Pilot Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Run a paired short pilot that tests whether right-orthogonal rotation of each penalized spline design block improves NUTS efficiency without changing the current BGAM.

**Architecture:** A focused R helper constructs the current half-max `ds x ps` brms formula and standata, rotates only matrices named `Zs_*`, and validates exact linear-predictor equivalence. One runner compiles the generated Stan program once and samples the baseline and rotated standata sequentially with identical settings; a summarizer reads CmdStan CSV diagnostics and reconstructs population-level linear predictors using each arm's own design matrices.

**Tech Stack:** R 4.4, brms 2.23, mgcv, cmdstanr, posterior, base R SVD, Slurm.

## Global Constraints

- Use half-max hard-cutoff high-heterogeneity seed 1 data with 20 participants.
- Do not modify, stop, or write into the formal Seed 22--100 experiment.
- Use the current `ds x ps`, `k = c(16, 12)`, `m = list(c(1, 0.5), c(2, 2))`, `full = FALSE` model.
- Apply right-orthogonal rotations only; never rescale singular values.
- Both arms use one chain, one core, no within-chain threading, 300 warmup, 100 sampling, `adapt_delta = 0.80`, `max_treedepth = 10`, `diag_e`, and the same sampling seed.
- Compile the Stan program once and exclude compilation from both arm runtimes.
- Store all new outputs under `experiments/block_orthogonalization_pilot/`.

---

### Task 1: Standata Rotation and Equivalence Checks

**Files:**
- Create: `experiments/block_orthogonalization_pilot/pilot_helpers.R`
- Test: `tests/test_block_orthogonalization_pilot.R`

**Interfaces:**
- Produces: `pilot_formula()`, `rotate_smooth_blocks(standata, tolerance)`, `validate_rotated_standata(original, rotated, rotations, tolerance)`, and `pilot_paths(project_root)`.
- `rotate_smooth_blocks()` returns `list(standata, rotations, diagnostics)`; each rotation entry contains `name`, `V`, singular values, and pre/post Gram off-diagonal diagnostics.

- [ ] **Step 1: Write failing structural and numerical tests**

Create a test that sources the helper and checks the formula text contains `bs = c("ds", "ps")`, `k = c(16, 12)`, and the specified `m`. Add a synthetic full-rank matrix test:

```r
set.seed(1)
Z <- matrix(rnorm(240), 40, 6)
sd <- list(Zs_1_1 = Z, untouched = 4L)
out <- rotate_smooth_blocks(sd)
beta <- rnorm(6)
rot_beta <- crossprod(out$rotations$Zs_1_1$V, beta)
stopifnot(max(abs(Z %*% beta - out$standata$Zs_1_1 %*% rot_beta)) < 1e-10)
stopifnot(identical(out$standata$untouched, sd$untouched))
stopifnot(max(abs(crossprod(out$rotations$Zs_1_1$V) - diag(6))) < 1e-10)
```

- [ ] **Step 2: Run the test and verify it fails because the helper is absent**

Run:

```bash
Rscript --vanilla tests/test_block_orthogonalization_pilot.R
```

Expected: non-zero exit with failure to source `pilot_helpers.R`.

- [ ] **Step 3: Implement the minimal rotation helper**

For every matrix whose name matches `^Zs_[0-9]+_[0-9]+$`, compute:

```r
decomp <- svd(Z, nu = 0L, nv = ncol(Z))
V <- decomp$v
Z_rotated <- Z %*% V
```

Keep all other standata byte-for-byte unchanged. Validate dimensions, finite values, orthogonality of `V`, equality of the transformed linear predictor, and near-diagonality of `crossprod(Z_rotated)`. Use a relative tolerance scaled by matrix norm.

- [ ] **Step 4: Run the helper test and verify it passes**

Run:

```bash
Rscript --vanilla tests/test_block_orthogonalization_pilot.R
```

Expected: exit 0 and an explicit pass message.

### Task 2: Shared-Compilation Paired Runner

**Files:**
- Create: `experiments/block_orthogonalization_pilot/run_pilot.R`
- Modify: `tests/test_block_orthogonalization_pilot.R`

**Interfaces:**
- Consumes: all Task 1 helper functions.
- Produces: arm-specific CmdStan CSV files, `standata_baseline.rds`, `standata_rotated.rds`, `rotations.rds`, `design_diagnostics.csv`, `runtime.csv`, and `run_config.rds` under the isolated experiment directory.

- [ ] **Step 1: Extend the test with a dry-run contract**

Source `run_pilot.R` with `PILOT_DRY_RUN=true` and assert that it constructs two arms named `baseline` and `rotated`, both with:

```r
list(chains = 1L, parallel_chains = 1L, iter_warmup = 300L,
     iter_sampling = 100L, adapt_delta = 0.80,
     max_treedepth = 10L, metric = "diag_e", threads_per_chain = 1L)
```

Assert that the data path resolves to `data/montecarlo_rhalfmax/sim_001.rds` and that no formal Seed 22--100 path appears in any output path.

- [ ] **Step 2: Run the test and verify the new dry-run assertions fail**

Run the same test command and expect a non-zero exit because `run_pilot.R` is absent.

- [ ] **Step 3: Implement the runner**

The runner must:

1. aggregate trials and apply `contrasts(condition) <- c(-0.5, 0.5)`;
2. call `make_standata()` and `make_stancode()` once;
3. create and validate the rotated copy;
4. compile one `cmdstan_model` before starting either timer;
5. sample baseline, then rotated, with separate `output_dir` values and the same seed;
6. time only `model$sample()` for each arm;
7. set `save_warmup = TRUE` so adaptation and warmup progress remain auditable;
8. stop immediately on non-zero CmdStan return codes or missing CSV files.

- [ ] **Step 4: Run all dry-run tests and parse the runner**

Run:

```bash
Rscript --vanilla tests/test_block_orthogonalization_pilot.R
Rscript -e 'parse(file="experiments/block_orthogonalization_pilot/run_pilot.R")'
```

Expected: both exit 0.

### Task 3: Sampler and Prediction Comparison

**Files:**
- Create: `experiments/block_orthogonalization_pilot/summarize_pilot.R`
- Modify: `tests/test_block_orthogonalization_pilot.R`

**Interfaces:**
- Consumes: saved arm CSV files and standata from Task 2.
- Produces: `sampler_comparison.csv`, `prediction_comparison.csv`, and a concise stdout conclusion.

- [ ] **Step 1: Add failing unit tests for sampler summaries and linear-predictor reconstruction**

Use small synthetic draws and design matrices to verify that the reconstruction function sums `Intercept`, `Xc * b`, `Xs * bs`, and every `Zs_* * s_*` contribution. Verify baseline coefficients and rotated coefficients related by `V' beta` yield identical predictions.

- [ ] **Step 2: Run tests and verify failure because summarizer functions are absent**

Run the same R test file and expect non-zero exit.

- [ ] **Step 3: Implement summary functions**

For each arm, calculate:

- elapsed seconds and seconds per 100 total iterations;
- adapted step size;
- mean, median, and maximum `n_leapfrog__`;
- mean and maximum `treedepth__`, plus count at 10;
- number of `divergent__ == 1` transitions.

Reconstruct population-level linear predictors for the 3,111 unique electrode-time-condition rows from saved transformed parameter draws. For the rotated arm, use its rotated `Zs_*` matrices directly. Compare posterior means using maximum absolute difference, RMSE, and CZ at 250 ms. Report Monte Carlo uncertainty and do not require draw-wise equality between independent chains.

- [ ] **Step 4: Run tests and parse all experiment scripts**

Run:

```bash
Rscript --vanilla tests/test_block_orthogonalization_pilot.R
Rscript -e 'for (f in Sys.glob("experiments/block_orthogonalization_pilot/*.R")) parse(file=f)'
```

Expected: exit 0.

### Task 4: TU Dortmund Slurm Submission

**Files:**
- Create: `experiments/block_orthogonalization_pilot/submit_pilot.sh`

**Interfaces:**
- Consumes: experiment scripts from Tasks 1--3 and the existing remote R/CmdStan environment.
- Produces: one Slurm job that runs both arms sequentially and then summarizes them.

- [ ] **Step 1: Implement the Slurm script**

Request one node, one CPU, sufficient memory for the existing standata, and a time limit long enough for two 400-iteration arms. Load the same R module and `R_LIBS_USER` used by the validated TU Dortmund formal jobs. Run `run_pilot.R`, then `summarize_pilot.R`, with logs isolated under the experiment directory.

- [ ] **Step 2: Validate shell and R syntax locally**

Run:

```bash
bash -n experiments/block_orthogonalization_pilot/submit_pilot.sh
Rscript --vanilla tests/test_block_orthogonalization_pilot.R
```

Expected: exit 0.

- [ ] **Step 3: Upload only the pilot files and submit one job**

Use `rsync` to upload the new experiment directory and required test-independent helper dependencies. Confirm the remote data file and R/CmdStan environment before `sbatch`. Do not upload over formal experiment scripts.

- [ ] **Step 4: Verify accepted queue state and paths**

Run `squeue` and `scontrol show job`; verify one CPU, the intended account/partition, an isolated work directory, and no `BadConstraints` or immediate failure. Record the job ID in the handoff response.
