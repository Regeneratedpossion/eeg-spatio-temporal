# DS Threading Throughput Comparison, Seeds 2--10

## Objective

Run a paired, full-sampling comparison of two CPU-allocation strategies for
the half-max, high-heterogeneity, `n_ppt = 20`, raw-coordinate `ds x ps` BGAM.
Run simulation seeds 2 through 10 under both strategies. Seed 1 is not rerun
and is excluded from the paired comparison.

The comparison answers whether cluster CPUs produce more completed simulation
fits when used for independent models and chains instead of within-chain
threading. No generation or fitting job is submitted until the user reviews
the implementation, resource requests, output paths, and exact `sbatch`
commands.

## Isolated simulation inputs

The current half-max directory contains only seed 1 and has an existing
manifest. It remains read-only. Seeds 2--10 are generated once into:

```text
data/threading_compare_halfmax_high_p20/
```

For each seed, the existing paired hard-radius generator creates the complete
paired replicate in memory; only its half-max `sim_<seed>.rds` and
`gt_<seed>.rds` objects are saved in the isolated directory. The manifest fixes
the scenario, participant count, trial count, signal mode, radius rule,
generation parameters, and replicate seeds 2--10.

Generation refuses to overwrite any existing simulation, truth, or manifest.
Validation checks all nine seed pairs, their metadata, dimensions, finite
values, expected nine active electrodes, and exact zero signal outside the
half-max support before any fit is eligible to start.

## Fixed statistical model

Both arms use the same data and model:

```r
eeg_mean | se(eeg_sem, sigma = TRUE) ~
  condition +
  (1 + condition | participant) +
  t2(
    xproj, yproj, time,
    by = condition,
    bs = c("ds", "ps"),
    k = c(16, 12),
    d = c(2, 1),
    m = list(c(1, 0.5), c(2, 2)),
    full = FALSE
  )
```

Both arms retain raw predictors, sum contrasts `cond1 = -0.5` and
`cond2 = 0.5`, Gaussian likelihood with known SEM and estimated residual
sigma, brms default priors, four chains, 3,000 iterations per chain, 1,000
warmup iterations, `adapt_delta = 0.8`, `max_treedepth = 10`, `diag_e`, and
saved warmup draws. The MCMC base seed equals the simulation seed in both arms.

## Comparison arms

### Arm A: no within-chain threading

Each Slurm array task fits one simulation seed on one node:

```text
chains             = 4
parallel chains    = 4
within-chain       = disabled
CPUs per task      = 4
```

The `threads` argument is omitted from `brm()`, and `STAN_NUM_THREADS` is not
set. The Stan model is compiled without `reduce_sum` threading. Array indices
are seeds 2--10.

### Arm B: four within-chain threads

Each Slurm array task fits one simulation seed on one node:

```text
chains             = 4
parallel chains    = 4
threads per chain  = 4
CPUs per task      = 16
```

The fit uses `threads = brms::threading(threads = 4)` and
`STAN_NUM_THREADS = 4`. Array indices are seeds 2--10.

Both arrays use the `ultralong` partition and the same CPU-family constraint so
that hardware type is held fixed. They use separate arm-and-seed temporary
directories, compilation caches, logs, sampler directories, and fit outputs.

## Outputs and overwrite protection

All fitting outputs live under:

```text
results/threading_compare_ds_full/
  off/seed_002/ ... seed_010/
  threads4/seed_002/ ... seed_010/
```

Each seed-arm directory contains the `brmsfit`, raw sampler artifacts,
resolved configuration, source-input checksums, runtime summary, NUTS summary,
posterior convergence summary, and a completion flag. A non-empty seed-arm
directory or existing completion flag aborts that task rather than overwriting
it.

Generation, Arm A, and Arm B are separate jobs. The two fit arrays are eligible
only after the generation-and-validation job exits successfully. A failed seed
does not overwrite or invalidate completed paired seeds, but the aggregate
comparison excludes a seed unless both arms completed successfully.

## Comparison metrics

For every paired seed and arm, record:

- wall-clock hours;
- allocated CPUs and CPU-hours;
- warmup and sampling seconds;
- seconds per sampling iteration and per leapfrog;
- divergences, maximum-treedepth hits, step size, and leapfrog counts;
- maximum R-hat and minimum bulk/tail ESS;
- ESS per wall hour and ESS per CPU-hour;
- hostname and CPU feature.

The aggregate report uses seeds with both arms complete and reports paired
ratios for wall time, CPU-hours, ESS efficiency, and model throughput. It also
compares selected posterior summaries to detect material numerical differences
between the threaded and non-threaded executions.

## Verification and submission gate

Before submission, the implementation must pass source-level tests for model
identity, arm isolation, seed mapping, resource mapping, manifest guards, and
aggregation rules. All R files must parse, all shell files must pass `bash -n`,
and a remote preflight must verify packages, CmdStan, input/output paths,
cluster constraints, and resolved Slurm directives.

After preflight, work stops for user review. The review includes every new
file, the exact model calls for both arms, requested nodes/CPUs/memory/time,
estimated maximum CPU-hours, output locations, dependency structure, and the
exact generation and array `sbatch` commands. Submission requires a new,
explicit user approval.
