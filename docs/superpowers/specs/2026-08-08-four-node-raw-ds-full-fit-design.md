# Four-Node Raw DS Full-Fit Design

## Objective

Re-run the existing half-max, high-heterogeneity, `n_ppt = 20`, seed-1
unstandardized `ds x ps` BGAM without changing the statistical model. Change
only execution topology: replace four chains sharing one node with four
single-chain R processes placed on four distinct Slurm nodes.

This is a full fit, not a short benchmark. Each chain runs 3,000 iterations,
including 1,000 warmup and 2,000 retained sampling iterations. No job is
submitted until the user has reviewed the implementation diff and explicitly
approved the exact `sbatch` command.

## Fixed scientific inputs and model

The source objects remain:

- `data/montecarlo_rhalfmax/sim_001.rds`
- `data/montecarlo_rhalfmax/gt_001.rds`

The R script aggregates trials with the existing `aggregate_trials()` helper,
uses raw `xproj`, `yproj`, and `time`, and assigns the existing sum contrast:

```text
cond1 = -0.5
cond2 =  0.5
```

The formula remains exactly:

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

The family, default brms priors, `adapt_delta = 0.8`,
`max_treedepth = 10`, diagonal Euclidean metric, saved warmup, and four
within-chain threads remain unchanged.

## Execution architecture

One Slurm allocation requests four nodes and four tasks with exactly one task
per node:

```text
nodes             = 4
ntasks             = 4
ntasks-per-node    = 1
cpus-per-task      = 4
threads-per-chain  = 4
```

An `srun` step launches one R process per task. Each process reads
`SLURM_PROCID`, maps task IDs `0:3` to CmdStan chain IDs `1:4`, and fits one
chain with `chains = 1` and `cores = 1`.

Every process uses base seed 1 together with its unique `chain_ids` value.
This preserves CmdStan's distinct chain-specific random streams and prevents
four separately launched R processes from generating identical chains.

Each process uses a chain-specific temporary directory, CmdStan cache, sampler
directory, log, and result directory. Concurrent compilation and output writes
therefore cannot race on the same files.

## Outputs and combination

All artifacts live under the isolated experiment directory:

```text
experiments/unstandardized_ds_full_halfmax_seed001_four_nodes/
```

Each `results/chain_<id>/` directory contains:

- the single-chain `brmsfit`;
- raw CmdStan CSV, metric, and configuration artifacts;
- resolved model configuration;
- runtime including hostname and chain ID;
- NUTS summary;
- a chain completion flag.

After the four-task `srun` step exits successfully, a separate one-task step
validates all four completion flags and checks that the recorded hostnames are
four distinct nodes. It then calls `brms::combine_models()` with data checking
enabled, saves the combined four-chain `brmsfit`, and writes combined runtime
and NUTS summaries.

The experiment-level completion flag is written only after the combined fit
and summaries have been saved successfully.

## Failure and overwrite protection

Preparation or fitting fails if any of these conditions holds:

- source data or truth is absent;
- scenario, radius rule, participant count, channel count, or time count is
  wrong;
- raw predictors contain missing or non-finite values;
- the condition contrast differs from `-0.5/+0.5`;
- a chain ID is outside `1:4`;
- a chain result directory is non-empty;
- any chain process exits non-zero;
- any chain completion flag is missing;
- hostnames are not four distinct nodes;
- model data differ when `combine_models()` checks the fits;
- an experiment-level completed output already exists.

`srun --kill-on-bad-exit=1` stops the distributed step when a chain fails.
Partial chain artifacts remain available for diagnosis, but neither combined
output nor the experiment completion flag is produced.

The currently running one-node standardized and unstandardized jobs are
read-only and are not cancelled, modified, or overwritten.

## Verification and review gates

Before submission, the implementation must pass:

1. a source-level R test proving the formula, raw predictors, sampler controls,
   chain-ID mapping, output isolation, and completion guards;
2. an R parse check for the fitting and combination scripts;
3. `bash -n` for the Slurm script;
4. local and remote checks that the two input RDS files exist;
5. a dry inspection of the resolved Slurm directives and exact `srun` command;
6. a diff review showing every new file and confirming no existing fit script
   was modified.

After these checks, work stops. The user receives the implementation diff, the
verified file list, the exact proposed `sbatch` command, and a concise account
of all changes. Submission occurs only after a new explicit approval.

## Success criteria

The run is successful when four chains finish on four distinct nodes, combine
into one valid four-chain `brmsfit`, and produce complete runtime and NUTS
summaries with no overwrite of existing experiments. The primary performance
comparison is per-chain wall time and sampling seconds per iteration against
the currently running one-node raw fit; convergence diagnostics remain
scientific acceptance criteria rather than being relaxed for speed.
