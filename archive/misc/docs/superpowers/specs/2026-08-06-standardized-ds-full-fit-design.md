# Standardized DS x PS Full-Fit Design

## Objective

Prepare one full BGAM fit for the existing half-max, high-heterogeneity,
`n_ppt = 20`, seed-1 pilot replicate. Standardize the continuous predictors,
retain the selected `ds x ps` smooth specification, and stop for a documented
pre-fit audit before submitting any sampling job.

## Scope and invariants

The source objects remain:

- `data/montecarlo_rhalfmax/sim_001.rds`
- `data/montecarlo_rhalfmax/gt_001.rds`

The source RDS files will not be modified or duplicated. Trial aggregation,
known standard errors, likelihood, priors, random-effects structure, seed,
basis dimensions, sampler control, and evaluation rules remain unchanged.

The model uses:

```r
eeg_mean | se(eeg_sem, sigma = TRUE) ~
  condition +
  (1 + condition | participant) +
  t2(
    x_z, y_z, time_z,
    by = condition,
    bs = c("ds", "ps"),
    k = c(16, 12),
    d = c(2, 1),
    m = list(c(1, 0.5), c(2, 2)),
    full = FALSE
  )
```

The intended full sampling configuration is four chains, 3,000 iterations per
chain, 1,000 warmup iterations, seed 1, `adapt_delta = 0.8`,
`max_treedepth = 10`, and diagonal Euclidean metric. The within-chain thread
count will be recorded explicitly in the audit and will not be inferred from
machine defaults.

## Preprocessing

Condition is explicitly assigned sum contrasts:

```text
cond1 = -0.5
cond2 =  0.5
```

Scaling constants are calculated once from the unique design points, not from
the repeated aggregated rows:

- separate means for `xproj` and `yproj`;
- one shared spatial scale,
  `sqrt((var(xproj) + var(yproj)) / 2)`, to preserve the electrode geometry;
- the mean and standard deviation of the unique time points.

The transformed predictors are:

```r
x_z    = (xproj - x_center) / space_scale
y_z    = (yproj - y_center) / space_scale
time_z = (time - time_center) / time_scale
```

The constants are saved as an RDS object and must be reused unchanged for all
posterior prediction data.

## Pre-fit audit gate

Preparation produces a human-readable audit without sampling. The audit must
verify:

1. input paths, scenario, radius rule, seed, participant count, channel count,
   and time-point count;
2. source-file checksums so the fitted replicate is identifiable;
3. original and transformed means, standard deviations, ranges, missingness,
   and finite-value checks;
4. the exact condition contrast matrix;
5. preservation of pairwise spatial distances up to the single common scale;
6. saved scaling constants and their output path;
7. formula, priors, family, chain/iteration/control settings, requested CPUs,
   and within-chain threads;
8. a new isolated output directory with an overwrite guard;
9. successful R parsing and a zero-chain Stan compilation check.

The workflow stops after producing this audit. Sampling is submitted only
after the user reviews and approves it.

## Outputs

All new artifacts use an isolated experiment directory and do not overwrite
the completed A--F pilot arms or the original 61.6-hour fit. The directory
contains the preprocessing constants, audit tables, resolved configuration,
compiled-model cache, sampler CSV files, final fit, runtime summary, NUTS
summary, and completion flag.

## Failure handling

Preparation fails before submission if an invariant is false, any transformed
value is missing or non-finite, contrasts differ from `-0.5/+0.5`, spatial
distance preservation fails tolerance, required input files are absent, or a
completed output already exists. A failed or incomplete run must not write a
completion flag.
