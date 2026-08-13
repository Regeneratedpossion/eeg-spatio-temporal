# Block-Orthogonalization Pilot Design

## Goal

Test whether a model-preserving rotation of the penalized spline design
blocks improves NUTS efficiency for the current half-max BGAM.

## Scope

- Dataset: half-max hard cutoff, high-heterogeneity scenario, `n_ppt = 20`,
  seed 1.
- The formal Seed 22--100 simulation is not modified or interrupted.
- The pilot changes no estimand, likelihood, prior, spline basis, penalty,
  basis dimension, decision rule, or diagnostic rule.

## Paired Arms

Both arms use the current `ds x ps` BGAM:

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

- Arm A: unmodified `brms::make_standata()` design matrices.
- Arm B: for every penalized block `Zs_*`, compute the thin SVD
  `Z = U D V'` and replace it with `Z V = U D`.

Only the right-orthogonal rotation is applied. Singular values are not
rescaled. Because each corresponding latent coefficient vector has an
isotropic standard-normal prior, this rotation preserves the prior and the
implied smooth-function distribution.

## Sampling Configuration

Both arms use:

- one chain;
- 300 warmup iterations and 100 sampling iterations;
- one CPU core and no within-chain threading;
- `adapt_delta = 0.80`;
- `max_treedepth = 10`;
- the same seed and initialization policy;
- diagonal Euclidean metric.

The generated Stan program is shared between arms; only the standata design
matrices differ.

## Validation

Before sampling:

1. Verify each rotation matrix is orthogonal.
2. Verify `Z %*% beta` equals `(Z %*% V) %*% (t(V) %*% beta)` numerically.
3. Verify the two arms have identical dimensions and all non-`Zs_*` standata.
4. Verify the rotated blocks have orthogonal columns to numerical tolerance.

After sampling, compare:

- wall-clock time and seconds per 100 iterations;
- adapted step size;
- mean, median, and maximum leapfrog count;
- treedepth distribution and maximum-treedepth hits;
- divergences;
- posterior predictions after applying the required inverse coefficient
  rotation.

Prediction agreement is a model-equivalence check, not an expectation that
two short independent MCMC runs produce identical posterior summaries.

## Decision Rule

The rotation is considered computationally promising only if it reduces
leapfrog work and wall time without divergences or a clear discrepancy in
posterior predictions. A short one-chain pilot is not evidence to replace the
formal model; any promising result requires a longer paired validation before
future simulations could use it.

## Outputs

Use an isolated experiment directory containing the arm configuration,
rotation metadata, sampler CSV files, runtime summary, and comparison table.
Do not write into the formal Seed 22--100 result directories.
