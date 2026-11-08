# Laplacian K=20 + Constant Basis Result

Date recorded: 2026-06-20

## Model

Dataset:

- `montecarlo_p20_gaussian_hardlabel`
- `ppt = 20`
- continuous Gaussian spatial signal with hard active labels available for secondary threshold diagnostics

Spatial model:

- Graph-Laplacian low-rank spatial basis
- `K = 20`
- Includes the constant/global spatial eigenvector
- This means the basis contains:
  - 1 constant/global component
  - 19 non-trivial graph-Laplacian spatial eigenvectors

Bayesian temporal model per basis score:

```r
score | se(score_sem) ~
  condition +
  s(time, by = condition, k = 10) +
  (1 + condition | participant)
```

Sampling:

- `iter = 3000`
- `warmup = 1000`
- `adapt_delta = 0.999`
- `max_treedepth = 13`

Remote result directory:

```text
/work/smjifan1/bgam_eeg_project/results/laplacian_basis_brms_p20_gaussian_hardlabel_constant_long/k20
```

## Main Diagnostics

Hard-label diagnostics:

```text
convergence_pass        = FALSE
temporal_pass           = TRUE
spatial_top_pass        = TRUE
spatial_F1              = 0.3396
spatiotemporal_F1       = 0.5451
peak_recovery_pct       = 81.4%
mean_recovery_pct       = 95.1%
```

Continuous recovery diagnostics:

```text
spatiotemporal_cor_all              = 0.962
spatiotemporal_cor_effect_window    = 0.964
spatiotemporal_rmse                 = 0.0446
spatiotemporal_mae                  = 0.0278
spatiotemporal_cosine               = 0.965
spatiotemporal_soft_dice            = 0.722

peak_spatial_cor                    = 0.976
peak_spatial_cosine                 = 0.985
peak_spatial_soft_dice              = 0.850

truth_peak_uV                       = 1.500
estimate_peak_uV                    = 1.221
peak_recovery_pct                   = 81.4%

truth_weighted_mean_recovery_pct    = 104.5%
truth_top25_mean_recovery_pct       = 102.7%
```

## Convergence Details

The model did not fully pass convergence because of a small number of divergent transitions.

Original constant K=20 long run:

```text
total divergent = 16
lap11 divergent = 2
lap19 divergent = 14
max Rhat <= 1.0097
min ESS bulk ~= 956
treedepth hits = 0
```

Repair attempt:

```text
lap11: 2 divergent -> 0 divergent
lap19: 14 divergent -> 29 divergent
```

Simplified `lap19` tests:

```text
random intercept only: n_divergent = 15
no random effect:      n_divergent = 9
```

Interpretation:

- `lap11` can be fixed by targeted rerun.
- `lap19` remains unstable even under simpler random-effect structures.
- `lap19` is a weak/noisy high-order basis rather than a major signal component.

Score diagnostics for `lap19`:

```text
score_sd                  = 0.226
sem_mean                  = 0.223
z_abs_p99                 = 2.64
max_abs_condition_diff    = 0.137
participant_diff_sd       = 0.030
```

This suggests the divergent behavior is not caused by an obvious data outlier, but by weak-signal posterior geometry for a high-order spatial basis.

## Why This Result Matters

Compared with the previous K=20 Laplacian basis without the constant component, amplitude recovery improved substantially.

Previous K=20 without constant:

```text
peak recovery                 ~= 69.4%
truth-weighted mean recovery  ~= 67.8%
top-25% mean recovery         ~= 75.9%
```

K=20 with constant:

```text
peak recovery                 ~= 81.4%
truth-weighted mean recovery  ~= 104.5%
top-25% mean recovery         ~= 102.7%
```

Oracle projection showed that the earlier amplitude shrinkage was largely caused by omitting the constant/global spatial component:

```text
K=20 without constant:
  theoretical peak recovery ~= 68.9%

K=20 with constant:
  theoretical peak recovery ~= 81.1%
```

Therefore, the amplitude shrinkage was mostly a basis-representation issue, not primarily a Bayesian sampling issue.

## Current Interpretation for Thesis

The `K=20 + constant` model is the current best main result for the low-rank Laplacian approach:

- It is computationally feasible.
- It recovers temporal localization.
- It identifies the correct top spatial channel/region.
- It substantially improves continuous amplitude recovery.
- It preserves the sensor-level posterior reconstruction workflow.

Main limitation:

- It still has a small unresolved convergence issue in one weak high-order basis (`lap19`).

Recommended phrasing:

> Including the constant graph-Laplacian component substantially improved amplitude recovery, indicating that the previous shrinkage was primarily due to the restricted spatial basis rather than the temporal Bayesian GAMM itself. The remaining sampler issue was localized to a weak high-order basis component and did not materially affect the reconstructed continuous effect map.

## Next Sensitivity Checks

Recommended next checks:

1. `K=19 + constant`, excluding unstable `lap19`.
2. `K=30 + constant` quick screen, to test whether additional spatial detail improves peak recovery.
3. Compare continuous recovery metrics rather than relying on hard active/inactive labels.

