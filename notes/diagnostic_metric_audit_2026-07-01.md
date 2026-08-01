# BGAM Diagnostic Metric Audit, 2026-07-01

## Question 1: Do F1 and onset/offset diagnostics use model output?

The main diagnostic script is `code/archive/11_run_bgam_diagnostics.R`.

- `posterior_condition_diff()` calls `add_epred_draws(fit, re_formula = NA)` and computes `cond2 - cond1` from posterior expected values.
- `onset_bias_ms` and `offset_bias_ms` are computed from clusters in posterior odds at the seed electrode.
- `spatial_F1` is computed at peak time from channels whose posterior lower 95% interval is above zero.
- `spatiotemporal_F1` is computed from channel-time cells whose posterior odds exceed `ODDS_THRESHOLD`.

Therefore these diagnostics do read the fitted model. They are not directly copied from grid coordinates or a fixed P-only constant.

However, `spatial_F1` is a coarse thresholded count. Two different fits can produce the same count if they cross the same active-channel threshold set.

## P20 Manual Recheck

Two P20/B_ps fits were re-evaluated into `results/diagnostics/audit_p20_f1_recheck/`.

| Recheck job | Original run | Configuration | Result |
|---|---:|---|---|
| 39983284 | 39961266_4 | P20 B_ps default-control | completed in 00:01:38 |
| 39983285 | 39962984_4 | P20 B_ps sigmaTRUE default-control | completed in 00:01:40 |

Summary:

| Configuration | onset_bias_ms | offset_bias_ms | spatial_F1 | spatiotemporal_F1 | peak_recovery_pct | mean_recovery_pct |
|---|---:|---:|---:|---:|---:|---:|
| P20 B_ps default-control | +10 | -12 | 0.3396226415 | 0.5275862069 | 84.0484 | 99.5168 |
| P20 B_ps sigmaTRUE | +10 | -12 | 0.3396226415 | 0.5284974093 | 83.8867 | 99.4014 |

The full posterior-derived files are not identical:

- spatial max absolute posterior mean difference: 0.008463968
- spatial max absolute lower-CI difference: 0.007123108
- spatiotemporal max absolute posterior mean difference: 0.01071224
- spatiotemporal max absolute odds difference: 7200

Interpretation: the diagnostic code uses model output, but `spatial_F1` is insensitive here because both P20/B fits produce the same peak-time confusion counts: TP=9, FP=35, FN=0.

## Question 2: Is P5 structurally overparameterized because it has too few electrodes?

No. In these datasets `P=5` and `P=20` are participant counts, not electrode counts. Both datasets have 61 channels and 51 time points.

The smooth basis design over unique channel-time grid points is identical for P5 and P20:

| Dataset | Basis | Rows | Columns | Rank | Condition number | Max abs column correlation |
|---|---|---:|---:|---:|---:|---:|
| P5 hard-label | tp k=12 | 3111 | 143 | 143 | 4285.3 | 0.9912 |
| P5 hard-label | ps k=12 | 3111 | 143 | 143 | 5698.9 | 0.9881 |
| P20 hard-label | tp k=12 | 3111 | 143 | 143 | 4285.3 | 0.9912 |
| P20 hard-label | ps k=12 | 3111 | 143 | 143 | 5698.9 | 0.9881 |

Interpretation: k=(12,12) is highly collinear, but not uniquely because of P5 electrode count. The P5 weakness is more plausibly fewer participants / less replication supporting the same high-dimensional smooth and participant random-effect structure.

## Question 3: Tensor-product A treedepth explosion

The A/tp k=12 smooth basis has high condition number and near-duplicate directions. A quick k sweep on P5 shows that reducing k sharply improves conditioning:

| Basis | k | Columns | Rank | Condition number | Max abs column correlation |
|---|---:|---:|---:|---:|---:|
| tp | 6 | 35 | 35 | 287.2 | 0.9812 |
| tp | 8 | 63 | 63 | 616.6 | 0.9828 |
| tp | 10 | 99 | 99 | 1150.3 | 0.9867 |
| tp | 12 | 143 | 143 | 4285.3 | 0.9912 |
| ps | 6 | 35 | 35 | 424.9 | 0.9805 |
| ps | 8 | 63 | 63 | 726.4 | 0.9766 |
| ps | 10 | 99 | 99 | 1456.7 | 0.9764 |
| ps | 12 | 143 | 143 | 5698.9 | 0.9881 |

Interpretation: the A/tp treedepth explosion is consistent with weakly identified smooth directions at k=12. For the next round, test k=8 or k=10 before more sampler tuning.

## Question 4: Pending recovery diagnostics

The three previously pending rows have been filled in `results/tables/ssh_model_run_inventory.csv`.

| Run | Recovery diagnostic status |
|---|---|
| 39961266_2 | filled from existing remote diagnostics |
| 39962984_1 | filled after diagnostic job 39983282 |
| 39965039_2 | filled after diagnostic job 39983283 |

This enables the P20/B_ps comparison:

| Run | Sigma mode | adapt_delta | spatial_F1 | spatiotemporal_F1 | mean recovery |
|---|---|---:|---:|---:|---:|
| 39961266_4 | default-control | 0.80 | 0.340 | 0.528 | 99.52% |
| 39962984_4 | sigmaTRUE | 0.80 | 0.340 | 0.528 | 99.40% |
| 39965039_2 | sigmaTRUE | 0.95 | 0.340 | 0.528 | 99.41% |

Interpretation: for P20/B_ps, sigma handling and adapt_delta mainly changed sampler diagnostics/runtime, not recovery quality.

## Question 5: Wall-clock structure

`results/tables/ssh_model_run_inventory.csv` now includes a numeric `model_elapsed_sec` column next to the human-readable `wall_clock` field.
