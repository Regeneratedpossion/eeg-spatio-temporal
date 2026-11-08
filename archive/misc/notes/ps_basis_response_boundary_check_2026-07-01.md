# P-spline Basis Response Boundary Check, 2026-07-01

Purpose: check whether high-correlation P-spline tensor-product basis columns are spatially concentrated near scalp boundaries.

## Result

For `t2(xproj, yproj, time, bs = c("ps", "ps"))`, high-correlation basis pairs are not randomly distributed.

Using `abs(correlation) > 0.97`:

| basis | k | high-corr pairs | dominant max-response channels |
|---|---:|---:|---|
| ps | 8 | 6 | FPZ / AF1 |
| ps | 10 | 8 | FPZ / AF1 |
| ps | 12 | 11 | FPZ / AF1 |

The pattern is consistent across k. The paired columns are almost perfect sign reversals on the relevant time traces:

| basis | k | time-profile correlation at FPZ | time-profile correlation at AF1 |
|---|---:|---:|---:|
| ps | 8 | -1.000 | -1.000 |
| ps | 12 | -1.000 | -1.000 |

## Interpretation

The residual high correlation in the `ps` basis is a structured boundary-response issue, not a random pair of close physical electrodes.

Important correction: `FPZ/AF1` are not physically nearest neighbours in the electrode coordinate set. Their physical distance ranks only 127th among all electrode pairs. The issue is that transformed tensor-product P-spline basis columns have very similar response shapes on the discrete EEG grid, with maxima around the anterior scalp boundary.

Therefore:

- Do not merge/delete `FPZ` and `AF1` as a local electrode-level patch.
- If `B_ps k=8` still has poor sampler geometry, inspect the P-spline basis/penalty/boundary behavior rather than deleting electrodes.
- The graph-Laplacian/electrode-graph route remains the cleaner alternative because it encodes scalp topology directly rather than relying on a projected-plane tensor spline boundary.

## Current k=8 job

`B_ps k=8` job `39983288` is still running. At last check it was in warmup around 200 / 3000 iterations.
