# Method choice rationale

A short record of *why* each modelling decision was made, useful when
writing methods sections and answering reviewer questions.

---

## Why BGAM (not pure cluster permutation)

- Cluster permutation gives only cluster-level inference; cannot
  pinpoint precise onset/offset or amplitude.
- Bayesian GAM gives full posterior over (electrode, time) — every
  point has a credible interval and a posterior odds ratio.
- Reference: Maris & Oostenveld (2007) acknowledge this limitation;
  Nalborczyk & Bürkner (2025) demonstrate the BGAM solution in 1D.

## Why 3D (xproj, yproj, time) joint, not 1D × N channels

- Pure 1D-per-channel approach (Nalborczyk's recipe applied to all
  electrodes) loses spatial smoothness information.
- Joint modelling lets neighbouring electrodes share information
  (shrinkage), useful when n_participants is small.
- Trade-off: spatial smoothing causes amplitude shrinkage and bump
  diffusion — discussed in progress_report.md §3.3, §3.6.

## Why `t2()` (not `te()` or `s+s+ti()`)

- `te()`: single global penalty — too restrictive for 3D ST data.
- `t2()`: marginal penalties per dimension — flexible, identifiable.
- `s(x,y) + s(t) + ti(x,y,t)`: would give cleaner ANOVA decomposition
  but **brms does not support `ti()`** (only `t2()`). Tried, failed.
- Therefore: `t2()` is the only viable joint smooth in brms.

## Why `bs = c("tp", "tp")` (not `"gp"` or `"cr"`)

- `tp`: thin-plate basis, mgcv default, well-tested, no extra params.
- `gp`: Gaussian process basis — would give Matérn-style covariance,
  potentially more appropriate for EEG. Pending exploration.
- `cr`: cubic regression spline — only valid for 1D, can't be used in
  the 2D spatial component of `t2()`.

## Why `k = c(10, 12)` for `t2()`

- Spatial dim k = 10: with 61 electrodes, k ≤ n_unique_locations / 5
  is a common heuristic. 10 is conservative.
- Time dim k = 12: time has 51 sample points; k = 12 is loose enough
  to capture an ERP-shape without overfitting noise.
- Tried k = 20 (`v4_k20`) — gave very similar results, slightly slower.
  k = 12 chosen as the sweet spot for ESS efficiency.

## Why `(1 + condition | participant)` random slopes

- Each participant has their own baseline EEG level (intercept) and
  may respond differently to the condition (slope).
- Without random slope, condition effect would be over-confident
  (Type-I inflation in classical terms).
- Reference: Barr et al. (2013) "Keep It Maximal", standard practice
  in psycholinguistic / cognitive-neuroscience GLMM literature.

## Why `eeg_mean | se(eeg_sd)` (measurement-error)

- Aggregating trials to per-cell mean discards information about
  trial-to-trial variability. Passing the per-cell SD as known
  measurement error preserves this.
- Allows fitting to ~ N_ppt × N_ch × N_time rows instead of trial-level
  rows — much faster, similar inference (Nalborczyk's "meta-analytic"
  trick).
- Side effect: model has no `sigma` parameter to estimate (fixed by
  `eeg_sd`). Hence prior on `class = sigma` would error (we removed it).

## Why sum-to-zero contrast on `condition`

- Default treatment contrast: Intercept = mean of cond1.
- Sum-to-zero (-0.5, +0.5): Intercept = grand mean.
- Makes parametric `conditioncond2` interpretable as the global mean
  difference (cond2 − cond1).
- Required to match the paper's set-up.

## Why `adapt_delta = 0.99`, `max_treedepth = 13`

- Earlier `v3_fast3` used `max_treedepth = 10` and saturated on
  every draw (8000/8000 hits) — sampler was step-size-bound.
- Bumping to 13 → 0 hits in `v4_k12`. Sampler now operates freely.
- `adapt_delta = 0.99` is high but `t2` smooths can have funnels;
  cheap insurance.

## Why `ndraws(meta_gam) = 8000` (4 chains × 2000 sampling)

- 8000 posterior draws gives effective sample size ~ 2000+ per
  parameter, sufficient for stable interval estimation.
- More draws would help only marginally; runtime cost is large.

## Why `n_permutations = 4096` for cluster comparison (planned)

- Matches Nalborczyk paper exactly (their value: 2^12 = 4096).
- Allows ~ 0.0002 minimum p-value resolution.
- Sufficient for one-sample test; bigger would be overkill.
