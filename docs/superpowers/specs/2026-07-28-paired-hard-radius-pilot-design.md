# Paired Hard-Radius Pilot Design

## Goal

Move spatial hard thresholding into EEG data generation and produce a paired
pilot with three active radii. The three radius variants for each replicate
must share exactly the same noise draws, participant heterogeneity, temporal
signals, and untruncated Gaussian spatial weights. Only the deterministic
spatial truncation radius may differ.

The pilot uses:

- `n_ppt = 20`;
- `scenario = "high"`;
- `n_trials = 20` per participant and condition, inherited from the existing
  generator;
- five replicates with seeds `1:5`;
- `signal_mode = "hard_signal"` only.

This change is confined to simulation generation, ground truth, pilot
validation, and pilot truth figures. Model formulas, priors, SESOI, posterior
odds, cluster-mass, TFCE, detection rules, and evaluation algorithms remain
unchanged.

## Radius Rules

The raw circular Gaussian spatial weight remains

\[
w_{\mathrm{raw}}(d) = \exp\left(-\frac{d^2}{2\sigma_{\mathrm{space}}^2}\right).
\]

Each final spatial weight is derived using an explicit `r_active` parameter:

```r
spatial_weight = if_else(
    spatial_distance <= r_active,
    spatial_weight_raw,
    0
)
```

The three pilot rules are:

| Rule | `radius_rule` | `r_active` | Expected active electrodes |
|---|---|---:|---:|
| One spatial SD | `1sd` | `sigma_space` | 5 |
| Half-maximum contour | `half-max` | `sigma_space * sqrt(2 * log(2))` | 9 |
| Two spatial SDs | `2sd` | `sigma_space * 2` | 21 |

The old `active_threshold = mean(spatial_weight_raw) +
sd(spatial_weight_raw)` quantity is removed. Active electrodes have no
independent threshold parameter and are always derived from the truncated
weight:

```r
active_electrodes = sensors |>
    filter(spatial_weight > 0) |>
    pull(channel)
```

## Paired Generation Architecture

Each replicate calls `set.seed(seed)` exactly once and then creates one shared
random base. The base contains both conditions' noise draws, participant onset
latencies, participant amplitudes, participant noise scales, temporal signal
templates, sensor distances, and `spatial_weight_raw`.

Radius materialization is deterministic. For each radius rule it:

1. derives `spatial_weight` from the same sensor base;
2. adds the radius-specific signal to the shared condition-2 noise;
3. leaves condition-1 and condition-2 noise unchanged;
4. downsamples into the existing final `sim_data` schema;
5. builds radius-specific continuous and binary ground truth.

The base is generated once and materialized three times. The implementation
must not run the complete random generator three times, reset the seed for each
radius, or depend on saved and restored `.Random.seed` states.

The final saved `sim_data` objects retain the existing public columns so model
and evaluation code do not need changes. Any internal signal template or noise
component used to materialize the variants is not added to the saved analysis
schema.

## Generator Interfaces and Defaults

`code/00_generate_data.R` separates stochastic base construction from
deterministic radius materialization. `r_active` and `radius_rule` are explicit
arguments at the materialization/truth boundary. `signal_mode` defaults to and
is validated as `"hard_signal"`.

`code/00_generate_montecarlo.R` defines the three radius configurations and
writes each paired variant to its designated directory:

- `data/montecarlo_r1sd/`;
- `data/montecarlo_rhalfmax/`;
- `data/montecarlo_r2sd/`.

All Slurm generator assignments of `MC_SIGNAL_MODE=continuous_label` are
changed to `MC_SIGNAL_MODE=hard_signal`. Existing unrelated local edits in
those scripts are preserved. The pilot runs locally through a dedicated
five-replicate entry point and does not submit or run a full Monte Carlo job.

Each output directory receives a scenario manifest that records at least the
scenario, participant and trial counts, fixed signal mode, radius value,
radius rule, replicate seeds, and paired-design identity. Existing manifests
with incompatible settings cause an explicit failure rather than mixed output.

## Ground Truth

Every saved ground-truth object includes:

```r
signal_mode   = "hard_signal"
active_radius = r_active
radius_rule   = radius_rule
```

`continuous_effect` is built from the truncated `spatial_weight`, so every
inactive electrode has an exactly zero stored effect at all times. Both
`active_cells` and `active_electrodes` are derived using
`spatial_weight > 0`; neither uses a separate weight threshold.

## Automated Tests

Tests are written before production changes and exercise real generator
functions. They establish:

1. the three radius formulas and expected active counts of 5, 9, and 21 on
   `data/sim_data_pilot.rds` sensor coordinates;
2. exact zero spatial weights and exact zero continuous effects outside each
   radius;
3. ground-truth `signal_mode`, `active_radius`, and `radius_rule` fields;
4. a single stochastic base can be deterministically materialized into all
   three variants;
5. the shared per-row condition noise is bitwise identical across variants;
6. generator defaults and Slurm generator settings use `hard_signal`;
7. saved data retain the schema expected by current model and evaluation code.

Existing tests that contain the removed `active_threshold` fixture are updated
to express active status through `spatial_weight > 0`, without altering their
unrelated truth assertions.

## Pilot Validation

The pilot generates only seeds `1:5` for the single `high` heterogeneity
scenario. For every replicate and radius, validation reports:

```r
length(g$active_electrodes)
max(abs(ce$effect[inact]))
sum(ce$effect[inact] != 0)
```

Required results are respectively 5/9/21 active electrodes, exact maximum
inactive effect `0`, and zero nonzero inactive cells.

Pairing is checked at two levels:

1. the shared base noise arrays are compared with exact equality;
2. empirical-minus-truth residual arrays are compared across the three saved
   datasets and their maximum pairwise absolute differences are reported.

The first comparison is the structural hard guarantee that the stochastic
realization is paired. The second is the requested end-to-end diagnostic.
IEEE-754 aggregation and cancellation can theoretically introduce differences
at machine precision when computing `mean(noise + signal) - truth`; therefore
the report distinguishes exact stochastic equality from floating-point
roundoff and prints the observed maximum difference rather than silently
rounding it.

For each radius, the validation also reports the outermost active electrode or
tied electrodes and their raw Gaussian weights. Expected weights are
approximately 0.70, 0.51, and 0.19 for `1sd`, `half-max`, and `2sd` on the
current discrete sensor layout.

## Truth Topomaps

A pilot-only plotting script reads `gt_001.rds` from each radius directory and
plots `continuous_effect` at the analysis-grid time nearest 250 ms. It produces
one truth topomap per radius and one three-panel comparison, in both PNG and
PDF formats. The maps use the existing sensor coordinates, inverse-distance
interpolation, and scalp decorations, and require no fitted model.

The figures are written under a pilot-specific results directory rather than
overwriting existing slide or model-diagnostic figures.

## Error Handling and Scope Protection

Generation fails clearly when `r_active` is missing, non-finite, non-positive,
or when `radius_rule` is not one of `1sd`, `half-max`, and `2sd`. It also fails
if an output directory contains an incompatible manifest.

No full Monte Carlo run or Slurm submission is part of this work. No model,
prior, inferential threshold, SESOI, posterior-odds, cluster-mass, TFCE,
detection, or evaluation implementation is modified.
