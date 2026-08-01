# Half-Max BGAM Geometry Pilot Design

## Objective

Measure the separate computational effects of lowering the NUTS target
acceptance rate and replacing the current `t2()` P-spline construction with a
better-conditioned `te()` P-spline construction. This is a computational
screening experiment, not part of the three-radius scientific comparison.

## Input and fixed model settings

- Dataset: `data/montecarlo_rhalfmax/sim_001.rds`
- Ground truth: `data/montecarlo_rhalfmax/gt_001.rds`
- Scenario: `high`
- Radius rule: `half-max`
- Participants: 20
- Seed: 1
- Spatial basis dimension: `k_xy = 16`
- Temporal basis dimension: `k_time = 12`
- Metric: `diag_e`
- Chains: 1
- Threads per chain: 4
- Iterations: 1250
- Warmup: 1000
- Maximum treedepth: 10
- Progress refresh: 25
- Save warmup: true
- Likelihood, random-effects structure, measurement-error treatment, and
  prior families: identical across all arms

The three arms must run sequentially inside one Slurm allocation on the same
node. Model compilation happens before the timed fitting section so compilation
time is not attributed to sampling.

## Experimental arms

| Arm | Smooth construction | `adapt_delta` | Isolated intervention |
|---|---|---:|---|
| A | `t2(..., bs = c("ps", "ps"), full = FALSE)` | 0.80 | Same-node short baseline |
| B | `t2(..., bs = c("ps", "ps"), full = FALSE)` | 0.75 | Lower target acceptance |
| C | `te(..., bs = c("ps", "ps"), np = FALSE)` | 0.80 | Better-conditioned tensor basis |

All other settings are held fixed. Arms A and B have identical parameter
dimensions and use the same seed. Arm C changes the tensor penalty structure,
so the same seed cannot guarantee identical parameter-level initialization.

## Outputs

Each arm writes to its own directory beneath
`experiments/geometry_pilot_halfmax/results/` and must refuse to overwrite an
existing completed result. Each arm retains:

- the serialized `brmsfit`;
- raw CmdStan CSV files;
- saved warmup and adapted metric artifacts;
- model configuration;
- wall-clock timing split into compilation, warmup, and sampling where
  available;
- NUTS diagnostic summary;
- a compact prediction summary on a fixed electrode-time grid.

A combined comparison table reports:

- warmup seconds per iteration;
- sampling seconds per iteration;
- step size;
- median, 90th percentile, and maximum `n_leapfrog`;
- count and rate of `treedepth__ == 10`;
- divergences;
- mean acceptance statistic;
- posterior mean at CZ and 250 ms;
- posterior-mean map difference relative to Arm A.

R-hat is not used because this screening pilot has one chain. ESS estimates may
be retained as descriptive output but are not used as the primary selection
criterion with only 250 post-warmup draws.

## Interpretation rule

This pilot is descriptive and does not automatically select a winner. Any arm
with sampling divergences is excluded from immediate four-chain confirmation.
For the remaining arms, the report presents sampling seconds per iteration,
leapfrog summaries, treedepth saturation, and fixed-grid posterior-mean
differences side by side; the evidence is reviewed before any full rerun.

If Arms B and C both show useful computational changes without divergences or
an unacceptable prediction change, a later fourth screening arm may combine
`te(..., np = FALSE)` with `adapt_delta = 0.75` to measure their interaction.
The combined arm is not submitted automatically and is outside the initial
three-arm job.

## Operational safeguards

- The pilot uses an experiment-specific output and cache directory.
- Existing paired-radius fits and diagnostics are read-only.
- The currently running `dense_halfmax` job is neither cancelled nor modified.
- Any non-zero arm exit stops the sequential job and preserves logs and raw
  artifacts for diagnosis.
- Before submission, the R script is parsed locally and remotely, the Slurm
  script passes `bash -n`, and input files are checked on TU Dortmund.
