# Project Guide

## Formal analysis contract

The thesis analysis uses one configuration defined by
`code/config/paper_config.R`:

- scenario: `high`;
- participants: 20;
- trials per participant and condition: 20;
- paired replicates: 20;
- signal mode: `hard_signal`;
- radius rules: `1sd`, `half-max`, and `2sd`;
- active electrodes: 5, 9, and 21;
- reference onset and offset: 0.160 and 0.342 s;
- BGAM margins: P-splines;
- basis dimensions: `k_xy = 16`, `k_t = 12`;
- measurement error: `se(eeg_sem, sigma = TRUE)`;
- four chains, 3,000 iterations, and 1,000 warmup iterations.

Radius changes do not change noise, participant draws, temporal envelopes, or
untruncated spatial weights within a replicate.

## Read the formal code in this order

1. `code/config/paper_config.R`
   Defines the paper settings without source-time side effects.

2. `code/utils.R`
   Provides project paths, resources, thresholds, trial aggregation, and
   shared evaluation helpers.

3. `code/pipeline/simulation.R`
   Generates a common stochastic base and applies hard spatial radii. Spatial
   weights outside a radius are exactly zero.

4. `code/pipeline/generate_montecarlo.R`
   Generates all three radius datasets from one base for each seed.

5. `code/truth/ground_truth.R`
   Provides population-level ground-truth reconstruction used by paired-radius
   diagnostics.

6. `code/pipeline/fit_bgam.R`
   Defines the common BGAM fitting function.

7. `code/pipeline/fit_paired_radius.R`
   Applies the same model configuration to a selected radius and replicate.

8. `code/pipeline/compare_methods.R`
   Implements BGAM, pointwise, multiple-testing, cluster-mass, and TFCE
   comparison helpers.

9. `code/pipeline/aggregate_results.R`
   Aggregates completed per-replicate method results.

## Truth interpretation

Spatial support is exact: outside the selected radius, injected signal is zero.
The active spatial set is derived from positive truncated spatial weights.

Temporal support is participant-specific. Each participant receives a finite,
zero-padded temporal waveform shifted by that participant's onset. There is no
common 160--342 ms clipping operation after shifting. The population mean can
therefore have small tails outside the fixed timing-scoring interval.

The fixed reference interval is represented at 100 Hz by times 0.16 through
0.34 s. It is a scoring definition, not a claim that population truth is
exactly zero at every other time.

The project distinguishes:

- theoretical population truth for cross-replicate recovery and map targets;
- realised finite-sample truth for conditional residual diagnostics.

## Supporting directories

- `code/reporting/`: thesis figures, slides, posterior maps, and topomaps.
- `code/diagnostics/`: fit diagnostics, residual checks, and methodological
  investigations.
- `code/sensitivity/`: prior, basis, threshold, alpha, SESOI, and Laplacian
  analyses.
- `code/pilot/`: reduced-scope paired generation, validation, and aggregation.
- `code/archive/`: superseded single-radius and one-off historical scripts.

Supporting scripts do not define the formal paper configuration.

## Slurm entry points

```text
slurm/formal/generate_paired_radius.sh
slurm/formal/fit_paired_radius_hpc2n.sh
slurm/formal/fit_paired_radius_dortmund.sh
slurm/formal/aggregate_paired_radius.sh
slurm/pilot/run_paired_radius_pilot.sh
```

The formal fitting scripts map their 60 array tasks to three radius rules and
20 replicate seeds. Umea and Dortmund launchers differ only in scheduler,
module, memory, and filesystem configuration.

## Local-only outputs

```text
data/montecarlo_r1sd/
data/montecarlo_rhalfmax/
data/montecarlo_r2sd/
results/
logs/
cmdstan_cache/
```

These directories are not part of the GitHub source tree. Repository cleanup
must never delete or rewrite their contents.

## Verification commands

Parse all R scripts:

```bash
find code -name '*.R' -print0 | \
  xargs -0 -n1 Rscript -e 'parse(file=commandArgs(TRUE)[1])'
```

Check Slurm syntax:

```bash
for script in slurm/formal/*.sh slurm/pilot/*.sh; do
  bash -n "$script" || exit 1
done
```

Run regression tests:

```bash
for test_file in tests/test_*.R; do
  PROJECT_ROOT="$PWD" Rscript --vanilla "$test_file" || exit 1
done
```
