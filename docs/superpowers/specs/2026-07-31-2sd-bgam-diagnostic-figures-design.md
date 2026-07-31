# 2 SD BGAM Diagnostic Figures Design

## Objective

Generate three diagnostic PNG figures for the existing 2 SD BGAM fit. Use the
data-generating superpopulation mean as truth and do not run cluster-mass or
TFCE.

## Shared Inputs and Conventions

- Fit: `results/paired_radius_bgam_pilot/2sd/fit_seed001.rds`.
- Data: `data/montecarlo_r2sd/sim_001.rds`.
- Truth: reconstructed with `paired_radius_superpopulation_truth()`.
- Display channels: `CZ`, `FCZ`, `CPZ`, `C2`, `PZ`, and `FP1`.
- Fit and posterior-probability landmarks: 200, 250, and 300 ms, matching the
  saved temporal active-window onset, peak, and offset.
- Topomap times: 160, 250, and 340 ms.
- Practical-effect probability: `P(cond2 - cond1 > 0.1 uV | data)`.
- Detection threshold: posterior odds greater than 10, equivalent to
  probability greater than `10 / 11 = 0.9091`.

## Figure 1: BGAM Fit

Create six electrode facets over time. Each facet contains:

- participant-level empirical condition differences as faint grey curves;
- empirical grand-mean difference as a dotted black curve;
- BGAM posterior mean as a solid coloured curve;
- BGAM 95% credible interval as a ribbon;
- superpopulation truth as a distinct dashed curve;
- vertical markers at 200, 250, and 300 ms and a horizontal zero line.

This figure evaluates whether the fitted effect follows the superpopulation
effect in central, boundary, and inactive spatial locations.

## Figure 2: Posterior Probability

Create the same six electrode facets. Plot
`P(cond2 - cond1 > 0.1 uV | data)` over time, with horizontal reference lines
at 0.5 and 0.9091. Shade contiguous periods exceeding 0.9091 in green. Mark
200, 250, and 300 ms vertically.

## Figure 3: Topomap Diagnostic

Create one multi-panel figure with columns at 160, 250, and 340 ms and rows:

1. ground truth (superpopulation);
2. empirical condition difference;
3. noise residual: empirical minus ground truth;
4. BGAM estimate;
5. empirical minus BGAM.

Use one shared symmetric colour scale across all five rows and all three times,
with limits set by the largest absolute interpolated value over all 15 panels.
Mark active and inactive electrode locations with different point shapes.

## Outputs and Verification

Write PNG files under `results/paired_radius_bgam_pilot/2sd/figures/`:

- `bgam_fit_superpopulation.png`;
- `posterior_probability_sesoi.png`;
- `topomap_diagnostics_superpopulation.png`.

Also save the plotted numerical data as RDS files for reproducibility. Tests
must verify the probability definition, 0.9091 threshold, superpopulation
truth use, requested channels, topomap times, unified colour limit, non-empty
output files, and valid PNG dimensions.
