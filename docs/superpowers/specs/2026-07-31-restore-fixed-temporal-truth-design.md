# Restore Fixed Temporal Ground Truth

## Objective

Restore the pre-specified temporal ground-truth boundaries for every
heterogeneity scenario: onset 0.160 s and offset 0.342 s. On the 100 Hz
analysis grid, active time points are 0.16 through 0.34 s.

## Scope

- `build_ground_truth()` must derive `true_onset`, `true_offset`,
  `temporal_active_times`, and `active_cells` from the pre-specified
  parameters.
- The `mean(reference_waveform) + sd(reference_waveform)` threshold must no
  longer define temporal truth.
- Existing paired-radius `gt_001.rds` files must receive the corrected label
  fields without changing `continuous_effect` or any simulated EEG file.
- Spatial weights, active electrodes, model fits, estimation rules, and model
  outputs remain unchanged.

## Verification

- A regression test must fail against the current high-heterogeneity branch
  and pass after the fix.
- All three ground-truth files must report onset 0.160, offset 0.342, and the
  literal grid `seq(0.16, 0.34, 0.01)`.
- Before/after serialized hashes of `continuous_effect` must match.
- Hashes of all three `sim_001.rds` files must remain unchanged.

