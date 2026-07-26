# Continuous Ground-Truth Diagnostics Design

## Objective

Ensure that diagnostic plots and amplitude-recovery metrics use the exact
population-level continuous effect saved by the data generator. Diagnostic
code must not reconstruct a substitute Gaussian curve from hard-label onset
and offset boundaries.

## Scope

This change covers two diagnostic paths:

1. `code/08_slide_topomap_n20.R`, which produces the ground-truth and
   empirical-minus-ground-truth topomap layers.
2. `code/02_diagnostics.R`, which computes peak and window-mean amplitude
   recovery.

Regenerating figures, updating thesis text, changing the data-generating
process, and recomputing the historical summary table are outside this code
change.

## Design

Add small ground-truth accessors to `code/utils.R`:

- Validate that `ground_truth$continuous_effect` exists and contains
  `channel`, `time`, and `effect`.
- Select exact channel-time effects from the saved analysis grid.
- Fail with an explicit error when exact continuous truth is unavailable.

The topomap script will obtain its truth layer by joining the requested
channel-time grid to the saved continuous effect. It will no longer infer a
temporal waveform from `true_onset` and `true_offset`.

The amplitude diagnostic will compute:

- peak truth at the same seed-electrode and peak-time grid cell used for the
  posterior estimate;
- window-mean truth over exactly the same discrete time points used for the
  posterior window mean.

It will no longer use the analytic `truth_window_mean()` approximation or
assume that `bump_amplitude` is the group-level peak after onset averaging.

## Legacy-data policy

Ground-truth files without `continuous_effect` will stop with a clear error.
There will be no fallback to the historical Gaussian approximation, because
that would silently reproduce the metric mismatch. Historical outputs remain
unchanged and can be labelled as legacy outside this code path.

## Error handling

Diagnostics stop when:

- `continuous_effect` is missing;
- required columns are missing;
- a requested channel-time cell is absent or duplicated;
- a requested time is not on the saved analysis grid.

Error messages identify the missing requirement and explain that the dataset
must be regenerated with the current generator.

## Tests

Tests will be written before implementation and will verify:

1. exact lookup returns saved values, including non-zero boundary effects;
2. missing `continuous_effect` fails explicitly;
3. missing or duplicate requested cells fail;
4. peak and window means use the requested discrete grid cells;
5. the topomap truth construction no longer depends on threshold-defined
   onset and offset values.

## Acceptance criteria

- No diagnostic truth value is reconstructed from `true_onset` and
  `true_offset`.
- For the high-heterogeneity seed-1 truth, CZ effects at 0.16, 0.25, and
  0.34 seconds equal the saved values (approximately 0.1211, 1.4574, and
  0.1211 microvolts).
- Legacy truth objects fail rather than returning an approximate recovery.
- New regression tests and existing relevant checks pass.
