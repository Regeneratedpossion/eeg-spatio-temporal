# Paired-Radius Superpopulation Truth Design

## Objective

Diagnose whether the fitted BGAM recovers the data-generating
superpopulation mean effect, rather than the realised mean effect of the 20
participants drawn for seed 1.

## Estimand

For channel `c` and analysis time `t`, define

```text
theta_pop(c, t) = spatial_weight(c) * E_A[A]
                  * mean_tau(f0(t - tau)),
```

where `tau` is discrete uniform on the prespecified onset support
`140, 142, ..., 180` ms, `E_A[A]` is `bump_amplitude = 1.5`, and `f0` is the
normalised ERP template used by the generator. The hard-radius spatial weight
is zero outside the active set, so population truth is exactly zero there.

The calculation must not use seed 1 participant amplitudes, participant
onsets, or the existing `continuous_effect`, all of which encode the realised
20-participant sample. The existing `population_mean_waveform` field is also
not used because the current generator populated it with the realised sample
waveform.

## Implementation

Add a focused helper that reconstructs the superpopulation temporal waveform
from the original ERP template, the saved onset support, downsampling factor,
and expected amplitude. Multiply it by the saved deterministic hard-truncated
spatial weights to obtain the complete channel-by-time truth surface.

Update the paired-radius diagnostic to use this surface for:

- `truth_cz_250` and `peak_ratio`;
- projection recovery;
- truth values attached to false-positive locations.

Retain the realised CZ truth as `realized_truth_cz_250` for auditability and
write `truth_estimand = "superpopulation"` in the diagnostic summary.

Do not change the active set, SESOI, posterior odds, detection masks,
cluster-mass, TFCE, FPR, leakage, convergence diagnostics, EDF companion fit,
or the empirical-minus-BGAM topomap.

## Verification

Tests must fail before implementation and then establish that:

1. the 2 SD CZ truth at 250 ms is approximately `1.457449` microvolts;
2. inactive-electrode truth is exactly zero;
3. changing realised participant draws or `continuous_effect` does not change
   superpopulation truth;
4. the helper rejects malformed or incomplete inputs.

After the tests pass, rerun only the existing 2 SD fitted model and verify all
expected diagnostic outputs, their estimand label, and the unchanged
non-truth-dependent metrics.
