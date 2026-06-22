# Low-rank Laplacian basis model

This folder contains the current low-rank spatial-basis approximation for the
spatiotemporal Bayesian GAMM.

## Model idea

Instead of fitting one full sensor x time BGAMM directly, the 61-channel scalp
map is first projected onto graph-Laplacian spatial basis functions derived from
the sensor adjacency matrix.

For each participant x condition x time cell:

1. Build a graph Laplacian from the sensor adjacency matrix.
2. Take the first 20 smooth eigenvectors, including the constant/global spatial
   component.
3. Project each 61-channel scalp map onto these spatial basis vectors.
4. Fit one Bayesian temporal GAMM to each basis score using `brms`.
5. Reconstruct the posterior sensor x time effect map by summing the posterior
   basis-score effects multiplied by their spatial basis vectors.

The workflow is split into:

```text
00_config.R                 fixed K=20 configuration
01_fit_basis_scores.R        fit temporal GAMMs to Laplacian basis scores
02_reconstruct_effect.R      reconstruct sensor x time posterior effects
03_diagnose_hard.R           hard-threshold active/inactive diagnostics
04_diagnose_continuous.R     continuous Gaussian-field diagnostics
```

## Current configuration

```text
Dataset: P20 gaussian hardlabel
Spatial basis: graph-Laplacian eigenvectors
K: 20 + constant/global component
Temporal model: Bayesian temporal GAMM on each basis score
iter: 3000
warmup: 1000
adapt_delta: 0.999
max_treedepth: 13
```

## Current diagnostic results

### Hard-threshold diagnostic

```text
convergence_pass = FALSE
temporal_pass = TRUE
spatial_top_pass = TRUE
peak_ci_covers_truth = FALSE
mean_ci_covers_truth = FALSE
spatial_F1 = 0.340
spatiotemporal_F1 = 0.545
peak_recovery = 81.4%
mean_recovery = 95.1%
```

`convergence_pass = FALSE` is caused by residual divergent transitions in two
basis models:

```text
basis 11: 2 divergent transitions
basis 19: 14 divergent transitions
```

However, all basis models have:

```text
Rhat < 1.01
ESS > 900
treedepth warnings = 0
```

So the model is not globally failing, but the convergence is not completely
clean under the strict criterion that requires zero divergent transitions.

### Continuous Gaussian-field diagnostic

```text
spatiotemporal_cor_all = 0.962
spatiotemporal_cor_effect_window = 0.964
peak_spatial_cor = 0.976
peak_recovery = 81.4%
truth_weighted_mean_recovery = 104.5%
top25_mean_recovery = 102.7%
```

## Interpretation

The continuous diagnostic is the more appropriate primary evaluation for this
simulation, because the simulated spatial effect is a smooth Gaussian field
rather than a binary active/inactive mask. Under this evaluation, the low-rank
Laplacian model recovers the main spatiotemporal structure well: the temporal
effect window is detected, the strongest spatial effect is located around CZ,
and the reconstructed sensor x time effect map is highly correlated with the
true continuous effect.

Adding the constant/global spatial component also substantially improves
amplitude recovery. The peak effect is recovered at about 81% of the true value,
while the truth-weighted and top-25% mean recoveries are close to 100%. This
suggests that the remaining approximation error is no longer dominated by the
strong amplitude shrinkage seen in the smaller-basis version.

The main limitation is that two basis-score models still show divergent
transitions: basis 11 has 2 divergences and basis 19 has 14. Therefore, the
model should not be described as fully clean under strict convergence criteria,
even though Rhat, ESS, and treedepth diagnostics are otherwise acceptable. The
hard-threshold F1 scores should be treated as secondary diagnostics because
they impose a binary active/inactive evaluation on a continuous Gaussian signal.
