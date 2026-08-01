# Low-rank Laplacian basis model

This folder contains the current low-rank spatial-basis approximation for the
spatiotemporal Bayesian GAMM.

## Model idea

Instead of fitting one full sensor x time BGAMM directly, the 61-channel scalp
map is first projected onto graph-Laplacian spatial basis functions derived from
the sensor adjacency matrix. This follows the broader idea of using low-rank
basis functions to represent spatial processes efficiently (Cressie et al.,
2022). Here, the sensors are treated as nodes in an adjacency graph. The
eigenvectors of the graph Laplacian then provide smooth spatial basis functions
over this sensor graph (Chung, 1997; Belkin & Niyogi, 2003).

This is currently an approximation: each basis score is fitted with a separate
brms model, so the reconstruction does not yet represent one joint posterior
over all spatial basis coefficients.

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
code/20_laplacian_sensitivity.R         numbered workflow entry point
code/laplacian/config.R                 fixed K=20 configuration
code/laplacian/fit_basis_scores.R       fit temporal GAMMs to basis scores
code/laplacian/reconstruct_effect.R     reconstruct sensor x time effects
code/laplacian/diagnose_hard.R          hard-threshold diagnostics
code/laplacian/diagnose_continuous.R    continuous-field diagnostics
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

Runtime on the cluster was about 2.5 hours for the 20 basis-score models. The
slowest single basis model took about 17 minutes. The diagnostic reconstruction
steps finished within a few minutes.

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

The continuous diagnostic is the better primary evaluation here, because the
simulated spatial effect is a smooth Gaussian field rather than a binary
active/inactive mask. This is also closer to spatial statistical models that
treat the signal as a continuous latent spatial field, often approximated using
structured low-rank or Markov representations (Lindgren et al., 2011; Cressie
et al., 2022). Under this evaluation, the low-rank Laplacian model recovers the
main spatiotemporal pattern well: the time window is detected, the strongest
spatial effect is around CZ, and the reconstructed sensor x time map is highly
correlated with the true effect.

Adding the constant/global spatial component improves amplitude recovery. The
peak effect is recovered at about 81% of the true value, and the truth-weighted
and top-25% mean recoveries are close to 100%.

The main limitation is that basis 11 and basis 19 still show divergent
transitions. Therefore, this result is promising but not fully clean under
strict convergence criteria. The hard-threshold F1 scores should be treated as
secondary diagnostics, because they turn a continuous Gaussian signal into a
binary active/inactive problem.

## References

- Cressie, N., Sainsbury-Dale, M., & Zammit-Mangion, A. (2022).
  *Basis-Function Models in Spatial Statistics*.
  https://arxiv.org/abs/2202.03660
- Lindgren, F., Rue, H., & Lindström, J. (2011). An explicit link between
  Gaussian fields and Gaussian Markov random fields: the stochastic partial
  differential equation approach. *Journal of the Royal Statistical Society:
  Series B*, 73(4), 423-498.
- Chung, F. R. K. (1997). *Spectral Graph Theory*. American Mathematical
  Society.
- Belkin, M., & Niyogi, P. (2003). *Laplacian Eigenmaps and Spectral Techniques
  for Embedding and Clustering*. Advances in Neural Information Processing
  Systems.
