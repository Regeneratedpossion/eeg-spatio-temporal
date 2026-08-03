# Half-Max Posterior Geometry Figure Design

## Purpose

Create one compact diagnostic figure for discussion with the supervisors. The
figure should show whether the strongest pairwise marginal posterior relation
among condition-2 BGAM smooth coefficients resembles a straight tilted ellipse,
as motivated by Betancourt's discussion of matching kinetic energy to posterior
geometry.

## Input

- Full half-max fit:
  `results/paired_radius_bgam_pilot/half-max/fit_seed001.rds`
- Use all four chains and all 8,000 post-warmup draws.
- Select the condition-2 smooth-coefficient pair with the largest absolute
  Pearson correlation from the complete fit. For the current file this is
  `s_t2xprojyprojtimeconditioncond2_1[89]` and
  `s_t2xprojyprojtimeconditioncond2_3[18]`.

## Figure

- Standardize both displayed coefficients so the axes are comparable.
- Plot posterior draws as small, translucent points.
- Use one color per chain and include a compact chain legend.
- Overlay a 95% covariance ellipse.
- Use equal coordinate scaling.
- Annotate the observed correlation and the following full-fit diagnostics:
  zero divergences, median `n_leapfrog = 1023`, and max-treedepth saturation
  `3802 / 8000 = 47.5%`.
- Title: `Half-max BGAM: posterior geometry of the most correlated smooth coefficients`.
- Subtitle: `A straight tilted ellipse is consistent with a Euclidean metric mismatch (Betancourt, 2018, Sections 4.2 and 6.1)`.

## Outputs

- PNG for quick viewing and email.
- PDF for slides or supervisor discussion.
- Save both under `results/diagnostics/halfmax_posterior_geometry/`.

## Validation

- Confirm that exactly 8,000 post-warmup draws are plotted.
- Confirm that all four chains are present.
- Recompute the selected pair and correlation from the fit rather than
  hard-coding the numerical correlation.
- Confirm that the generated PNG and PDF are non-empty.

## Scope

This is a pairwise marginal diagnostic, not a visualization of the entire
high-dimensional posterior and not evidence by itself that a low-rank metric
will solve the sampling problem.
