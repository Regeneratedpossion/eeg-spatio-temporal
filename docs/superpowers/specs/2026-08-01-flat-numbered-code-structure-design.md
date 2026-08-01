# Flat Numbered Code Structure

## Goal

Restore the familiar numbered layout under `code/`, while improving it so
that every executable script has one unique, continuous number. Preserve the
latest implementation from GitHub `main`; this is a path and naming change,
not a rollback to older code.

## Scope

- Flatten active scripts from `pipeline/`, `pilot/`, `diagnostics/`,
  `reporting/`, `sensitivity/`, `truth/`, and `config/` into `code/`.
- Keep `code/archive/` for obsolete and one-time scripts.
- Keep sourced helper modules unnumbered.
- Update all internal `source()` calls, Slurm entrypoints, tests, README, and
  project-guide references.
- Do not change simulation logic, model formulas, priors, sampler settings,
  SESOI, posterior odds, detection rules, cluster-mass, TFCE, or metrics.
- Do not generate data, refit models, or create result files.

## Final Layout

```text
code/
├── 00_check_environment.R
├── 01_simulate_eeg_data.R
├── 02_generate_montecarlo.R
├── 03_fit_bgam.R
├── 04_compare_methods.R
├── 05_aggregate_results.R
├── 06_generate_paired_radius_pilot.R
├── 07_validate_paired_radius_pilot.R
├── 08_fit_paired_radius_bgam.R
├── 09_diagnose_paired_radius_bgam.R
├── 10_aggregate_paired_radius_bgam.R
├── 11_run_bgam_diagnostics.R
├── 12_validate_heterogeneity_bgam.R
├── 13_compute_delta_sesoi.R
├── 14_plot_sesoi_posterior.R
├── 15_diagnose_practical_null_fpr.R
├── 16_plot_paired_radius_truth.R
├── 17_plot_2sd_bgam_diagnostics.R
├── 18_plot_simulated_eeg.R
├── 19_plot_bgam_fit.R
├── 20_plot_posterior_probability.R
├── 21_plot_topomaps.R
├── 22_prior_predictive_check.R
├── 23_basis_screening.R
├── 24_basis_screening_full.R
├── 25_laplacian_sensitivity.R
├── 26_posterior_odds_sensitivity.R
├── 27_significance_alpha_sensitivity.R
├── utils.R
├── paper_config.R
├── ground_truth.R
├── laplacian_config.R
├── laplacian_fit_basis_scores.R
├── laplacian_reconstruct_effect.R
├── laplacian_diagnose_hard.R
├── laplacian_diagnose_continuous.R
└── archive/
```

## Numbering Rules

- `00-05`: main Monte Carlo pipeline.
- `06-10`: paired-radius pilot pipeline.
- `11-15`: BGAM and SESOI diagnostics.
- `16-21`: thesis and slide figures.
- `22-27`: prior, basis, Laplacian, and threshold sensitivity analyses.
- Only directly executable scripts receive numbers.
- Helper modules remain unnumbered and are sourced by numbered scripts.
- Files in `archive/` do not participate in the active numbering sequence.

## Archive Policy

Move superseded Monte Carlo drivers, one-time ground-truth migration scripts,
and unused utilities to `code/archive/`. Preserve their existing filenames
inside the archive so their history remains recognisable. Do not archive an
active Slurm entrypoint or a file required by an active numbered script.

## Implementation Safety

The current local checkout contains unrelated uncommitted changes and differs
from GitHub `main`. Perform the reorganisation in an isolated worktree based
on the latest `origin/main`. Do not reset, delete, overwrite, or include the
user's existing local changes.

Use Git-aware moves so history remains traceable. After moving files, update
every path reference before deleting empty category directories.

## Verification

1. Confirm all expected numbered scripts exist exactly once and numbers
   `00-27` are continuous.
2. Confirm no active reference points to the removed category directories.
3. Parse every active R file without executing analyses.
4. Run path- and orchestration-related tests.
5. Run the remaining non-model unit tests that do not generate results or fit
   BGAMs.
6. Confirm the diff contains only moves, path updates, and documentation
   updates; statistical constants and model code must remain unchanged.

## Out of Scope

The separate inconsistency between the fixed `0.10` microvolt ROPE and the
baseline-calibrated high-scenario SESOI of approximately `0.033` microvolts is
not corrected in this structural change. It requires a separate reviewed
methods change.
