# Flat Numbered Code Structure Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the category-based `code/` tree on GitHub `main` with one flat, uniquely and continuously numbered set of executable R scripts while preserving the latest statistical implementation.

**Architecture:** Perform Git-aware moves from the latest `origin/main` in an isolated worktree. Number only executable entrypoints (`00-27`); keep sourced helper modules unnumbered and preserve obsolete or one-time files under `code/archive/`. Update every caller and verify that the diff changes paths and documentation only, not model or evaluation behavior.

**Tech Stack:** Git, R 4.x, base R parsing, repository R tests, Slurm shell scripts.

## Global Constraints

- Base implementation work on the latest available `origin/main`, not the dirty local checkout.
- Preserve all current simulation, model, prior, sampler, SESOI, posterior-odds, cluster-mass, TFCE, detection, and metric logic byte-for-byte except where a path or filename literal must change.
- Do not generate data, fit models, execute posterior prediction, or create result files.
- Keep `code/archive/`; remove the empty category directories after all moves.
- Do not include or overwrite unrelated changes from `/Users/jiahuifan/Desktop/bgam_eeg_project`.
- The `0.10` versus `0.033` SESOI methods correction remains out of scope.

---

## Target File Map

| Current path on `origin/main` | Final path |
|---|---|
| `code/check_env.R` | `code/00_check_environment.R` |
| `code/pipeline/simulation.R` | `code/01_simulate_eeg_data.R` |
| `code/pipeline/generate_montecarlo.R` | `code/02_generate_montecarlo.R` |
| `code/pipeline/fit_bgam.R` | `code/03_fit_bgam.R` |
| `code/pipeline/compare_methods.R` | `code/04_compare_methods.R` |
| `code/pipeline/aggregate_results.R` | `code/05_aggregate_results.R` |
| `code/pilot/00_generate_paired_radius_pilot.R` | `code/06_generate_paired_radius_pilot.R` |
| `code/pilot/00_validate_paired_radius_pilot.R` | `code/07_validate_paired_radius_pilot.R` |
| `code/pipeline/fit_paired_radius.R` | `code/08_fit_paired_radius_bgam.R` |
| `code/diagnostics/35_diagnose_paired_radius_bgam_pilot.R` | `code/09_diagnose_paired_radius_bgam.R` |
| `code/pilot/36_aggregate_paired_radius_bgam_pilot.R` | `code/10_aggregate_paired_radius_bgam.R` |
| `code/diagnostics/02_diagnostics.R` | `code/11_run_bgam_diagnostics.R` |
| `code/diagnostics/29_validate_heterogeneity_bgam.R` | `code/12_validate_heterogeneity_bgam.R` |
| `code/diagnostics/28_delta_sesoi_probe.R` | `code/13_compute_delta_sesoi.R` |
| `code/reporting/31_plot_all_electrode_sesoi_posterior.R` | `code/14_plot_sesoi_posterior.R` |
| `code/diagnostics/19_practical_null_fpr_rebound.R` | `code/15_diagnose_practical_null_fpr.R` |
| `code/reporting/07_plot_paired_radius_truth_pilot.R` | `code/16_plot_paired_radius_truth.R` |
| `code/reporting/37_plot_2sd_bgam_diagnostics.R` | `code/17_plot_2sd_bgam_diagnostics.R` |
| `code/reporting/05_slide_simulated_eeg_n20.R` | `code/18_plot_simulated_eeg.R` |
| `code/reporting/06_slide_bgam_fit_n20.R` | `code/19_plot_bgam_fit.R` |
| `code/reporting/07_slide_posterior_probability_n20.R` | `code/20_plot_posterior_probability.R` |
| `code/reporting/08_slide_topomap_n20.R` | `code/21_plot_topomaps.R` |
| `code/sensitivity/10_prior_predictive_check.R` | `code/22_prior_predictive_check.R` |
| `code/sensitivity/11_basis_screening.R` | `code/23_basis_screening.R` |
| `code/sensitivity/18_basis_screening_full_ab_sigma_true_control.R` | `code/24_basis_screening_full.R` |
| `code/sensitivity/12_laplacian_basis_brms_fixed.R` | `code/25_laplacian_sensitivity.R` |
| `code/sensitivity/32_threshold_sensitivity.R` | `code/26_posterior_odds_sensitivity.R` |
| `code/sensitivity/33_alpha_sensitivity.R` | `code/27_significance_alpha_sensitivity.R` |
| `code/config/paper_config.R` | `code/paper_config.R` |
| `code/truth/ground_truth.R` | `code/ground_truth.R` |
| `code/sensitivity/laplacian/00_config.R` | `code/laplacian_config.R` |
| `code/sensitivity/laplacian/01_fit_basis_scores.R` | `code/laplacian_fit_basis_scores.R` |
| `code/sensitivity/laplacian/02_reconstruct_effect.R` | `code/laplacian_reconstruct_effect.R` |
| `code/sensitivity/laplacian/03_diagnose_hard.R` | `code/laplacian_diagnose_hard.R` |
| `code/sensitivity/laplacian/04_diagnose_continuous.R` | `code/laplacian_diagnose_continuous.R` |
| `code/truth/migrate_ground_truth.R` | `code/archive/migrate_ground_truth.R` |
| `code/sensitivity/laplacian/README.md` | `code/archive/laplacian_README.md` |

Existing `code/archive/09_run_montecarlo_sem.R` and
`code/archive/compute_trial_effective_sample_size.R` remain in place.

### Task 1: Lock Down the Final Inventory

**Files:**
- Modify: `tests/test_code_inventory.R`

**Interfaces:**
- Consumes: the exact final paths in the target file map.
- Produces: a deterministic pass/fail check for the flat active code inventory.

- [ ] **Step 1: Replace the expected inventory with the flat target list**

Define `expected_numbered` as the 28 paths from
`code/00_check_environment.R` through
`code/27_significance_alpha_sensitivity.R`. Define `expected_helpers` as:

```r
c(
    "code/utils.R",
    "code/paper_config.R",
    "code/ground_truth.R",
    "code/laplacian_config.R",
    "code/laplacian_fit_basis_scores.R",
    "code/laplacian_reconstruct_effect.R",
    "code/laplacian_diagnose_hard.R",
    "code/laplacian_diagnose_continuous.R"
)
```

Add these structural assertions:

```r
stopifnot(all(file.exists(c(expected_numbered, expected_helpers))))

numbered_files <- list.files(
    "code", pattern = "^[0-9]{2}_.+[.]R$", full.names = TRUE
)
observed_numbers <- as.integer(sub("^code/([0-9]{2})_.*$", "\\1", numbered_files))
stopifnot(identical(sort(observed_numbers), 0:27))

removed_dirs <- file.path(
    "code",
    c("pipeline", "pilot", "diagnostics", "reporting", "sensitivity", "truth", "config")
)
stopifnot(!any(dir.exists(removed_dirs)))
```

- [ ] **Step 2: Run the inventory test against the category-based tree**

Run: `Rscript --vanilla tests/test_code_inventory.R`

Expected: FAIL because the new flat paths do not exist yet.

- [ ] **Step 3: Commit the failing structural test**

```bash
git add tests/test_code_inventory.R
git commit -m "test: require flat numbered code inventory"
```

### Task 2: Move Entry Points and Helper Modules

**Files:**
- Move: every file in the target file map.
- Preserve: `code/archive/09_run_montecarlo_sem.R`
- Preserve: `code/archive/compute_trial_effective_sample_size.R`

**Interfaces:**
- Consumes: the current category-based `origin/main` tree.
- Produces: the complete flat file inventory expected by Task 1.

- [ ] **Step 1: Move the 28 executable scripts with Git-aware renames**

Run:

```bash
git mv code/check_env.R code/00_check_environment.R
git mv code/pipeline/simulation.R code/01_simulate_eeg_data.R
git mv code/pipeline/generate_montecarlo.R code/02_generate_montecarlo.R
git mv code/pipeline/fit_bgam.R code/03_fit_bgam.R
git mv code/pipeline/compare_methods.R code/04_compare_methods.R
git mv code/pipeline/aggregate_results.R code/05_aggregate_results.R
git mv code/pilot/00_generate_paired_radius_pilot.R code/06_generate_paired_radius_pilot.R
git mv code/pilot/00_validate_paired_radius_pilot.R code/07_validate_paired_radius_pilot.R
git mv code/pipeline/fit_paired_radius.R code/08_fit_paired_radius_bgam.R
git mv code/diagnostics/35_diagnose_paired_radius_bgam_pilot.R code/09_diagnose_paired_radius_bgam.R
git mv code/pilot/36_aggregate_paired_radius_bgam_pilot.R code/10_aggregate_paired_radius_bgam.R
git mv code/diagnostics/02_diagnostics.R code/11_run_bgam_diagnostics.R
git mv code/diagnostics/29_validate_heterogeneity_bgam.R code/12_validate_heterogeneity_bgam.R
git mv code/diagnostics/28_delta_sesoi_probe.R code/13_compute_delta_sesoi.R
git mv code/reporting/31_plot_all_electrode_sesoi_posterior.R code/14_plot_sesoi_posterior.R
git mv code/diagnostics/19_practical_null_fpr_rebound.R code/15_diagnose_practical_null_fpr.R
git mv code/reporting/07_plot_paired_radius_truth_pilot.R code/16_plot_paired_radius_truth.R
git mv code/reporting/37_plot_2sd_bgam_diagnostics.R code/17_plot_2sd_bgam_diagnostics.R
git mv code/reporting/05_slide_simulated_eeg_n20.R code/18_plot_simulated_eeg.R
git mv code/reporting/06_slide_bgam_fit_n20.R code/19_plot_bgam_fit.R
git mv code/reporting/07_slide_posterior_probability_n20.R code/20_plot_posterior_probability.R
git mv code/reporting/08_slide_topomap_n20.R code/21_plot_topomaps.R
git mv code/sensitivity/10_prior_predictive_check.R code/22_prior_predictive_check.R
git mv code/sensitivity/11_basis_screening.R code/23_basis_screening.R
git mv code/sensitivity/18_basis_screening_full_ab_sigma_true_control.R code/24_basis_screening_full.R
git mv code/sensitivity/12_laplacian_basis_brms_fixed.R code/25_laplacian_sensitivity.R
git mv code/sensitivity/32_threshold_sensitivity.R code/26_posterior_odds_sensitivity.R
git mv code/sensitivity/33_alpha_sensitivity.R code/27_significance_alpha_sensitivity.R
```

Do not edit file contents during this step.

- [ ] **Step 2: Move the eight helper modules with Git-aware renames**

Run:

```bash
git mv code/config/paper_config.R code/paper_config.R
git mv code/truth/ground_truth.R code/ground_truth.R
git mv code/sensitivity/laplacian/00_config.R code/laplacian_config.R
git mv code/sensitivity/laplacian/01_fit_basis_scores.R code/laplacian_fit_basis_scores.R
git mv code/sensitivity/laplacian/02_reconstruct_effect.R code/laplacian_reconstruct_effect.R
git mv code/sensitivity/laplacian/03_diagnose_hard.R code/laplacian_diagnose_hard.R
git mv code/sensitivity/laplacian/04_diagnose_continuous.R code/laplacian_diagnose_continuous.R
```

`code/utils.R` is already at its final path and is not moved.

- [ ] **Step 3: Archive the migration script and Laplacian README**

```bash
git mv code/truth/migrate_ground_truth.R code/archive/migrate_ground_truth.R
git mv code/sensitivity/laplacian/README.md code/archive/laplacian_README.md
```

- [ ] **Step 4: Confirm all category directories are empty and remove them**

Run:

```bash
find code/pipeline code/pilot code/diagnostics code/reporting \
     code/sensitivity code/truth code/config -type f -print
```

Expected: no output. Remove only those seven now-empty directories.

- [ ] **Step 5: Run the inventory test**

Run: `Rscript --vanilla tests/test_code_inventory.R`

Expected: PASS.

- [ ] **Step 6: Commit the moves**

```bash
git add code tests/test_code_inventory.R
git commit -m "refactor: restore flat numbered code layout"
```

### Task 3: Update R Dependency Paths

**Files:**
- Modify: `code/00_check_environment.R`
- Modify: `code/01_simulate_eeg_data.R`
- Modify: `code/02_generate_montecarlo.R`
- Modify: `code/05_aggregate_results.R`
- Modify: `code/06_generate_paired_radius_pilot.R`
- Modify: `code/08_fit_paired_radius_bgam.R`
- Modify: `code/09_diagnose_paired_radius_bgam.R`
- Modify: `code/12_validate_heterogeneity_bgam.R`
- Modify: `code/17_plot_2sd_bgam_diagnostics.R`
- Modify: `code/25_laplacian_sensitivity.R`
- Modify: `code/27_significance_alpha_sensitivity.R`
- Modify: `code/utils.R`
- Modify: path-sensitive tests under `tests/`

**Interfaces:**
- Consumes: final paths created in Task 2.
- Produces: scripts whose `source()` and resource checks resolve only final paths.

- [ ] **Step 1: Update main-pipeline source paths**

Apply these exact replacements in active files and tests:

```text
code/config/paper_config.R                  -> code/paper_config.R
code/truth/ground_truth.R                   -> code/ground_truth.R
code/pipeline/simulation.R                  -> code/01_simulate_eeg_data.R
code/pipeline/generate_montecarlo.R         -> code/02_generate_montecarlo.R
code/pipeline/fit_bgam.R                    -> code/03_fit_bgam.R
code/pipeline/compare_methods.R             -> code/04_compare_methods.R
code/pipeline/aggregate_results.R           -> code/05_aggregate_results.R
code/pipeline/fit_paired_radius.R           -> code/08_fit_paired_radius_bgam.R
code/diagnostics/35_diagnose_paired_radius_bgam_pilot.R
                                             -> code/09_diagnose_paired_radius_bgam.R
```

- [ ] **Step 2: Update Laplacian helper paths**

Replace the five `code/sensitivity/laplacian/*.R` source paths in
`code/25_laplacian_sensitivity.R` with the five final unnumbered helper paths.

- [ ] **Step 3: Update self-references in comments and usage errors**

Update file header names, `Rscript code/...` examples, and path-bearing error
messages to the final filename of each active script. Do not rewrite analysis
descriptions or statistical comments.

- [ ] **Step 4: Prove that active code contains no old category paths**

Run:

```bash
rg -n 'code/(pipeline|pilot|diagnostics|reporting|sensitivity|truth|config)/' \
  code tests --glob '!code/archive/**'
```

Expected: no output.

- [ ] **Step 5: Parse every active R file**

Run:

```bash
Rscript --vanilla -e 'files <- c(list.files("code", "[.]R$", full.names=TRUE), list.files("tests", "[.]R$", full.names=TRUE)); bad <- files[!vapply(files, function(x) {tryCatch({parse(x); TRUE}, error=function(e) FALSE)}, logical(1))]; stopifnot(length(bad) == 0L)'
```

Expected: exit status 0.

- [ ] **Step 6: Commit dependency-path updates**

```bash
git add code tests
git commit -m "fix: update dependencies for numbered code paths"
```

### Task 4: Update Slurm Entrypoints and Tests

**Files:**
- Modify: `slurm/formal/aggregate_paired_radius.sh`
- Modify: `slurm/formal/fit_paired_radius_dortmund.sh`
- Modify: `slurm/formal/fit_paired_radius_hpc2n.sh`
- Modify: `slurm/formal/generate_paired_radius.sh`
- Modify: `slurm/pilot/run_paired_radius_pilot.sh`
- Modify: every test containing an old `code/...` path literal.

**Interfaces:**
- Consumes: final executable filenames from Task 2.
- Produces: cluster jobs and tests that invoke the flat entrypoints.

- [ ] **Step 1: Update formal Slurm commands**

Use these commands:

```text
Rscript code/02_generate_montecarlo.R "$N_SIM"
Rscript code/08_fit_paired_radius_bgam.R "$RADIUS_RULE" "$REPLICATE_SEED"
Rscript code/10_aggregate_paired_radius_bgam.R
```

- [ ] **Step 2: Update pilot Slurm commands**

Use these commands:

```text
Rscript code/08_fit_paired_radius_bgam.R "$RADIUS_RULE" 1
Rscript code/09_diagnose_paired_radius_bgam.R "$RADIUS_RULE"
```

- [ ] **Step 3: Update test source paths**

Replace old category paths using the same exact path map from Task 3. Keep test
assertions and fixtures unchanged except for `tests/test_code_inventory.R`.

- [ ] **Step 4: Validate shell syntax**

Run:

```bash
for file in slurm/formal/*.sh slurm/pilot/*.sh; do bash -n "$file"; done
```

Expected: exit status 0.

- [ ] **Step 5: Run path- and orchestration-focused tests**

Run:

```bash
Rscript --vanilla tests/test_code_inventory.R
Rscript --vanilla tests/test_paired_radius_orchestration.R
Rscript --vanilla tests/test_paired_radius_validation.R
Rscript --vanilla tests/test_paired_radius_bgam_pilot.R
Rscript --vanilla tests/test_paired_radius_topomap.R
Rscript --vanilla tests/test_paired_radius_superpopulation_truth.R
Rscript --vanilla tests/test_2sd_bgam_diagnostic_figures.R
```

Expected: all commands exit successfully without fitting a model.

- [ ] **Step 6: Commit Slurm and test path updates**

```bash
git add slurm tests
git commit -m "fix: point jobs and tests to numbered scripts"
```

### Task 5: Rewrite the Reading Guide Around Execution Order

**Files:**
- Modify: `README.md`
- Modify: `PROJECT_GUIDE.md`
- Modify: `notes/diagnostic_metric_audit_2026-07-01.md`
- Modify: `notes/low_rank_laplacian_vs_full_bgam_summary.md`

**Interfaces:**
- Consumes: the final code inventory and Slurm entrypoints.
- Produces: documentation that lets a reader find the formal workflow by number.

- [ ] **Step 1: Replace the category-directory explanation**

Document the six number ranges exactly:

```text
00-05 main Monte Carlo pipeline
06-10 paired-radius pilot
11-15 BGAM and SESOI diagnostics
16-21 thesis and slide figures
22-27 prior, basis, Laplacian, and threshold sensitivity
unnumbered sourced helper modules
archive obsolete and one-time scripts
```

- [ ] **Step 2: Update every documented command and file link**

Apply the target file map to all `code/...` references in the four files. Keep
the existing scientific description and parameter values unchanged.

- [ ] **Step 3: Confirm no documentation points to removed directories**

Run:

```bash
rg -n 'code/(pipeline|pilot|diagnostics|reporting|sensitivity|truth|config)/' \
  README.md PROJECT_GUIDE.md notes
```

Expected: no output, except historical text explicitly marked as archived;
prefer updating those references to `code/archive/` as well.

- [ ] **Step 4: Commit documentation updates**

```bash
git add README.md PROJECT_GUIDE.md notes
git commit -m "docs: document numbered execution order"
```

### Task 6: Full Structural Regression and Statistical Diff Audit

**Files:**
- Verify: all changed files.

**Interfaces:**
- Consumes: completed Tasks 1-5.
- Produces: evidence that reorganisation is complete and behavior-preserving.

- [ ] **Step 1: Run every non-archived test script**

Run:

```bash
for file in tests/test_*.R; do Rscript --vanilla "$file"; done
```

Expected: all tests exit successfully. If a test attempts a model fit or writes
results, stop it and report it rather than expanding the authorized scope.

- [ ] **Step 2: Re-run R parsing and Slurm syntax checks**

Run the R parse command from Task 3 and the `bash -n` loop from Task 4.

Expected: both exit successfully.

- [ ] **Step 3: Verify continuous numbering and removed directories**

Run:

```bash
Rscript --vanilla tests/test_code_inventory.R
find code -mindepth 1 -maxdepth 1 -type d -print
```

Expected: inventory test passes and the only first-level code directory is
`code/archive`.

- [ ] **Step 4: Audit the diff for unintended statistical changes**

Run:

```bash
git diff origin/main -- code slurm tests README.md PROJECT_GUIDE.md notes
git diff --check origin/main
rg -n 'ROPE_UV|ODDS_THRESHOLD|adapt_delta|max_treedepth|k_xy|k_t|cluster|TFCE' code
```

Review every non-path content change. Revert any change to constants, formulas,
thresholds, priors, model controls, or metric definitions.

- [ ] **Step 5: Commit verification-only corrections if needed**

If verification requires path-only corrections:

```bash
git add code slurm tests README.md PROJECT_GUIDE.md notes
git commit -m "fix: complete numbered path migration"
```

If no corrections are needed, do not create an empty commit.

- [ ] **Step 6: Report the final branch and evidence**

Report the branch name, final `git status -sb`, numbered file list, tests run,
and any tests intentionally skipped. Do not push or open a pull request unless
the user explicitly requests publication.
