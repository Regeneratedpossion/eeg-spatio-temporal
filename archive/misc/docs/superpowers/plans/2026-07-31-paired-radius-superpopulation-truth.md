# Paired-Radius Superpopulation Truth Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Recompute the paired-radius diagnostic truth as the data-generating superpopulation mean and rerun only the fitted 2 SD model.

**Architecture:** Add a small pure R helper that reconstructs the expected ERP waveform from the original template, the prespecified discrete onset distribution, and `E(A) = bump_amplitude`. The existing diagnostic consumes the resulting complete spatial-temporal truth surface while retaining the realised seed-1 truth only as an audit field.

**Tech Stack:** R 4.4, base R, dplyr, brms/posterior, MNE through reticulate, test scripts using `stopifnot`.

## Global Constraints

- The estimand is the data-generating superpopulation mean effect.
- Use the exact onset support and ERP template used by the generator.
- Do not use realised participant amplitudes or onset draws in population truth.
- Keep active labels, SESOI, posterior odds, cluster-mass, TFCE, FPR, leakage, EDF, convergence checks, and the residual topomap unchanged.
- Rerun only the 2 SD fitted model; do not refit BGAM.
- For this rerun, compute BGAM diagnostics only; do not execute cluster-mass or TFCE permutations.

---

### Task 1: Superpopulation Truth Helper

**Files:**
- Create: `code/paired_radius_superpopulation_truth.R`
- Create: `tests/test_paired_radius_superpopulation_truth.R`

**Interfaces:**
- Consumes: a saved paired-radius `ground_truth` list and the path to `resources/brms_meeg/code/erp_template.R`.
- Produces: `paired_radius_superpopulation_truth(ground_truth, template_path, channels = NULL, times = NULL) -> data.frame(channel, time, effect)`.

- [ ] **Step 1: Write the failing test**

Create a regression test that sources the not-yet-existing helper, loads `data/montecarlo_r2sd/gt_001.rds`, and asserts:

```r
population <- paired_radius_superpopulation_truth(gt, template_path)
cz <- subset(population, channel == "CZ" & abs(time - 0.25) < 1e-12)$effect
stopifnot(isTRUE(all.equal(cz, 1.457449, tolerance = 1e-6)))
stopifnot(all(subset(population, !channel %in% gt$active_electrodes)$effect == 0))
```

Clone the object, replace participant draws and realised continuous effects,
and assert the returned population surface is identical. Remove
`onset_support_ms` and assert a clear validation error.

- [ ] **Step 2: Run test to verify it fails**

Run:

```bash
PROJECT_ROOT="$PWD" Rscript --vanilla tests/test_paired_radius_superpopulation_truth.R
```

Expected: failure because `code/paired_radius_superpopulation_truth.R` does not exist.

- [ ] **Step 3: Implement the minimal helper**

Load the ERP template into an isolated environment, reconstruct the generator's
onset-shifted unit waveforms, average over `ground_truth$onset_support_ms`,
multiply by `ground_truth$bump_amplitude`, downsample with
`ground_truth$downsample_factor`, and multiply by each channel's unique saved
hard-truncated `spatial_weight`. Validate required fields, unique channel
weights, requested channels, and requested times.

- [ ] **Step 4: Run test to verify it passes**

Run the command from Step 2. Expected: exit 0 and
`paired-radius superpopulation truth tests passed`.

### Task 2: Connect the Diagnostic to the New Estimand

**Files:**
- Modify: `code/35_diagnose_paired_radius_bgam_pilot.R`
- Modify: `tests/test_paired_radius_superpopulation_truth.R`

**Interfaces:**
- Consumes: `paired_radius_superpopulation_truth()` from Task 1.
- Produces: diagnostic summaries with `truth_estimand`, population
  `truth_cz_250`, audit-only `realized_truth_cz_250`, population `peak_ratio`,
  and population projection `recovery`.

- [ ] **Step 1: Add failing static integration assertions**

Assert that the diagnostic sources the helper, calls
`paired_radius_superpopulation_truth`, writes
`truth_estimand = "superpopulation"`, and reports
`realized_truth_cz_250` separately.

- [ ] **Step 2: Run test to verify it fails**

Run the Task 1 test command. Expected: failure because the diagnostic has not
yet been connected to the helper.

- [ ] **Step 3: Make the minimal diagnostic changes**

Source the helper. Construct both `truth_at_peak` from the superpopulation
helper and `realized_truth_at_peak` from `continuous_truth_cells`. Use only the
population surface for peak ratio, projection recovery, and FP truth values.
Add the estimand label and realised audit value to `diagnostic_summary.csv`.

- [ ] **Step 4: Run focused and existing regression tests**

Run:

```bash
PROJECT_ROOT="$PWD" Rscript --vanilla tests/test_paired_radius_superpopulation_truth.R
PROJECT_ROOT="$PWD" Rscript --vanilla tests/test_paired_radius_bgam_pilot.R
PROJECT_ROOT="$PWD" Rscript --vanilla tests/test_paired_radius_topomap.R
```

Expected: all three commands exit 0.

### Task 3: Rerun and Verify the 2 SD Diagnostic

**Files:**
- Produce: `results/paired_radius_bgam_pilot/2sd/diagnostic_summary.csv`
- Produce: `results/paired_radius_bgam_pilot/2sd/method_metrics_250ms.csv`
- Produce: `results/paired_radius_bgam_pilot/2sd/fp_details_250ms.csv`
- Produce: `results/paired_radius_bgam_pilot/2sd/edf_summary.csv`
- Produce: `results/paired_radius_bgam_pilot/2sd/edf_components.csv`
- Produce: `results/paired_radius_bgam_pilot/2sd/posterior_map_seed001.rds`
- Produce: `results/paired_radius_bgam_pilot/2sd/empirical_minus_bgam_250ms.csv`
- Produce: `results/paired_radius_bgam_pilot/2sd/empirical_minus_bgam_250ms.png`
- Produce: `logs/diagnose_paired_radius_2sd_superpopulation_20260731.log`

**Interfaces:**
- Consumes: the existing 2 SD `fit_seed001.rds`, its model config, and the unchanged simulation dataset.
- Produces: one complete 2 SD diagnostic result set under the new estimand.

- [ ] **Step 1: Run only the 2 SD diagnostic**

```bash
LC_ALL=en_US.UTF-8 LANG=en_US.UTF-8 PYTHONIOENCODING=UTF-8 \
PROJECT_ROOT="$PWD" Rscript --vanilla \
  code/35_diagnose_paired_radius_bgam_pilot.R 2sd bgam-only \
  2>&1 | tee logs/diagnose_paired_radius_2sd_bgam_only_superpopulation_20260731.log
```

Expected: exit 0 with one printed 2 SD summary and a one-row BGAM method table,
with no cluster-mass or TFCE permutation output.

- [ ] **Step 2: Verify outputs and invariant metrics**

Read all expected files and assert:

```r
summary$truth_estimand == "superpopulation"
all.equal(summary$truth_cz_250, 1.457449, tolerance = 1e-6)
all.equal(summary$realized_truth_cz_250, 1.118975, tolerance = 1e-6)
summary$n_active == 21
summary$n_inactive == 40
```

Compare the new output with
`results/paired_radius_bgam_pilot/2sd/pre_rerun_20260731_dortmund_code` and
assert that BGAM peak, FP/FPR, leakage, convergence statistics, EDF values,
method metrics, residual values, and the posterior map are unchanged. Only the
estimand label, truth values, peak ratio, and recovery may differ.

- [ ] **Step 3: Inspect the final log**

Confirm the log contains no `Error`, `UnicodeEncodeError`, or interrupted
permutation run, and report any nonfatal package-version warnings separately.
