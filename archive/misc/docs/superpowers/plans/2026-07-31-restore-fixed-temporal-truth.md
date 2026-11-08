# Restore Fixed Temporal Ground Truth Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Restore fixed 160/342 ms temporal truth without changing generated EEG or continuous effect fields.

**Architecture:** Keep temporal truth construction inside `build_ground_truth()` and add one narrow migration script for existing paired-radius RDS files. Protect the behavior with an integration-style R regression test using the real generator functions.

**Tech Stack:** R 4.x, base R, dplyr, RDS files.

## Global Constraints

- Do not modify simulated EEG data, model code, spatial truth, or fitted models.
- The 100 Hz active grid is exactly `seq(0.16, 0.34, by = 0.01)`.
- Preserve every `continuous_effect` object exactly.

---

### Task 1: Lock the fixed temporal truth contract

**Files:**
- Modify: `tests/test_heterogeneity_continuous_truth.R`
- Modify: `code/00_generate_data.R`

**Interfaces:**
- Consumes: `build_ground_truth(params, sensors, participant_parameters)`.
- Produces: fixed `true_onset`, `true_offset`, `temporal_active_times`, and `active_cells`.

- [ ] Add literal assertions for onset 0.160, offset 0.342, and active grid 0.16:0.34.
- [ ] Run `Rscript --vanilla tests/test_heterogeneity_continuous_truth.R` and verify failure because the current branch returns 0.20/0.30.
- [ ] Replace threshold-derived temporal labels with parameter-derived labels in `build_ground_truth()`.
- [ ] Re-run the test and verify it passes.

### Task 2: Migrate existing paired-radius truth labels

**Files:**
- Create: `code/38_repair_paired_radius_temporal_truth.R`
- Modify in place: `data/montecarlo_r1sd/gt_001.rds`
- Modify in place: `data/montecarlo_rhalfmax/gt_001.rds`
- Modify in place: `data/montecarlo_r2sd/gt_001.rds`

**Interfaces:**
- Consumes: each current ground-truth list.
- Produces: the same list with corrected temporal label fields only.

- [ ] Record hashes for each `sim_001.rds` and serialized `continuous_effect`.
- [ ] Implement a migration that validates required fields, sets fixed boundaries and grid, rebuilds `active_cells`, removes the obsolete `temporal_threshold`, and writes the RDS.
- [ ] Run the migration once for the three pilot directories.
- [ ] Confirm all continuous-effect and simulated-data hashes are unchanged.

### Task 3: Full verification

**Files:**
- Test: `tests/test_heterogeneity_continuous_truth.R`
- Test: `tests/test_paired_hard_radius_generator.R`
- Test: `tests/test_paired_radius_validation.R`

**Interfaces:**
- Consumes: current source and migrated RDS files.
- Produces: verification evidence only.

- [ ] Parse the modified R scripts.
- [ ] Run the three relevant regression tests.
- [ ] Run a direct three-file audit of onset, offset, grid, active-cell counts, spatial zeros, and continuous-effect preservation.
- [ ] Inspect `git diff` to ensure no unrelated file was changed by this operation.

