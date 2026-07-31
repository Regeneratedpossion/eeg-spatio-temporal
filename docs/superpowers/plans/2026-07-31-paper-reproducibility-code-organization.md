# Paper Reproducibility and Compatible Code Organisation Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox syntax for tracking.

**Goal:** Make the paired-radius high-heterogeneity analysis the reproducible thesis pipeline, correct truth and documentation bugs, and organise supporting scripts without breaking old entry points.

**Architecture:** A single configuration module defines the paper contract. Explicit ground-truth interfaces separate theoretical population truth, finite-sample realised truth, and timing labels. Formal implementations move behind compatibility entry points; supporting analyses move by purpose while old paths remain executable adapters.

**Tech Stack:** R 4.x, dplyr, tidyr, brms/cmdstanr, base-R RDS migration, Slurm, LaTeX.

## Global Constraints

- Formal simulation: paired 1sd, half-max, 2sd; scenario high; P=20; R=20.
- Formal BGAM: P-splines, k=c(16,12), full=FALSE, four chains, 3,000 iterations, 1,000 warmup.
- Preserve model formula, priors, adapt_delta, max_treedepth, SESOI, posterior odds, cluster-mass, TFCE, and detection rules.
- Spatial signal is exactly zero outside each radius.
- Temporal waveforms are finite and participant-shifted; do not impose a second common 160--342 ms clip.
- Reference timing estimands are 160 and 342 ms; the 100 Hz scoring grid is 0.16 through 0.34 s.
- Do not modify simulated EEG observations or fitted BGAM objects.
- Do not touch running Umea or Dortmund jobs.
- Preserve unrelated user changes in the dirty worktree.
- Do not push without explicit approval.

---

### Task 1: Lock the formal paper configuration

**Files:**
- Create: code/config/paper_config.R
- Create: tests/test_paper_config_contract.R
- Modify: code/00_generate_data.R
- Modify: code/34_fit_paired_radius_bgam_pilot.R

**Interfaces:**
- Produces paper_analysis_config(), returning simulation, radii, model, and evaluation sections.
- Uses base R only and has no source-time side effects.

- [ ] Write a failing test asserting scenario high, P=20, R=20, hard_signal, onset 0.160, offset 0.342, k_xy=16, k_t=12, P-splines, sigma_true=TRUE, and active counts 5/9/21.
- [ ] Run PROJECT_ROOT="$PWD" Rscript --vanilla tests/test_paper_config_contract.R and confirm failure because the module is missing.
- [ ] Implement paper_analysis_config() with the exact values from the approved design.
- [ ] Make radius_specifications() and paired_radius_bgam_config() read this contract while preserving their existing return shapes.
- [ ] Run test_paper_config_contract.R, test_paired_radius_bgam_pilot.R, and test_paired_hard_radius_generator.R; require all to pass.
- [ ] Commit only these files with message: refactor: centralize formal paper configuration.

The test must include literal assertions equivalent to:

    cfg <- paper_analysis_config()
    stopifnot(
        identical(cfg$simulation$scenario, "high"),
        cfg$simulation$n_ppt == 20L,
        cfg$simulation$n_trials == 20L,
        identical(cfg$model$k_xy, 16L),
        identical(cfg$model$k_t, 12L),
        identical(cfg$radii$expected_active, c(5L, 9L, 21L))
    )

---

### Task 2: Make theoretical and realised truth explicit

**Files:**
- Create: code/truth/ground_truth.R
- Create: tests/test_explicit_truth_estimands.R
- Modify: code/00_generate_data.R
- Modify: code/paired_radius_superpopulation_truth.R
- Modify: tests/test_paired_radius_superpopulation_truth.R
- Modify: tests/test_heterogeneity_continuous_truth.R

**Interfaces:**
- Produces theoretical_population_effect(), realised_effect_field(), truth_effect_cells(), and add_explicit_truth_fields().
- paired_radius_superpopulation_truth() remains as a compatibility adapter.

- [ ] Write a failing fixture test requiring theoretical_effect, realised_effect, realised_mean_waveform, and compatibility aliases.
- [ ] Change participant onset and amplitude draws in the fixture; assert theoretical truth is unchanged and realised truth changes.
- [ ] Run the test and confirm failure on missing explicit fields.
- [ ] Move the tested superpopulation calculation into code/truth/ground_truth.R.
- [ ] Make build_ground_truth() store both explicit fields.
- [ ] Keep continuous_effect byte-identical to realised_effect.
- [ ] Keep population_mean_waveform byte-identical to realised_mean_waveform during compatibility.
- [ ] Convert code/paired_radius_superpopulation_truth.R into an adapter to theoretical_population_effect().
- [ ] Run explicit-truth, superpopulation-truth, heterogeneity-truth, and noninterference tests.
- [ ] Commit with message: fix: expose theoretical and realised ground truth.

The formal selector has this interface:

    truth_effect_cells <- function(
        ground_truth,
        estimand = c("theoretical", "realised"),
        template_path = NULL,
        channels = NULL,
        times = NULL
    )

Formal callers must pass estimand explicitly.

---

### Task 3: Lock participant-specific temporal support

**Files:**
- Create: tests/test_temporal_support_semantics.R
- Modify: code/00_generate_data.R

**Interfaces:**
- Consumes temporal_signal_at_onset(params, onset_ms).
- Produces exact shifted finite support and separate fixed reference labels.

- [ ] Test onset values 140, 160, and 180 ms.
- [ ] Assert each participant waveform is exactly zero outside its own shifted finite support.
- [ ] Assert the 140 ms waveform has positive values before 160 ms.
- [ ] Assert the 180 ms waveform has positive values after 342 ms.
- [ ] Assert the theoretical population field has small non-zero tails outside the reference interval.
- [ ] Run the characterisation test against current code; it must pass before any edit.
- [ ] Add comments distinguishing participant signal support from reference timing labels; do not add temporal clipping.
- [ ] Re-run and commit with message: test: lock participant-specific temporal support.

---

### Task 4: Safely migrate all 15 paired-radius truth files

**Files:**
- Create: code/truth/migrate_ground_truth.R
- Modify: code/38_repair_paired_radius_temporal_truth.R
- Modify: tests/test_repair_paired_radius_temporal_truth.R
- Modify locally: gt_001.rds through gt_005.rds in each of the three paired-radius directories

**Interfaces:**
- Produces migrate_ground_truth_file(path, backup_dir, template_path).
- The old repair script delegates to the new migration module.

- [ ] Extend the migration test before implementation.
- [ ] Require explicit truth fields, fixed timing labels, no temporal_threshold, and unchanged continuous_effect.
- [ ] Run the test and observe the expected missing-field failure.
- [ ] Implement discovery of exactly the three paired directories and all gt_*.rds files.
- [ ] Validate scenario high, radius rule, and expected active count before writing.
- [ ] Use a non-overwriting backup, temporary RDS, read-back validation, and atomic rename.
- [ ] Record SHA-256 for all 15 sim_*.rds and digests for all existing continuous_effect fields.
- [ ] Run the migration once for all 15 GT files.
- [ ] Compare every stable field with backups and confirm all simulated-data hashes are unchanged.
- [ ] Commit migration code and tests only, with message: fix: migrate paired truth to explicit estimands.

Each migrated file must satisfy:

    identical(gt$true_onset, 0.160)
    identical(gt$true_offset, 0.342)
    identical(gt$temporal_active_times, seq(0.16, 0.34, 0.01))
    identical(gt$continuous_effect, gt$realised_effect)

---

### Task 5: Require explicit truth in formal evaluation

**Files:**
- Modify: code/utils.R
- Modify: code/35_diagnose_paired_radius_bgam_pilot.R
- Modify: code/37_plot_2sd_bgam_diagnostics.R
- Create: tests/test_truth_metric_semantics.R
- Modify: tests/test_continuous_ground_truth.R
- Modify: tests/test_2sd_bgam_diagnostic_figures.R

**Interfaces:**
- Produces zero_truth_leakage(theta_hat, theta).
- continuous_truth_cells() remains a compatibility adapter returning realised truth.

- [ ] Write a failing test with small non-zero temporal tails.
- [ ] Demonstrate that label-based FP can exceed exact-zero FP.
- [ ] Test leakage as energy where the selected continuous truth equals zero.
- [ ] Implement explicit theoretical selection for formal paired diagnostics.
- [ ] Implement explicit realised selection for empirical-minus-injected-signal diagnostics.
- [ ] Preserve current spatial_fp_metrics() and spatial_leakage() at 250 ms.
- [ ] Add metric_scope="spatial_250ms" to paired summaries.
- [ ] Run truth-metric, continuous-truth, 2sd-figure, and paired-validation tests.
- [ ] Commit with message: fix: distinguish label and zero-truth metrics.

The leakage regression assertion is:

    expected <- sum(theta_hat[theta == 0]^2) / sum(theta_hat^2)
    stopifnot(zero_truth_leakage(theta_hat, theta) == expected)

---

### Task 6: Correct the thesis, README, and project guide

**Files:**
- Modify: /Users/jiahuifan/Desktop/thesis/thesis_script.tex
- Modify: README.md
- Modify: PROJECT_GUIDE.md
- Create: tests/test_thesis_contract.R

**Interfaces:**
- Consumes paper_analysis_config() and THESIS_TEX.
- Produces a local code--paper consistency check.

- [ ] Write a failing thesis-contract test that rejects the m_tau+s_tau temporal definition.
- [ ] Assert formal settings and the theoretical/realised distinction appear in the thesis source.
- [ ] Run with THESIS_TEX=/Users/jiahuifan/Desktop/thesis/thesis_script.tex and confirm failure.
- [ ] Replace the mean-plus-SD subsection with the fixed timing-scoring set.
- [ ] Explain that 160/342 ms are reference estimands and the 100 Hz grid ends at 340 ms.
- [ ] Explain participant-specific finite support and population-mean tails.
- [ ] Describe 250 ms FP/FPR as exact-zero spatial metrics.
- [ ] Describe full-grid outside-window detections as label-based metrics.
- [ ] Define continuous leakage by exact-zero cells of the selected continuous truth.
- [ ] Update README and PROJECT_GUIDE to paired radii, high, P=20, k=(16,12), P-splines, sigma=TRUE, and the new module map.
- [ ] Mark old P=5 and single-radius results as historical if retained.
- [ ] Run the thesis-contract test and compile with latexmk using halt-on-error.
- [ ] Commit repository docs with message: docs: align repository with formal paired analysis.
- [ ] If the thesis directory is a Git repository, commit thesis_script.tex and regenerated PDF separately; otherwise leave verified modifications and report this.

The replacement timing definition is:

    \mathcal{T}_{\mathrm{score}}
    =\{t\in\mathcal{G}:0.160\leq t\leq0.342\}.

---

### Task 7: Move formal implementations behind compatibility entry points

**Files:**
- Create: code/pipeline/simulation.R
- Create: code/pipeline/generate_montecarlo.R
- Create: code/pipeline/fit_bgam.R
- Create: code/pipeline/fit_paired_radius.R
- Create: code/pipeline/compare_methods.R
- Create: code/pipeline/aggregate_results.R
- Modify old entry points: 00_generate_data.R, 00_generate_montecarlo.R, 01_fit_bgam.R, 03_compare_methods.R, 04_aggregate_results.R, 34_fit_paired_radius_bgam_pilot.R
- Create: tests/test_compatibility_entrypoints.R

**Interfaces:**
- Formal modules define functions without running work when sourced.
- Old numbered files preserve public function names and direct-execution behaviour.

- [ ] Write a failing source-safety and function-formals parity test.
- [ ] Assert sourcing a formal module creates no data or result files.
- [ ] Move simulation and fitting function bodies first.
- [ ] Replace old files with thin source adapters.
- [ ] Expose main(args=commandArgs(trailingOnly=TRUE)) only for executable modules.
- [ ] Move comparison and aggregation implementations.
- [ ] Move MNE imports and directory creation out of source time and into explicit loaders or main().
- [ ] Generalise run_paired_radius_bgam_fit(radius_rule, seed=1L).
- [ ] Read sim_%03d.rds and write fit_seed%03d.rds using seed.
- [ ] Verify the old one-argument command remains seed 1.
- [ ] Run compatibility, orchestration, and paired-fit tests.
- [ ] Commit with message: refactor: deepen formal paired analysis pipeline.

---

### Task 8: Classify supporting scripts without breaking old paths

**Files:**
- Create directories: code/reporting, code/diagnostics, code/sensitivity, code/pilot, code/exploratory, code/archive
- Move supporting implementations
- Retain adapters at old documented or Slurm-referenced paths
- Create: tests/test_code_inventory.R

**Interfaces:**
- Supporting scripts remain executable from old and canonical paths.
- The inventory test is the authoritative classification list.

- [ ] Write a failing inventory test with this exact mapping.

Reporting:
- 05_slide_simulated_eeg_n20.R
- 06_slide_bgam_fit_n20.R
- 07_slide_posterior_probability_n20.R
- 08_slide_topomap_n20.R
- 07_plot_paired_radius_truth_pilot.R
- 31_plot_all_electrode_sesoi_posterior.R
- 37_plot_2sd_bgam_diagnostics.R

Diagnostics:
- 02_diagnostics.R
- 19_practical_null_fpr_rebound.R
- 28_delta_sesoi_probe.R
- 29_validate_heterogeneity_bgam.R
- 35_diagnose_paired_radius_bgam_pilot.R

Sensitivity:
- 10_prior_predictive_check.R
- 11_basis_screening.R
- 12_laplacian_basis_brms_fixed.R
- 18_basis_screening_full_ab_sigma_true_control.R
- 32_threshold_sensitivity.R
- 33_alpha_sensitivity.R
- laplacian directory

Pilot:
- 00_generate_paired_radius_pilot.R
- 00_validate_paired_radius_pilot.R
- 36_aggregate_paired_radius_bgam_pilot.R

Exploratory:
- compute_trial_effective_sample_size.R

Archive:
- 09_run_montecarlo_sem.R implementation

Keep utils.R and check_env.R at top level.

- [ ] Move one category at a time and add old-path adapters.
- [ ] Do not alter constants or output paths during mechanical moves.
- [ ] Update internal source references to canonical paths.
- [ ] Parse every R file recursively.
- [ ] Run bash -n on every Slurm script.
- [ ] Run inventory, topomap, threshold-sensitivity, and alpha-sensitivity tests.
- [ ] Commit each category separately using message: refactor: classify CATEGORY scripts.

---

### Task 9: Full verification and handoff

**Files:**
- Verify all modified code, tests, documentation, thesis, and local GT files.
- Do not create or modify fitted models.

**Interfaces:**
- Produces verification evidence and a final code--paper comparison table.

- [ ] Run every tests/test_*.R with PROJECT_ROOT and THESIS_TEX set.
- [ ] Generate one high/P20 paired replicate in a temporary directory.
- [ ] Assert active counts 5/9/21, exact spatial zeros, identical paired residuals, explicit truth fields, and fixed timing labels.
- [ ] Re-audit all 15 migrated GT files against backups.
- [ ] Compile the thesis and inspect the Estimands and ground truth pages.
- [ ] Inspect Git status and scope-specific diffs.
- [ ] Confirm no fitted model, results file, SSH file, or unrelated user file was staged.
- [ ] Report each formal setting, configuration location, pipeline consumer, thesis location, and passing test.
- [ ] Report explicitly that no BGAM was refitted and no cluster task was modified.

The full test loop is:

    task_exit_code=0
    for test_file in tests/test_*.R; do
      PROJECT_ROOT="$PWD" \
      THESIS_TEX=/Users/jiahuifan/Desktop/thesis/thesis_script.tex \
      Rscript --vanilla "$test_file" || task_exit_code=1
    done
    exit "$task_exit_code"

