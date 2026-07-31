# 2 SD BGAM Diagnostic Figures Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Generate the approved BGAM fit, SESOI posterior-probability, and multi-layer topomap figures for the existing 2 SD fit.

**Architecture:** Create one focused R plotting script that reads the existing fit, simulation data, posterior map, and superpopulation-truth helper. Separate data preparation functions from plotting functions so probability and landmark definitions can be unit tested without rerunning posterior prediction.

**Tech Stack:** R 4.4, brms, posterior, dplyr, tidyr, ggplot2, patchwork, base R PNG output.

## Global Constraints

- Use only the existing 2 SD BGAM fit; do not refit it.
- Do not run cluster-mass, TFCE, or other comparison methods.
- Use superpopulation truth from `paired_radius_superpopulation_truth()`.
- Define posterior probability as `P(cond2 - cond1 > 0.1 uV | data)`.
- Use odds greater than 10, equivalent to probability greater than `10 / 11`.
- Plot `CZ`, `FCZ`, `CPZ`, `C2`, `PZ`, and `FP1`.
- Use 200, 250, and 300 ms as diagnostic landmarks.

---

### Task 1: Testable Diagnostic-Figure Data Functions

**Files:**
- Create: `code/37_plot_2sd_bgam_diagnostics.R`
- Create: `tests/test_2sd_bgam_diagnostic_figures.R`

**Interfaces:**
- Produces: `diagnostic_channels(available_channels) -> character`.
- Produces: `diagnostic_landmarks(ground_truth, available_times) -> numeric`.
- Produces: `summarise_difference_draws(diff_draws, grid, sesoi, odds_threshold) -> data.frame` containing posterior mean, 95% CrI, `p_sig`, and `above_threshold`.

- [ ] **Step 1: Write the failing test**

Source the not-yet-existing plotting script. Assert the fixed channel order and
landmarks. Supply a small known posterior-difference matrix and assert that
`p_sig` counts draws strictly greater than `0.1`, the credible limits are draw
quantiles, and `above_threshold` uses `10 / 11`.

- [ ] **Step 2: Verify the test fails**

Run:

```bash
PROJECT_ROOT="$PWD" Rscript --vanilla tests/test_2sd_bgam_diagnostic_figures.R
```

Expected: failure because `code/37_plot_2sd_bgam_diagnostics.R` is absent.

- [ ] **Step 3: Implement the data functions**

Add the three interfaces above. Validate draw-column count against the grid,
require all six channels, and map requested landmarks to exact available
analysis times.

- [ ] **Step 4: Verify the test passes**

Run the command from Step 2. Expected: exit 0 and
`2sd BGAM diagnostic figure tests passed`.

### Task 2: Generate and Verify the Three Figures

**Files:**
- Modify: `code/37_plot_2sd_bgam_diagnostics.R`
- Modify: `tests/test_2sd_bgam_diagnostic_figures.R`
- Produce: `results/paired_radius_bgam_pilot/2sd/figures/bgam_fit_superpopulation.png`
- Produce: `results/paired_radius_bgam_pilot/2sd/figures/posterior_probability_sesoi.png`
- Produce: `results/paired_radius_bgam_pilot/2sd/figures/topomap_diagnostics_superpopulation.png`
- Produce: `results/paired_radius_bgam_pilot/2sd/figures/diagnostic_plot_data.rds`

**Interfaces:**
- Consumes: `fit_seed001.rds`, `sim_001.rds`, `gt_001.rds`, and `posterior_map_seed001.rds`.
- Produces: `run_2sd_bgam_diagnostic_figures()` and the four output files.

- [ ] **Step 1: Add failing integration assertions**

Assert statically that the script calls `posterior_epred` with
`re_formula = NA`, uses `ROPE_UV`, calls
`paired_radius_superpopulation_truth`, and writes exactly the approved output
filenames.

- [ ] **Step 2: Verify the integration assertions fail**

Run the Task 1 test command. Expected: failure because figure generation is
not implemented.

- [ ] **Step 3: Implement figure generation**

Construct six-channel posterior condition differences once and reuse the
draws for fit and probability figures. Construct the topomap layers from the
saved all-channel posterior mean, participant-aggregated empirical difference,
and reconstructed superpopulation truth. Use `idw_interp()` and
`scalp_decorations()` for the same scalp geometry as existing project plots.
Save all plotted numerical data in one RDS list.

- [ ] **Step 4: Run the plotting script**

```bash
LC_ALL=en_US.UTF-8 LANG=en_US.UTF-8 PROJECT_ROOT="$PWD" \
Rscript --vanilla code/37_plot_2sd_bgam_diagnostics.R
```

Expected: exit 0 and four non-empty output files.

- [ ] **Step 5: Verify files and data**

Check all PNG files are readable, have nonzero dimensions, and exceed 10 KB.
Read the RDS and assert six channels, three topomap times, SESOI `0.1`,
probability threshold `10 / 11`, superpopulation CZ truth at 250 ms near
`1.457449`, and no conventional-method objects.

- [ ] **Step 6: Run the full R test suite**

Execute every `tests/test_*.R` file with `PROJECT_ROOT="$PWD"`. Expected: all
commands exit 0; package-build-version warnings may remain nonfatal.
