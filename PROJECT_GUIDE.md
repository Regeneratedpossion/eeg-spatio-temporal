# Project Guide

This file is the compact guide to the current hard-threshold + SEM BGAM
pipeline. It replaces the older inventory, reading guide, structure note, and
data index files.

## Main Question

Can a spatio-temporal BGAM recover a known simulated EEG effect over the full
channel x time grid?

The current simulation injects a condition effect around electrode `CZ` during
the true signal window `160-342 ms`.

## Read The Code In This Order

1. `code/utils.R`
   Shared paths, thresholds, resources, trial aggregation, topomap helpers.

2. `code/00_generate_data.R`
   Defines one simulated dataset and the hard-threshold spatial ground truth.

3. `code/00_generate_montecarlo.R`
   Generates 30 Monte Carlo datasets with `P = 20`, `R = 20`.

4. `code/01_fit_bgam.R`
   Fits the current BGAM:

   ```r
   eeg_mean | se(eeg_sem) ~
       condition +
       (1 + condition | participant) +
       t2(xproj, yproj, time,
          by = condition,
          bs = c("tp", "tp"),
          k = c(k_xy, k_t),
          d = c(2, 1),
          full = full)
   ```

5. `code/09_run_montecarlo_sem.R`
   Runs one Monte Carlo simulation: load data, fit BGAM, save posterior summary,
   compare methods, write one CSV.

6. `code/03_compare_methods.R`
   Builds time x channel masks for raw p, BH, Holm, cluster, TFCE, and BGAM.
   The largest connected blob gives onset and offset for one simulation.

7. `code/04_aggregate_results.R`
   Aggregates all simulations into Bias, MAE, RMSE, SD, F1_ST, and FP/FN
   summaries.

8. `code/05_slide_simulated_eeg_n20.R` to `code/08_slide_topomap_n20.R`
   Generate slide figures for simulated data, model fit, posterior probability,
   and topomaps.

9. `code/10_prior_predictive_check.R`
   Prior predictive check using the current SEM BGAM formula.

## Current Data Story

The active Monte Carlo data are generated, not stored in GitHub.

```text
data/montecarlo_n20_hard/
├── sim_001.rds ... sim_030.rds
└── gt_001.rds  ... gt_030.rds
```

Main settings:

- participants: `P = 20`;
- trials per participant per condition: `R = 20`;
- electrodes: 61;
- time grid: downsampled to 100 Hz;
- true signal window: `0.160-0.342 s`;
- seed electrode: `CZ`;
- spatial truth: Gaussian around `CZ`, set to zero outside one `sigma_space`.

Ground truth active electrodes are:

```text
spatial_distance <= sigma_space
```

## Current Result Folders

Generated results are local or cluster outputs and are ignored by Git.

```text
results/fits/         # fitted brms models
results/montecarlo/   # BGAM posterior summaries
results/comparisons/  # per-simulation CSVs and final tables
results/figures/      # slide and diagnostic figures
cmdstan_cache/        # compiled Stan code
```

## Current Slurm Files

Main scripts:

```text
slurm/gen_montecarlo_n20_hard.sh
slurm/run_montecarlo_n20_hard_sem.sh
slurm/run_montecarlo_n20_hard_sem_k1212_test.sh
```

The main run uses:

```text
MC_DATA_SUBDIR=montecarlo_n20_hard
MC_COMPARE_SUBDIR=montecarlo_n20_hard_sem
MC_SUMMARY_SUBDIR=montecarlo_n20_hard_sem
K_XY=20
K_T=20
```

The smaller basis test uses:

```text
K_XY=12
K_T=12
```

## Files Kept For Support

- `code/02_diagnostics.R`: detailed single-fit diagnostics.
- `code/check_env.R`: checks R/Python packages and project resources.
- `code/compute_trial_effective_sample_size.R`: trial-correlation check.
- `code/plot_pilot_preview.R`: quick pilot preview.
- `code/run_montecarlo.R`: legacy runner kept only for old scripts.

## Files Not Uploaded

The GitHub repository intentionally excludes:

- generated data;
- fitted models;
- figures;
- logs;
- CmdStan cache;
- old archive folders;
- large reference PDFs.
