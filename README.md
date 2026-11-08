# Bayesian GAM for spatio-temporal EEG inference

Code of the master thesis. It simulates EEG data with a known effect, fits a Bayesian generalised additive model (BGAM) with brms and cmdstanr, and compares it with cluster-based permutation tests and TFCE.

## What the thesis does in one paragraph

Four scenarios (A = amplitude, ω_A = between-participant variability of the amplitude) are simulated on a 61-channel layout:

| Scenario | Name in code | A | ω_A |
|---|---|---:|---:|
| 1 | `cond_I` | 0.5 | 0.26 |
| 2 | `cond_H` | 0.5 | 0.70 |
| 3 | `cond_F` | 1.5 | 0.26 |
| 4 | `cond_E` | 1.5 | 0.70 |

Each scenario has 100 simulated datasets (seeds 1–100). The BGAM is fitted to every dataset, and 362 of the 400 fits pass the convergence criterion (`resources/passed_fits.txt`). Onset, offset and spatial extent of the effect are then compared with the conventional methods.

## Structure

```text
eeg-spatio-temporal/
├── code/
│   ├── 00_config.R          paths and scenario map, loaded by every script
│   ├── 01_simulate/         simulate the datasets
│   ├── 02_fit/              fit the Bayesian model (one seed per job)
│   ├── 03_conventional/     cluster test, TFCE and other conventional methods
│   ├── 04_evaluate/         turn each fit into detection masks and measures
│   ├── 05_aggregate/        collect the results of all seeds
│   ├── 06_tables/           tables of the thesis
│   ├── 07_figures/          figures of the thesis and their data
│   └── R/                   helper functions
├── slurm/                   cluster job scripts, same numbers as code/
├── resources/               small input files (electrodes, adjacency, passed fits, brms template)
├── inventory.csv            which script makes which table or figure
└── README.md
```

Old code, notes and results from earlier versions of the project are not on this branch. They are on the branch `archive`.

## Where the parameters are

| Parameters | File and function |
|---|---|
| Simulation (participants, trials, noise, onsets, the four scenarios) | `code/01_simulate/01_generate_data.R`: `define_simulation_parameters()`, `heterogeneity_scenarios()` |
| BGAM (basis, k, chains, iterations, priors, contrast coding) | `code/R/fit_config.R`: `ds_cr_seed_config()` |
| Evaluation (SESOI baseline, posterior odds, active window, thresholds) | `code/R/evaluation_functions.R`: `evaluation_settings()` |
| Conventional methods (permutations, seed, alpha, TFCE) | `code/03_conventional/03_evaluate_conventional.R`: `conventional_settings()` |
| Convergence rule (R-hat < 1.01, ESS > 400) | `code/05_aggregate/05_make_passed_fits.R` |

## Order of the scripts

| Step | Script | Where it runs |
|---|---|---|
| 00 | `00_config.R` (paths and scenario map) | loaded by the other scripts |
| 01 | `01_simulate/01_generate_dataset.R` simulates one dataset (`COND` = I, H, F or E, seed 1–100) | cluster or local |
| 02 | `02_fit/02_fit_bgam.R` fits the BGAM of one seed (`COND` = I, H, F or E) | cluster, about 20 h per fit |
| 03 | `03_conventional/03_evaluate_conventional*.R` runs the conventional tests | cluster |
| 04 | `04_evaluate/04_evaluate_bgam.R` turns each fit into detection masks and measures | cluster |
| 05 | `05_aggregate/05_*.R` collects the per-seed results | cluster |
| 06 | `06_tables/06_tables.R` makes the result tables | local |
| 07 | `07_figures/07_figure_*.R` makes the figures | local |

Run every script from the repository root. Details, environment variables and the Slurm call are in [code/README.md](code/README.md).

## To reproduce only the tables and figures

The small per-seed results that steps 06 and 07 need are in `results/` and `data/*/gt_*.rds`, so a clone can make the tables and figures without fitting anything: `Rscript code/06_tables/06_tables.R passed`, then the scripts in `code/07_figures/`. Fits, simulated EEG data and logs are not stored.

## Software

R 4.4.2, brms 2.23.0, cmdstanr 0.8.0 with CmdStan 2.38.0, mgcv 1.9-4, posterior 1.7.0. The conventional tests use MNE-Python 0.23.4 through reticulate.

## Known issues

1. **Folder for tables and figures.** `THESIS_DIR` defaults to `results/` in the repository (see `code/00_config.R`). Tables go to `results/fig_data/tables/`, figure data to `results/fig_data/`, figures to `results/diagnostics/` and `results/figures_out/`.
2. **Summaries between steps 05 and 06.** `06_tables.R` reads `sesoi_rule_passed.csv` and `sesoi_rule_20261009.csv` from `results/fig_data/tables/`. They are made by `05_agg_sesoi_passed.R` and `05_agg_sesoi_all.R`, which now write to that folder (`AGG_OUT` overrides it). `05_agg_ring_fpr.R` and `05_agg_threshold_sweep.R` write to `results/fig_data/`. Run all four from the project root after the per-seed evaluations exist; they reproduce the stored files up to rounding in the last digit.
