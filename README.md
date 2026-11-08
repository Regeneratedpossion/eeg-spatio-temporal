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

## Folders

| Folder | Content |
|---|---|
| `code/` | all R and Python scripts, numbered in the order they are run |
| `code/R/` | helper functions that the numbered scripts load |
| `slurm/` | one cluster submit script per step, numbered like `code/` |
| `resources/` | small input files: electrode layout, adjacency matrix, list of passed fits, brms template |
| `inventory.csv` | every table and figure of the thesis, with the script that makes it |
| `archive/` | old code, notes and results, kept only for reference |

Generated data, fitted models, figures and logs are not in this repository.

## Order of the scripts

| Step | Script | Where it runs |
|---|---|---|
| 00 | `00_config.R` (paths and scenario map) | loaded by the other scripts |
| 01 | `01_generate_data.R` simulates the datasets | cluster or local |
| 02 | `02_fit_bgam.R` fits the BGAM of one seed (`COND` = I, H, F or E) | cluster, about 20 h per fit |
| 03 | `03_evaluate_conventional*.R` runs the conventional tests | cluster |
| 04 | `04_evaluate_bgam.R` turns each fit into detection masks and measures | cluster |
| 05 | `05_*.R` collects the per-seed results | cluster |
| 06 | `06_tables.R` makes the result tables | local |
| 07 | `07_figure_*.R` makes the figures | local |

Run every script from the repository root. Details, environment variables and the Slurm call are in [code/README.md](code/README.md).

## To reproduce only the tables and figures

Steps 06 and 07 need the aggregated results from step 05. Those are large and are not stored here.

## Software

R 4.4.2, brms 2.23.0, cmdstanr 0.8.0 with CmdStan 2.38.0, mgcv 1.9-4, posterior 1.7.0. The conventional tests use MNE-Python 0.23.4 through reticulate.
