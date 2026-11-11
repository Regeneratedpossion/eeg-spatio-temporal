# Spatio-temporal EEG modelling with Bayesian GAMs

Code for a master's thesis at TU Dortmund University. The study compares a Bayesian generalised additive model (BGAM) with six conventional methods for detecting **when and where** simulated EEG effects occur.

## Study

We simulate 100 datasets in each of four scenarios, with 20 participants and 61 electrodes per dataset.

| Scenario | Effect amplitude (µV) | Between-participant heterogeneity |
|---|---:|---:|
| 1 | 0.5 | 0.26 |
| 2 | 0.5 | 0.70 |
| 3 | 1.5 | 0.26 |
| 4 | 1.5 | 0.70 |

Of 400 BGAM fits, **362 meet the convergence criteria**. The analysis compares detection accuracy, timing error, false positives and effect recovery.

## Repository structure

```text
code/           Simulation, fitting, evaluation, tables and figures (steps 01–07)
code/R/         Shared R functions
slurm/          Jobs for the LIDO cluster
resources/      Electrode information, adjacency and the brms template
data/           Selected ground-truth data
results/        Selected results needed for figures and tables
inventory.csv   Which scripts produce each thesis figure or table
```

## Reproduce tables and figures

Run commands from the repository root. Using the stored results does **not** require fitting all BGAMs again.

```bash
Rscript code/06_tables/06_tables.R passed
```

The figure scripts are in `code/07_figures/`. See `inventory.csv` for their inputs and outputs. Full simulation and model fitting use `code/01_simulate/` through `code/05_aggregate/` and the scripts in `slurm/`.

See [Reproducibility notes](docs/REPRODUCIBILITY.md) and [Code instructions](code/README.md).

## Software and related work

The original analysis used **R 4.4.2**, `brms` 2.23.0, CmdStan 2.38.0 and **MNE-Python 0.23.4**.

This work builds on [Nalborczyk and Bürkner's M/EEG research](https://github.com/lnalborczyk/brms_meeg). Their [neurogam](https://github.com/lnalborczyk/neurogam) package provides related tools.
