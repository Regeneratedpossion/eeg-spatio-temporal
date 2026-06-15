# BGAM Spatio-Temporal EEG Inference

Bayesian generalised additive modelling for detecting simulated EEG effects
over electrode space and time. The current main pipeline uses a hard-threshold
spatial signal, trial-mean SEM measurement error, and a 2D spatio-temporal BGAM.

The repository is intentionally lightweight. Generated `.rds` data, fitted
models, figures, logs, and CmdStan cache files are not uploaded to GitHub.

---

## Repository Layout

```text
bgam_eeg_project/
├── PROJECT_GUIDE.md                    # compact guide to the current pipeline
│
├── code/
│   ├── 00_generate_data.R              # one simulated EEG dataset
│   ├── 00_generate_montecarlo.R        # 30 Monte Carlo datasets, P=20, R=20
│   ├── 01_fit_bgam.R                   # SEM BGAM model
│   ├── 02_diagnostics.R                # single-fit diagnostics
│   ├── 03_compare_methods.R            # single-simulation method comparison
│   ├── 04_aggregate_results.R          # final Monte Carlo aggregation
│   ├── 05_slide_simulated_eeg_n20.R    # slide figure: simulated data
│   ├── 06_slide_bgam_fit_n20.R         # slide figure: fitted curves
│   ├── 07_slide_posterior_probability_n20.R
│   ├── 08_slide_topomap_n20.R
│   ├── 09_run_montecarlo_sem.R         # one Slurm array task
│   ├── 10_prior_predictive_check.R
│   └── utils.R                         # paths, thresholds, shared helpers
│
├── resources/
│   ├── adjacency/                      # 61-channel adjacency matrix
│   └── brms_meeg/code/                 # noise generator and ERP template
│
├── slurm/                              # cluster submission scripts
└── refs.bib
```

Generated folders are local only:

```text
data/                 # generated simulations and ground truth
results/fits/         # fitted brms models
results/montecarlo/   # BGAM posterior summaries
results/comparisons/  # per-simulation and aggregated method results
results/figures/      # slide and diagnostic figures
cmdstan_cache/        # compiled Stan files
```

---

## Main Pipeline

Run commands from the project root. `PROJECT_ROOT` is optional because scripts
fall back to the current working directory.

```bash
# 1. Generate 30 hard-threshold Monte Carlo datasets, P=20 and R=20.
MC_DATA_SUBDIR=montecarlo_n20_hard Rscript code/00_generate_montecarlo.R

# 2. Run one simulation locally, or submit the Slurm array version.
MC_DATA_SUBDIR=montecarlo_n20_hard \
MC_COMPARE_SUBDIR=montecarlo_n20_hard_sem \
MC_SUMMARY_SUBDIR=montecarlo_n20_hard_sem \
Rscript code/09_run_montecarlo_sem.R 1

# 3. Aggregate all finished simulations.
MC_DATA_SUBDIR=montecarlo_n20_hard \
MC_COMPARE_SUBDIR=montecarlo_n20_hard_sem \
MC_SUMMARY_SUBDIR=montecarlo_n20_hard_sem \
MC_AGG_PREFIX=montecarlo_n20_hard_sem \
Rscript code/04_aggregate_results.R
```

On the cluster, use:

```bash
sbatch slurm/gen_montecarlo_n20_hard.sh
sbatch slurm/run_montecarlo_n20_hard_sem.sh
```

For a smaller basis-function test, use:

```bash
K_XY=12 K_T=12 Rscript code/09_run_montecarlo_sem.R 1
```

---

## Simulation

`code/00_generate_data.R` creates two-condition EEG data with:

- participants: `P = 20` in the Monte Carlo generator;
- trials per participant: `R = 20`;
- sampling rate: 500 Hz, downsampled to 100 Hz;
- true signal window: 160-342 ms;
- seed electrode: `CZ`;
- spatial signal: Gaussian around `CZ`, hard-thresholded at one spatial SD;
- noise: temporal EEG-like noise plus spatially correlated sensor noise.

Ground truth stores the true onset, offset, seed electrode, active electrodes,
signal amplitude, and simulation settings used for evaluation.

---

## Model Specification

The current BGAM in `code/01_fit_bgam.R` is:

```r
bgam_fit <- brm(
    formula = eeg_mean | se(eeg_sem) ~
        condition +
        (1 + condition | participant) +
        t2(xproj, yproj, time,
           by = condition,
           bs = c("tp", "tp"),
           k = c(k_xy, k_t),
           d = c(2, 1),
           full = full),
    data = sim_data_model,
    family = gaussian(),
    chains = 4,
    cores = 4,
    iter = 3000,
    warmup = 1000,
    seed = seed,
    threads = threading(threads = 4),
    backend = "cmdstanr",
    control = list(adapt_delta = 0.99, max_treedepth = 13)
)
```

Important current choices:

- The response is the trial-aggregated mean `eeg_mean`.
- Measurement error uses `se(eeg_sem)`, not `se(eeg_sd)`.
- There is no separate `s(time)` term in the current main model.
- The spatio-temporal smooth is a tensor product over `(xproj, yproj, time)`.
- Default basis size is `k = c(20, 20)`.
- `K_XY` and `K_T` can change the basis size without editing the script.

---

## Method Comparison

`code/03_compare_methods.R` compares six methods for each single simulation:

1. raw pointwise t-test;
2. Benjamini-Hochberg FDR;
3. Holm correction;
4. MNE spatio-temporal cluster permutation;
5. MNE TFCE;
6. BGAM posterior odds on the full channel x time grid.

Each method produces a time x channel detection mask. Small isolated blobs are
removed, then the largest connected component gives the estimated onset and
offset.

`code/04_aggregate_results.R` then summarises across simulations:

- onset/offset bias;
- MAE;
- RMSE;
- SD;
- spatio-temporal F1;
- TP/FP/FN diagnostics.

---

## Key Result

Current local 30-simulation summary, generated from the hard-threshold SEM
pipeline:

| Method | Onset MAE (s) | Offset MAE (s) | F1_ST |
|---|---:|---:|---:|
| BGAM_2D | 0.011 | 0.012 | 0.799 +/- 0.007 |
| raw_p | 0.010 | 0.013 | 0.462 +/- 0.014 |
| cluster | 0.016 | 0.019 | 0.514 +/- 0.020 |

Interpretation: raw pointwise tests can be competitive for temporal boundary
MAE in this controlled simulation, but BGAM gives much better spatio-temporal
localisation on the full channel x time grid.

---

## Slide Figures

The slide scripts use the same hard-threshold + SEM defaults via
`slide_inputs()` in `utils.R`.

```bash
Rscript code/05_slide_simulated_eeg_n20.R
Rscript code/06_slide_bgam_fit_n20.R
Rscript code/07_slide_posterior_probability_n20.R
Rscript code/08_slide_topomap_n20.R
```

They create:

- simulated EEG curves;
- BGAM fitted curves and condition differences;
- posterior probability per electrode;
- ground truth / empirical / BGAM / residual topomaps.

---

## Prior Predictive Check

```bash
Rscript code/10_prior_predictive_check.R
```

This runs `sample_prior = "only"` for the current hard-threshold + SEM BGAM
formula. A weak normal prior is supplied only because prior predictive sampling
requires a proper fixed-effect prior.

---

## Dependencies

Cluster environment used for the main Monte Carlo runs:

- TU Dortmund LIDO cluster;
- Slurm array jobs;
- `module load R/4.4.2-gcc141-base`;
- CmdStan backend through `cmdstanr`;
- 16 CPU threads and 128 GB memory in the current Slurm scripts.

R packages:

- `brms`
- `cmdstanr`
- `dplyr`
- `tidyr`
- `tidybayes`
- `posterior`
- `ggplot2`
- `patchwork`
- `reticulate`
- `igraph`

Python package:

- `mne` for cluster permutation and TFCE methods.

---

## References

- Nalborczyk, L. & Bürkner, P.-C. (2025). *Precise temporal localisation
  of M/EEG effects with Bayesian generalised additive multilevel models*.
  bioRxiv. doi:10.1101/2025.08.29.672336
