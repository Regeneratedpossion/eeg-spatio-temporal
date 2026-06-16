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
├── slides/
│   └── thesis_progress_slides.tex      # compact current thesis-progress slides
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

## Key Results

The completed slide-level results currently come from the smaller pilot
setting: `P = 5` participants and `R = 20` trials per condition. The full
`P = 20` hard-threshold SEM Monte Carlo runs are still running on the cluster.

### Pilot Model Comparison

The first pilot analysis compared measurement-error choices (`se(eeg_sd)` vs
`se(eeg_sem)`) and several basis-function settings. The best pilot fit was the
SEM model:

```text
k=(20,20), full=FALSE, se(eeg_sem)
onset bias  = +0.010 s
offset bias = -0.012 s
peak recovery = 89%
mean-window recovery = 106.1%
spatial F1 = 0.962
wall clock = 67 h
```

The SD-based fits were useful diagnostics, but they had larger onset/offset
bias and lower spatial F1. This is why the current main pipeline uses
`se(eeg_sem)`.

### Pilot Method Comparison

Using the selected BGAM specification, the pilot Monte Carlo comparison ran 30
independent `P = 5` simulations. BGAM_2D had the best offset recovery and the
best full spatio-temporal F1:

```text
BGAM_2D onset:  bias -0.011 s, MAE 0.011 s, RMSE 0.012 s, SD 0.005 s
BGAM_2D offset: bias +0.012 s, MAE 0.012 s, RMSE 0.014 s, SD 0.006 s
BGAM_2D F1_ST:  0.799 +/- 0.007
```

Raw pointwise p-values had competitive onset MAE, but produced many false
positives across the full time x channel mask:

```text
Raw p:   FP rate 0.102, FN rate 0.530, F1_ST 0.462, mean FP cells 268.9
Cluster: FP rate 0.003, FN rate 0.648, F1_ST 0.514, mean FP cells 8.0
TFCE:    FP rate 0.007, FN rate 0.608, F1_ST 0.539, mean FP cells 18.4
BGAM_2D: FP rate 0.042, FN rate 0.179, F1_ST 0.799, mean FP cells 111.3
```

Interpretation: BGAM_2D recovers most true signal cells with substantially
fewer false negatives than cluster/TFCE, but it also smooths beyond the strict
ground-truth mask, creating moderate false-positive spillover.

### Pilot Appendix Figures

The appendix pilot figures show the earlier SD-based diagnostic fit:

```text
P = 5, R = 20
k = (20,20), full = FALSE
measurement error = se(eeg_sd)
```

These figures show simulated EEG data and fitted BGAM curves against raw
participant means. They are kept for pilot diagnostics; the current main model
uses `se(eeg_sem)`.

### Alpha Threshold Note

The completed `P = 5` pilot FP/FN diagnosis used `alpha = 0.10`. With only five
participants, the permutation p-value resolution is coarse, so `0.05` is too
strict for a useful pilot comparison.

The current `P = 20` cluster runs use `alpha = 0.05`, which is appropriate once
the participant count is larger.

### Cluster Status Snapshot

Checked on `gw01` at `2026-06-16 14:36 CEST`. The active jobs are still in
`Step 2/3: fit BGAM`; no `sim_*.csv` comparison files or
`bgam_summary_*.rds` files had been written yet.

```text
bgam_mc_n20_hard_sem
  30 array tasks, P=20 hard-threshold SEM, k=(20,20), full=FALSE
  state: RUNNING on ultralong
  runtime: about 2 days 21 hours; time limit: 10 days
  sampled progress: still in warmup, around 100-200 / 3000 iterations (3-6%)

bgam_mc_n20_hard_sem_k1212
  3 array tasks, P=20 hard-threshold SEM, k=(12,12), full=FALSE
  state: RUNNING on ultralong
  runtime: about 2 days to 2 days 20 hours; time limit: 10 days
  sampled progress: still in warmup, around 300-600 / 3000 iterations (10-20%)

bgam_mc_p5_hard_sem
  1 task, P=5 hard-threshold SEM, k=(20,20), full=FALSE
  previous long-queue job was cancelled; first ultralong resubmission failed
  because the formula could not find k_xy/k_t
  code fixed and resubmitted as job 39882029
  state: PENDING on ultralong; time limit: 10 days

bgam_mc_p10_hard_sem
  1 task, P=10 hard-threshold SEM, k=(20,20), full=FALSE
  previous long-queue job timed out after 2 days at about 200 / 3000 warmup
  resubmitted as job 39882034
  state: PENDING on ultralong; time limit: 10 days

bgam_n20_r20_reff14
  single P=20, R=20 effective-R sensitivity fit with R_eff=14
  state: RUNNING on ultralong
  runtime: about 6 days 10 hours; time limit: 10 days
  sampled progress: still in warmup, around 600-700 / 3000 iterations (20-23%)
```

The SSH logs for the currently running Monte Carlo jobs show `alpha = 0.050`.
The pilot FP/FN slide above reports the completed pilot diagnosis at
`alpha = 0.10`.

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
