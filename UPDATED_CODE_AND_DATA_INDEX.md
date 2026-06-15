# Updated code and data index

This file collects the updated scripts and the data files needed to understand
the current BGAM EEG pipeline.

## 1. Core updated code

### Data generation

Path:

```text
code/00_generate_data.R
```

Purpose:

- Defines the 2D EEG simulation.
- Default setting is the pilot dataset: `P = 5`, `R = 20`.
- Generates:
  - condition 1: noise only
  - condition 2: noise + CZ-centred ERP signal
- Defines the ground truth:
  - temporal window: `0.160--0.342 s`
  - seed electrode: `CZ`
  - spatial Gaussian bump

Important function:

```r
make_dataset(seed = 42, save_rds = TRUE)
```

Direct run:

```bash
PROJECT_ROOT=$(pwd) Rscript code/00_generate_data.R
```

Outputs:

```text
data/sim_data_pilot.rds
data/ground_truth_pilot.rds
```

### Pilot preview plots

Path:

```text
code/plot_pilot_preview.R
```

Purpose:

- Reads or generates the pilot data.
- Plots quick preview figures for checking the simulated data.

Direct run:

```bash
PROJECT_ROOT=$(pwd) Rscript code/plot_pilot_preview.R
```

Outputs:

```text
figures/sim_data_topo_both.png
figures/spatial_ground_truth.png
```

## 2. Monte Carlo data

### P=20 generation

Path:

```text
code/00_generate_montecarlo.R
```

Purpose:

- Generate Monte Carlo datasets from the current hard-threshold data generator.
- Use `MC_N_PPT=20` for `P = 20`.
- Use `MC_DATA_SUBDIR` to choose the output directory.

Example:

```bash
PROJECT_ROOT=$(pwd) MC_N_PPT=20 MC_DATA_SUBDIR=montecarlo_n20_hard \
    Rscript code/00_generate_montecarlo.R 30
```

Main output files:

```text
data/montecarlo_n20_hard/sim_001.rds
data/montecarlo_n20_hard/gt_001.rds
```

These are the data used for the main model:

```text
P = 20 participants
R = 20 trials per participant per condition
C = 61 electrodes
K = 2 conditions
T = 51 time points
```

## 3. Model fitting

Path:

```text
code/01_fit_bgam.R
```

Purpose:

- Fits the spatiotemporal BGAM.
- Aggregates trials to participant x channel x condition x time.
- Uses the model:

```r
eeg_mean | se(...) ~
    condition +
    (1 + condition | participant) +
    t2(xproj, yproj, time,
       by = condition,
       k = c(k_xy, k_t),
       full = full)
```

Current main setting:

```text
P = 20
R = 20
k = c(20, 20)
SE mode = sd
full = FALSE
```

Main fit file on SSH:

```text
results/fits/bgam_fit_v5_k20_n20_r20_sd.rds
```

Slurm script:

```text
slurm/run_bgam_n20_r20_sd.sh
```

## 4. Diagnostics

Path:

```text
code/02_diagnostics.R
```

Purpose:

- Runs the full diagnosis for one fitted model.
- Produces:
  - convergence diagnostics
  - onset / offset recovery
  - spatial localisation
  - peak amplitude recovery
  - window-mean amplitude recovery
  - marginal condition difference check

Current P=20 diagnostic outputs:

```text
results/diagnostics/june_full/
```

Key result:

```text
P=20, R=20, k=(20,20), SD:
  0 divergences
  onset/offset around 150--350 ms
  top channel CZ
  spatial F1 high
  window-mean amplitude close to truth
  pointwise peak underestimated
```

## 5. P=20 topomap

Path:

```text
code/plot_topomap_n20_r20.R
```

Purpose:

- Plots topomaps at three time points:
  - `160 ms`: onset
  - `250 ms`: peak
  - `340 ms`: offset
- Rows:
  - ground truth
  - empirical difference: `cond2 - cond1`
  - BGAM estimate
  - residual: `empirical - model`

Output:

```text
results/figures/topomap_n20_r20_three_times.pdf
```

This is the figure currently used in the slide.

## 6. Monte Carlo / method comparison

Paths:

```text
code/03_compare_methods.R
code/run_montecarlo.R
code/04_aggregate_results.R
```

Purpose:

- Compares BGAM with:
  - raw p
  - BH
  - Holm
  - cluster permutation
  - TFCE
- BGAM detection uses posterior odds mask:

```text
posterior odds > 10
```

P=20 Monte Carlo Slurm script:

```text
slurm/run_montecarlo_n20.sh
```

Current P=20 30-MC job:

```text
Job ID: 39877493
Job name: bgam_mc_n20
Status: pending/running depending on cluster resources
```

## 7. Slides

Main slide source:

```text
thesis_progress_slides.tex
```

Compiled PDF:

```text
thesis_progress_slides.pdf
```

Current main narrative:

```text
1. Data setup
2. Ground truth
3. P=20 model fit
4. P=20 posterior probability
5. P=20 model result
6. P=20 topomaps over signal window
7. SE calibration / R_eff
8. Appendix: metrics, mask diagnostics, methods comparison
```
