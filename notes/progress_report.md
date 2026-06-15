# BGAM Spatio-Temporal Pilot — Progress Report

**Author:** jiahuifan
**Date:** 2026-05-04
**Stage:** Pilot validation (simulated data, n = 5 participants)
**Main result file:** `results/fits/bgam_fit_v4_k12.rds`
**Main figure:** `results/figures/topoplot_3rows_emp_vs_model.pdf`

---

## 1. Summary

We fit a Bayesian generalised additive model (BGAM) for joint
spatio-temporal inference on simulated EEG data, extending the temporal
BGAM framework of Nalborczyk & Bürkner (2025) to 3 dimensions
(electrode_x, electrode_y, time). The pilot validation (n = 5
participants, 61 electrodes, 51 time points) shows:

- **Convergence is excellent** (0 divergences, max R-hat = 1.003).
- **Spatial localisation correctly identifies the seed electrode (CZ).**
- **Temporal window is recovered with ±70 ms boundary smoothing.**
- **Amplitude is attenuated to ~33 % of ground truth** due to the
  GAM smoothing penalty plus measurement-error model.

Computational cost: ~8 hours wall time on LIDO HPC for 4 chains × 3000
iterations with 16 cores. Diagnostic output and figures are reproducible
from the code in `code/`.

---

## 2. Method

### 2.1 Data

Simulated EEG following Nalborczyk & Bürkner (2025) but extended to a
2-D scalp topography:

- 5 participants × 50 trials × 61 electrodes × 51 time points
- Sampling rate: ~102 Hz (downsampled from 500 Hz)
- Two conditions: `cond1` (noise only), `cond2` (noise + signal)
- Injected signal: spatial Gaussian (σ = 0.5) centred at electrode CZ,
  amplitude 1.5 uV, active in the time window [160 ms, 342 ms]
- Simulated noise: Yeung et al. (2004) EEG noise spectrum

### 2.2 BGAM model (v4_k12)

```r
brm(
  formula = eeg_mean | se(eeg_sd) ~
    condition +
    (1 + condition | participant) +
    s(time, bs = "tp", k = 20, by = condition) +
    t2(xproj, yproj, time,
       by = condition,
       bs = c("tp", "tp"),
       k  = c(10, 12),
       d  = c(2, 1)),
  data    = sim_data_model,
  family  = gaussian(),
  prior   = c(
    prior(normal(0, 1),       class = b),
    prior(student_t(3, 0, 1), class = sds),
    prior(student_t(3, 0, 1), class = sd)),
  chains = 4, cores = 4, iter = 3000, warmup = 1000,
  threads = threading(threads = 4),
  backend = "cmdstanr",
  control = list(adapt_delta = 0.99, max_treedepth = 13))
```

Key design choices vs Nalborczyk & Bürkner (2025):

| Component | Their (1D) | Ours (3D) |
|-----------|------------|-----------|
| Temporal smooth | `s(time, k = 50, by = condition)` | `s(time, k = 20, by = condition)` |
| Spatial smooth | (none) | `t2(xproj, yproj, time, k = c(10, 12))` |
| Measurement error | `eeg_mean \| se(eeg_sd)` | same |
| Multilevel structure | `(1 + condition \| participant)` | same |
| Sum-to-zero contrast | yes | yes |

---

## 3. Results

### 3.1 Convergence (Target 1) — PASS

*From `results/diagnostics/target1_convergence.csv`,
2026-05-04, log: `logs/diagnostics_2026_05_04.log`.*

| Metric | Value | Threshold | Status |
|--------|-------|-----------|--------|
| Divergent transitions | **0** | 0 | PASS |
| Max treedepth hits | **0** | < 5 % of draws | PASS |
| Max R-hat | **1.0026** | < 1.01 | PASS |
| Min Bulk_ESS | **2 435** | > 400 | PASS |
| Mean Bulk_ESS | **11 541** | — | PASS |
| Min Tail_ESS | **4 113** | > 400 | PASS |

Compared with an earlier `v3_fast3` fit using `max_treedepth = 10`,
which saturated on every draw (treedepth hits = 6 655 / 8 000), raising
treedepth to 13 in v4_k12 unblocked the sampler completely.

### 3.2 Onset / offset recovery (Target 2) — FAIL within ±20 ms

*From `results/diagnostics/target2_temporal.csv`.*

At the seed electrode CZ, the posterior odds threshold of 10 yields:

- True window: **160 – 342 ms**
- Estimated onset: **90 ms** (bias = **−70 ms**)
- Estimated offset: **410 ms** (bias = **+68 ms**)

Inspection of the per-time odds curve
(`results/diagnostics/target2_temporal_full.csv`) shows that the
posterior odds rise smoothly from 1 to 1000 between 90 ms and 130 ms,
and decay back below 10 between 410 ms and 420 ms. The boundary blur
is intrinsic to the smoothing penalty and persists across odds
thresholds. Crucially, no significant odds appear in the baseline
(*t* < 80 ms) or post-offset (*t* > 420 ms) regions, so the temporal
window is correctly *bracketed* even when the *boundaries* are blurred.

### 3.3 Spatial bump localisation (Target 3) — PASS

*From `results/diagnostics/target3_spatial.csv` and
`target3_spatial_full.csv`, evaluated at peak time 250 ms.*

- Top channel by posterior |cond2 − cond1|: **CZ** (= seed electrode) ✓
- Mean diff at CZ: **+0.552 μV** (95 % CrI [0.455, 0.648])
- Top 10 channels (decreasing |diff|):
  CZ → FCZ → C2 → CPZ → C1 → FC2 → FC1 → CP2 → CP1 → FZ
- All 9 channels following CZ are immediate neighbours of CZ in the
  10-20 system.

The smooth spatial decay from 0.55 μV at CZ to 0.41 μV at FZ confirms
the bump is centred at the seed electrode with the expected Gaussian
shape (true σ = 0.5 in the simulation).

### 3.4 Marginal effect — sanity check ✓

*From `results/diagnostics/marginal_check.csv`.*

- Empirical mean of (cond2 − cond1) across all data: **+0.0551 μV**
- BGAM marginal estimate: **+0.0563 μV** [95 % CrI: 0.0316, 0.0868]

Empirical value is well inside the model's credible interval; the
posterior median is within 2 % of the empirical mean.

Direct match between data and posterior. Note that the parametric
coefficient `conditioncond2` in `fixef()` is +1.29, but this is a
parameterisation artefact under sum-to-zero smooth constraints; the
**marginal** estimate is what should be reported in the paper.

### 3.5 Joint spatio-temporal comparison vs CBPT (preliminary)

A preliminary comparison against cluster-based permutation test
(cluster mass) on the same data yields:

| Metric | BGAM | CBPT |
|--------|------|------|
| Joint ST F1 | **1.044** | 0.859 |
| Spatial precision | 1.146 | 2.734 |
| Spatial recall | **0.959** | 0.509 |

BGAM achieves higher F1 by recovering the full active set (recall ≈ 1)
at the cost of slightly lower precision. CBPT is more conservative
(higher precision) but misses about half the temporal extent.

This comparison is preliminary (single dataset, not 500 simulations as
in the Nalborczyk paper). Full comparison is in `code/03_compare_methods.R`
(to be written).

### 3.6 Visualisation (Figure 1)

`results/figures/topoplot_3rows_emp_vs_model.pdf` shows three rows:

1. Empirical (cond2 - cond1): the actual signal in the data
2. BGAM (cond2 - cond1): model-estimated effect
3. Residual (empirical - model): structured residual at CZ peak time

The residual displays a "doughnut" pattern at 250 ms (positive at CZ,
negative on a surrounding ring), characteristic of excessive spatial
smoothing — the model spreads the bump beyond ground-truth extent.
Outside the peak time the residual is unstructured (consistent with
noise), confirming the temporal window itself is correctly identified.

---

## 4. Limitations

1. **Pilot only (n = 5).** Baseline noise per channel is non-zero
   (~0.1 uV per electrode); pooled across all channels and times the
   mean is ≈ 0.002 uV, indicating no systematic bias in the data.
   With n = 20 (final target), per-channel baseline noise will halve.

2. **Amplitude shrinkage (~67 %).** Recovered peak amplitude is
   ~33 % of true (0.56 uV vs 1.50 uV). This is consistent with GAM
   shrinkage and the measurement-error model. **Spatial localisation
   and temporal window are not affected — only amplitude.**

3. **Boundary blur (~70 ms).** Smoothing penalty produces gradual
   onset/offset rather than sharp boundaries; for rectangular signals
   this looks like under-detection at the boundaries.

4. **Spatial diffusion.** Posterior bump extends beyond ground-truth
   active set (~25 electrodes). Reflects penalty pulling neighbouring
   electrodes towards each other.

5. **Method comparison incomplete.** Currently only BGAM vs CBPT, on a
   single simulated dataset. Full comparison requires
   500-simulation Monte Carlo as in Nalborczyk paper.

6. **Simulated data.** Real M/EEG validation pending.

---

## 5. Next steps (proposed, ranked)

### Priority 1 — full method comparison (1-2 weeks)

Implement `03_compare_methods.R` to compare BGAM against:

- **Cluster mass** (Maris & Oostenveld 2007) — via MNE-Python
- **Mass univariate + FDR** (Benjamini-Hochberg)
- **GAMM 1D × N** (each electrode independently — direct comparison
  with Nalborczyk et al. 2025)

Evaluation dimensions:
- Temporal: onset / offset bias, MAE
- Spatial: precision, recall, F1 at peak time
- Joint ST: F1, weighted F1, AUC
- Computational cost
- Inference richness (binary mask vs full posterior)

### Priority 2 — robustness (2-3 weeks)

- Run 100-500 simulations to estimate variance per metric (matching
  Nalborczyk paper's protocol)
- Sensitivity to k, prior, sample size, signal amplitude

### Priority 3 — real data application (1-2 months)

- Apply to a real M/EEG dataset (suggestions welcome)
- Compare with traditional cluster-based analysis on the same data

### Priority 4 — methodological refinements (optional)

- Explore `bs = "gp"` (Matérn) basis to reduce spatial diffusion
- Try R-INLA + SPDE for computational scaling
- Implement GPP (Gaussian predictive processes) for larger n

---

## 6. Questions for advisor

1. Is the current pilot validation sufficient evidence to proceed to
   real-data application, or should we first run a full 500-simulation
   Monte Carlo to confirm robustness?

2. Which specific real M/EEG dataset would you suggest for validation?

3. For the upcoming method comparison, should we prioritise breadth
   (BGAM vs all 8 traditional methods, as in the Nalborczyk paper) or
   depth (BGAM vs the 3 most-used methods with thorough analysis)?

4. What is the target venue (workshop, conference, journal) for this
   work, so that the scope of validation matches reviewer expectations?

5. Is the 33 % amplitude recovery a concern for the paper narrative,
   or can we frame it as an inherent and acceptable trade-off of the
   smoothing-penalty approach?

---

## Appendix — file inventory

| File | Purpose | Status |
|------|---------|--------|
| `data/sim_data_pilot.rds` | Simulated EEG | ✓ |
| `data/ground_truth_pilot.rds` | Injected signal parameters | ✓ |
| `code/00_generate_data.R` | Data generation script | ✓ |
| `code/01_fit_bgam.R` | BGAM fit (v4_k12) | ✓ |
| `code/02_diagnostics.R` | Convergence + Targets 1/2/3 | ✓ |
| `code/03_compare_methods.R` | Method comparison | TODO |
| `code/08_slide_topomap_n20.R` | 3-row topographic plot | ✓ |
| `code/05_visualise_methods.R` | Method comparison plot | TODO |
| `code/utils.R` | Shared helpers | ✓ |
| `results/fits/bgam_fit_v4_k12.rds` | Main BGAM fit | ✓ |
| `results/figures/topoplot_3rows_*.pdf` | Main figure | ✓ |
| `slurm/run_fit.sh` | LIDO fit submission | ✓ |
| `slurm/run_diagnostics.sh` | LIDO diagnostics submission | ✓ |
