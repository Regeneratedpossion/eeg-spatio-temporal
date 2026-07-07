# Current BGAM Model Choice and Monte Carlo Plan

Source table: [`results/tables/ssh_model_run_summary.csv`](../../results/tables/ssh_model_run_summary.csv)

## Current SSH Model Run Summary

| Group | Run ID | P | Model | adapt | converged | divergent | treedepth | onset/offset | mean recovery | peak recovery | spatial F1 | ST F1 | hours | false positive |
|---|---:|---:|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| full BGAM | 39962984_1 | 5 | A_k1212_tp full BGAM | 0.80 | FALSE | 1 | 4285 | +10/-2 | 98.38 | 82.11 | 0.200 | 0.333 | 20.68 |  |
| full BGAM | 39963463_2 | 5 | B_k1212_ps full BGAM | 0.80 | FALSE | 2 | 0 | +10/-12 | 98.13 | 81.61 | 0.192 | 0.326 | 5.56 |  |
| full BGAM | 39962984_4 | 20 | B_k1212_ps full BGAM | 0.80 | TRUE | 0 | 244 | +10/-12 | 99.40 | 83.89 | 0.340 | 0.528 | 54.70 | 0.296 |
| adapt sensitivity | 39965039_1 | 5 | B_k1212_ps full BGAM | 0.95 | FALSE | 1 | 0 | +10/-12 | 98.15 | 81.61 | 0.192 | 0.325 | 10.77 |  |
| adapt sensitivity | 39965039_2 | 20 | B_k1212_ps full BGAM | 0.95 | TRUE | 0 | 140 | +10/-12 | 99.41 | 83.89 | 0.340 | 0.528 | 66.40 |  |
| basis sensitivity | 39986292 | 20 | **B_k1612_ps full BGAM** | 0.80 | TRUE | 0 | 22 | +10/-2 | 108.95 | 91.68 | 0.340 | 0.577 | 37.50 | 0.306 |
| low-rank Laplacian | 39963478-39963483 | 20 | Laplacian K20 + constant temporal BGAMM | 0.999 | FALSE | 16 | 0 | +10/-12 | 95.07 | 81.41 | 0.340 | 0.545 | 2.50 |  |

## Monte Carlo Repetition Plan

The current `B_k1612_ps` run took `37.50` wall-clock hours using `16` CPUs, or about `600` CPU-hours per fit. Based on this runtime and the current SSH job limits, a feasible Monte Carlo plan is:

```text
20 Monte Carlo repetitions, with at most 5 full BGAM jobs running at the same time.
```

### Monte Carlo Design Parameters

| Parameter | Current setting | Reason / note |
|---|---:|---|
| Full BGAM repetitions | 20 | Main planned batch |
| Seed set | 1--20 | One independent simulated dataset per seed |
| Maximum concurrent full BGAM jobs | 5 | Keeps the batch feasible under the current SSH / fair-share situation |
| Runtime per fit | about 37.50 wall-clock hours | Based on the current `B_k1612_ps` run |
| Expected batch duration | about 150 h plus queue time | If 5 jobs run in parallel |
| Fallback / extension | 10 or 30 repetitions | 10 if runtime becomes difficult; 30 only if the 20-run batch completes cleanly |

## BGAM Model Configuration

The planned Monte Carlo repetitions use the same full BGAM configuration as the current `B_k1612_ps` run.

```r
bgam_formula <- bf(
  eeg_mean | se(eeg_sem, sigma = TRUE) ~
    condition +
    (1 + condition | participant) +
    t2(
      xproj, yproj, time,
      by = condition,
      bs = c("ps", "ps"),
      k = c(16, 12),
      d = c(2, 1),
      full = FALSE
    )
)

fit <- brm(
  formula = bgam_formula,
  data = sim_data_model,
  family = gaussian(),
  chains = 4,
  cores = 4,
  iter = 3000,
  warmup = 1000,
  threads = threading(threads = 4),
  backend = "cmdstanr",
  control = list(
    adapt_delta = 0.80,
    max_treedepth = 10
  )
)
```

## Inference and Evaluation Parameters

The evaluation proceeds in three steps. First, the posterior condition effect is turned into a binary detection decision at each electrode-time cell. Second, a SESOI threshold calibrates this decision rule against baseline fluctuations, so that small noise-driven effects are not counted as detections. Third, recovery metrics compare the resulting detections and effect estimates with the known ground truth.

### Step 1: Posterior Detection Rule

This step turns the posterior distribution of the condition effect into a detected / not-detected decision.

$$
\Delta_{c,t} = \mathrm{cond2}_{c,t} - \mathrm{cond1}_{c,t}
$$

$$
p_{c,t} = \Pr(\Delta_{c,t} > \delta \mid \mathrm{data})
$$

$$
\mathrm{odds}_{c,t} = \frac{p_{c,t}}{1 - p_{c,t}}
$$

$$
\mathrm{detected}_{c,t} =
\begin{cases}
\mathrm{TRUE}, & \mathrm{if}\ \mathrm{odds}_{c,t} > 10,\\
\mathrm{FALSE}, & \mathrm{otherwise}.
\end{cases}
$$

Here, $\delta = 0$ is used for diagnostic plots. The threshold $\mathrm{odds}_{c,t} > 10$ corresponds to $p_{c,t} > 10/11$, or about $0.91$.

### Step 2: SESOI Correction

This step replaces the zero threshold with a baseline-derived SESOI threshold.

$$
\widehat{\delta}_{\mathrm{SESOI}}
= \max\left(0,\ \max_{c,\ t \in B}\widehat{\Delta}_{c,t}\right)
$$

where $B$ is the baseline window. The current setting uses $B = [0, 0.15]$ s.

The posterior probability is then computed as:

$$
p_{c,t}^{\mathrm{SESOI}}
= \Pr(\Delta_{c,t} > \widehat{\delta}_{\mathrm{SESOI}} \mid \mathrm{data})
$$

and the detection rule becomes:

$$
\mathrm{detected}_{c,t}^{\mathrm{SESOI}} = \mathrm{TRUE}
\quad \mathrm{if} \quad
\frac{p_{c,t}^{\mathrm{SESOI}}}{1 - p_{c,t}^{\mathrm{SESOI}}} > 10.
$$

The following check supports the current preference for the global SESOI rule:

| Rule | Baseline window | Practical-null FPR | Active post-offset rate | Active true-window rate |
|---|---:|---:|---:|---:|
| zero threshold | -- | 0.306 | 0.243 | 1.000 |
| global SESOI | 0.10 s | 0.060 | 0.000 | 0.994 |
| global SESOI | 0.15 s | 0.057 | 0.000 | 0.994 |
| electrode SESOI | 0.10 s | 0.125 | 0.208 | 1.000 |
| electrode SESOI | 0.15 s | 0.123 | 0.076 | 1.000 |

### Step 3: Recovery Metrics

This step compares the posterior detections and effect estimates with the known simulation truth.

| Recovery question | Code / metric | Truth used | Definition |
|---|---|---|---|
| Temporal onset/offset recovery | `target2_temporal` | true onset / offset | Estimated onset and offset bias at CZ; counted as a pass if the bias is <= `BIAS_TOLERANCE_MS = 20 ms` |
| Spatial recovery | `target3_spatial` | hard active-electrode label | F1 score for detected active electrodes at peak time |
| Spatio-temporal recovery | `target3b_spatiotemporal` | hard electrode-time mask | F1 score over electrode-time cells |
| Peak amplitude recovery | `target4_peak_amplitude` | continuous ground truth | Estimated peak effect at CZ divided by true peak amplitude |
| Window-mean recovery | `target4b_mean_amplitude` | continuous ground truth | Mean estimated effect in the true window divided by the true window mean |
| False-positive summary | `practical_null_fpr` | practical null region | Proportion of detected cells outside the practical signal region |

## Simulation Data-Generation Parameters

The current main dataset is `montecarlo_p20_gaussian_hardlabel`. Despite the folder name, the stored ground truth uses `signal_mode = continuous_label`: the injected signal is continuous, and hard labels are used only for evaluation.

### Data Structure

| Parameter | Symbol / code | Current value | Notes |
|---|---|---:|---|
| Participants | `n_ppt`, `P` | 20 | Main analysis sample size |
| Trials per participant per condition | `n_trials`, `R` | 20 | 20 trials for each participant and condition |
| Conditions | `condition` | 2 | `cond1` is noise only; `cond2` is signal plus independent noise |
| Electrodes | `C` | 61 | Reused from the pilot electrode layout |

### Time Axis

| Parameter | Symbol / code | Current value | Notes |
|---|---|---:|---|
| Original sampling rate | `srate` | 500 Hz | Data are generated on a 2 ms grid |
| Downsampling factor | `downsample_factor` | 5 | Analysis data are at 100 Hz |
| Analysis time grid | `time` | 0 to 0.50 s, step 0.01 s | 51 time points after downsampling |

### Injected Signal

| Parameter | Symbol / code | Current value | Notes |
|---|---|---:|---|
| Signal mode | `signal_mode` | `continuous_label` | Continuous Gaussian signal; hard labels are evaluation labels only |
| Signal centre | `seed_electrode` | CZ | Centre of the injected spatial signal |
| Peak signal amplitude | `bump_amplitude`, `A` | 1.5 | With baseline noise sd about 1.0, peak SNR is about 1.5 |
| Spatial signal width | `sigma_space` | 0.5 | Code and ground truth both use 0.5 |
| True onset | `true_onset` | 0.160 s | 160 ms |
| True offset | `true_offset` | 0.342 s | Last active analysis-grid point is 0.340 s |
| Participant onset latencies | `ronset` | Uniform `{150, 152, ..., 170}` ms | One onset latency per participant |

### Injected Noise

| Parameter | Symbol / code | Current value | Notes |
|---|---|---:|---|
| Temporal noise variance | `outvar` | 0.5 | EEG-like coloured temporal noise component |
| Spatial noise variance | `outvar` | 0.5 | Spatially correlated Gaussian noise component |
| Spatial noise correlation length | `rho` | 0.15 | Weak short-range spatial correlation |
| Total baseline noise scale | derived | about 1.0 sd | Two independent variance-0.5 components are added |

### Evaluation Ground Truth

| Parameter | Symbol / code | Current value | Notes |
|---|---|---:|---|
| Active-electrode threshold | `active_threshold` | mean + SD of raw spatial weights | Current threshold is about 0.397 |
| Active electrodes | `active_electrodes` | 9 electrodes | C1, C2, CP1, CP2, CPZ, CZ, FC1, FC2, FCZ |
| Continuous ground truth | `A * w_c * f(t)` | continuous field | Used for continuous recovery metrics |
| Hard-label ground truth | active electrodes x active time window | binary mask | Used for hard-label precision, recall, and F1 |

## Planned Model Comparison

The planned comparison uses three representative methods to compare different levels of spatio-temporal modelling.

| Method | Role | Reason |
|---|---|---|
| FDR-corrected mass-univariate test | Point-wise baseline | Controls multiple comparisons but does not borrow information across neighbouring electrodes or time points (Groppe et al., 2011) |
| Cluster-based permutation test | Standard EEG/MEG spatio-temporal baseline | Uses electrode-time adjacency to form clusters (Maris and Oostenveld, 2007) |
| Spatio-temporal BGAM | Proposed model-based approach | Estimates the full spatio-temporal effect surface with posterior uncertainty (Nalborczyk and Buerkner, 2025) |

This gives a simple comparison hierarchy: FDR does not model spatio-temporal neighbourhoods, cluster permutation uses neighbourhoods to define significant clusters, and BGAM directly models the condition effect as a smooth spatio-temporal surface.

For fair comparison, all methods will be converted to binary detected / not-detected maps over electrode-time cells. These maps will be compared with the hard-label ground truth using precision, recall, F1, false-positive rate, and onset/offset bias. Continuous recovery metrics, such as peak recovery and mean recovery, will be reported only for BGAM, because cluster-based and FDR methods are inference procedures rather than full effect-surface estimation models.
