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
