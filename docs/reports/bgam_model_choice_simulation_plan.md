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

| Category | Parameter | Symbol / code | Current value | Status | Notes |
|---|---|---|---:|---|---|
| Data size | Participants | `n_ppt`, `P` | 20 | confirmed | Main analysis data: `montecarlo_p20_gaussian_hardlabel` |
| Data size | Trials per participant per condition | `n_trials`, `R` | 20 | confirmed | 20 trials for each participant and condition |
| Data size | Electrodes | `C` | 61 | confirmed | Reused from the pilot electrode layout |
| Time axis | Original sampling rate | `srate` | 500 Hz | confirmed | Data are generated on a 2 ms grid |
| Time axis | Downsampling factor | `downsample_factor` | 5 | confirmed | Analysis data are at 100 Hz, with a 10 ms grid |
| Conditions | Number of conditions | `condition` | 2 | confirmed | `cond1`, `cond2` |
| Signal mode | Main data mode | `signal_mode` | `continuous_label` | confirmed | The injected signal is continuous; hard labels are used only for evaluation |
| Signal | Peak signal amplitude | `bump_amplitude`, `A` | 1.5 | confirmed | With baseline noise sd about 1.0, peak SNR is about 1.5 |
| Signal | Spatial signal width | `sigma_space` | 0.5 | confirmed | Code and ground truth both use 0.5 |
| Noise | Spatial noise correlation length | `rho` | 0.15 | confirmed | Weak short-range spatial correlation |
| Noise | Variance of each noise component | `outvar` | 0.5 | confirmed | Temporal and spatial noise components are independent and added together |
| Time truth | True onset | `true_onset` | 0.160 s | confirmed | 160 ms |
| Time truth | True offset | `true_offset` | 0.342 s | confirmed | Last active analysis-grid point is 0.340 s |
| Space truth | Seed electrode | `seed_electrode` | CZ | confirmed | Centre of the injected spatial signal |
| Participant variation | Onset-latency range | `ronset` | Uniform `{150, 152, ..., 170}` ms | confirmed | One onset latency per participant |
| Evaluation label | Active-electrode threshold | `active_threshold` | mean + SD of raw spatial weights | confirmed | Current threshold is about 0.397 |
| Evaluation label | Active electrodes | `active_electrodes` | 9 electrodes | confirmed | C1, C2, CP1, CP2, CPZ, CZ, FC1, FC2, FCZ |
