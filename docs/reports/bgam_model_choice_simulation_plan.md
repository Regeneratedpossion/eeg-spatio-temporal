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

## Current Preferred Model

My current preferred full BGAM specification is:

```text
B_k1612_ps full BGAM, P = 20, sigmaTRUE, adapt_delta = 0.80
```

The main reasons are:

- It converged without divergent transitions.
- It had fewer treedepth hits than the earlier P20 `B_k1212_ps` runs.
- It improved temporal recovery: onset/offset bias was `+10/-2` ms, compared with `+10/-12` ms for the earlier P20 `B_k1212_ps` run.
- It improved peak recovery: `91.68%`, compared with `83.89%` for `B_k1212_ps`.
- It improved spatio-temporal F1: `0.577`, compared with `0.528` for `B_k1212_ps`.
- It was faster than the earlier P20 `B_k1212_ps` runs: `37.50 h` versus `54.70-66.40 h`.

One remaining concern is false positives. The practical-null false-positive rate is similar to the earlier P20 model:

- `B_k1212_ps`: `0.296`
- `B_k1612_ps`: `0.306`

For this reason, the `B_k1612_ps` model should be used together with the SESOI-based posterior decision rule when reporting detected regions.

## Monte Carlo Repetition Plan

The SSH check suggests that the model is computationally feasible but still expensive. The current `B_k1612_ps` run took `37.50` wall-clock hours using `16` CPUs.

Approximate compute cost per full BGAM fit:

```text
37.50 hours x 16 CPUs = about 600 CPU-hours
```

The current SSH account allows many submitted jobs, but fair-share priority still matters. I therefore suggest:

```text
Run 20 Monte Carlo repetitions of the final full BGAM, with at most 5 running at the same time.
```

In Slurm array notation:

```bash
#SBATCH --array=1-20%5
```

This means:

- `20` total Monte Carlo repetitions.
- At most `5` full BGAM jobs running at once.
- Expected compute cost: about `20 x 600 = 12,000 CPU-hours`.
- Expected wall-clock duration if 5 jobs run in parallel: about `20 / 5 x 37.5 = 150` hours, plus queue time.

A smaller version would be `10` repetitions. A larger version, such as `30` repetitions, is possible in principle, but I would not start with it unless the queue situation remains favourable and the 20-run batch completes cleanly.

## Interpretation

The full BGAM repetitions should be used for the main recovery and false-positive results.

Cheaper proxy checks can still be run with many more seeds, for example `100-1000`, to test specific mechanisms such as coloured temporal noise, spatial asymmetry, and basis-size sensitivity. These proxy checks do not replace the full BGAM; they explain whether the observed behaviour is systematic across many random seeds.
