# BGAM Spatio-Temporal EEG Inference

Reproducible simulation and Bayesian generalised additive modelling for a
known EEG effect over electrode space and time.

The formal thesis analysis compares three hard spatial truncation radii using
the same stochastic base in every replicate:

| Radius rule | Radius | Active electrodes |
|---|---:|---:|
| `1sd` | `sigma_space` | 5 |
| `half-max` | `sigma_space * sqrt(2 * log(2))` | 9 |
| `2sd` | `2 * sigma_space` | 21 |

All three datasets use high participant heterogeneity, 20 participants, 20
trials per participant and condition, and 20 paired replicates. The formal
BGAM uses P-spline margins, `k = c(16, 12)`, four chains, 3,000 iterations,
1,000 warmup iterations, and `eeg_mean | se(eeg_sem, sigma = TRUE)`.

Generated data, fitted models, figures, logs, and CmdStan caches are not stored
in GitHub.

## Repository structure

```text
code/
|-- 00_...R to 05_...R   main Monte Carlo pipeline
|-- 06_...R to 10_...R   paired-radius pilot
|-- 11_...R to 15_...R   BGAM and SESOI diagnostics
|-- 16_...R to 21_...R   figures and topomaps
|-- 22_...R to 27_...R   sensitivity analyses
|-- *.R                  sourced helper modules
`-- archive/              superseded or one-off scripts

slurm/
|-- formal/          paired generation, fitting, and aggregation
`-- pilot/           one-replicate paired-radius run
```

See [code/README.md](code/README.md) for the numbered execution order. Older material is in `archive/`.

## Formal workflow

Run commands from the repository root.

Generate 20 paired replicates. Each seed is set once before the common base is
created; only the hard spatial radius differs across the three saved datasets.

```bash
MC_N_PPT=20 \
MC_SCENARIO=high \
Rscript code/01_generate_data.R 20
```

This writes:

```text
data/montecarlo_rhalfmax/
```

Fit one radius and replicate locally:

```bash
Rscript code/07_fit_paired_radius_bgam.R 1sd 1
Rscript code/07_fit_paired_radius_bgam.R half-max 1
Rscript code/07_fit_paired_radius_bgam.R 2sd 1
```

The second argument is the replicate seed. Omitting it preserves the pilot
default of seed 1.

Cluster entry points are:

```bash
sbatch slurm/formal/generate_paired_radius.sh
sbatch slurm/formal/fit_paired_radius_hpc2n.sh
sbatch slurm/formal/fit_paired_radius_dortmund.sh
sbatch slurm/formal/aggregate_paired_radius.sh
```

Both fit arrays contain 60 tasks: three radius rules times 20 replicates. The
model settings are identical across radii and clusters.

## Ground truth

The spatial Gaussian is hard-truncated. Electrodes outside the selected radius
have exactly zero weight and zero injected effect. Active electrodes are
derived as those whose truncated spatial weight is greater than zero.

Participant temporal waveforms are finite, zero-padded templates shifted by
participant onset. They are not clipped a second time to a common temporal
window. The fixed 160 and 342 ms values are reference timing estimands; on the
100 Hz grid the scoring times are 0.16 through 0.34 s.

Two continuous truth fields are distinguished:

- theoretical population truth, fixed across replicates;
- realised finite-sample truth, determined by participants in one dataset.

Spatial FP/FPR at 250 ms is an exact-zero metric. A full-grid detection outside
the timing-scoring window is a label-based false positive and may occur where
the population temporal tail is small but non-zero.

## Supporting analyses

Supporting scripts are retained for audit but are not part of the formal
execution sequence:

- `00`--`05` contain the main Monte Carlo pipeline;
- `06`--`10` contain the paired-radius pilot;
- `11`--`15` contain BGAM and SESOI diagnostics;
- `16`--`21` contain figures and topomaps;
- `22`--`27` contain prior, basis, Laplacian, and threshold sensitivity work;
- unnumbered files in `code/` are sourced helper modules;
- `code/archive/` contains superseded single-radius or one-off scripts.

Obsolete P=5/P=10, single-radius, `k = c(12, 12)`, basis-rescue, and Laplacian
Slurm launchers were removed. Their history remains available through Git.
