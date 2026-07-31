# Paper Repository Cleanup Design

## Status and precedence

This design supersedes the compatibility-retention portions of
`2026-07-31-paper-reproducibility-code-organization-design.md`. The approved
scientific contract from that document remains unchanged. In particular, the
formal analysis remains the paired-radius, high-heterogeneity, P=20 analysis
with BGAM `k = c(16, 12)`, P-spline margins, and `sigma = TRUE`.

The cleanup deliberately removes obsolete entry points after all live
references have been redirected. Git history remains the recovery mechanism
for deleted historical launchers.

## Goal

Make the GitHub repository present the thesis replication workflow first,
retain useful supporting research code in purpose-specific directories, and
delete obsolete or duplicated launch scripts. The cleanup must not modify
simulated observations, fitted models, scientific estimands, model formulas,
priors, SESOI, posterior-odds rules, cluster-mass, TFCE, or detection rules.

## Target code structure

```text
code/
|-- config/
|   `-- paper_config.R
|-- pipeline/
|   |-- simulation.R
|   |-- generate_montecarlo.R
|   |-- fit_bgam.R
|   |-- fit_paired_radius.R
|   |-- compare_methods.R
|   `-- aggregate_results.R
|-- truth/
|   |-- ground_truth.R
|   `-- migrate_ground_truth.R
|-- reporting/
|-- diagnostics/
|-- sensitivity/
|-- pilot/
|-- archive/
|-- utils.R
`-- check_env.R
```

The repository README starts with the formal workflow and links to one
canonical command for generation, fitting, evaluation, and aggregation.
Supporting code is documented separately and is not presented as part of the
paper's required execution path.

## Formal pipeline

The formal modules contain reusable functions and do not run work when
sourced. Their public contract is driven by `paper_analysis_config()`:

- paired radius rules `1sd`, `half-max`, and `2sd`;
- high between-participant heterogeneity;
- 20 participants, 20 trials, and 20 replications;
- one stochastic base shared by all three radii in each replication;
- hard spatial truncation with 5, 9, and 21 active electrodes;
- explicit theoretical and realised continuous ground truth;
- reference timing estimands of 160 and 342 ms, represented on the 100 Hz
  grid as 0.16 through 0.34 s;
- BGAM P-spline margins with `k = c(16, 12)`, four chains, 3,000 iterations,
  1,000 warmup iterations, and `sigma = TRUE`.

Small top-level entry points may remain only when they are the documented
formal commands. A moved diagnostic, pilot, sensitivity, reporting, or
historical implementation does not retain a duplicate top-level adapter.
All repository-owned callers are updated to the canonical path instead.

## Supporting-code classification

### Reporting

Move thesis figures, slide figures, topomaps, and posterior visualisations to
`code/reporting/`, including the current `05`--`08` slide scripts,
`07_plot_paired_radius_truth_pilot.R`,
`31_plot_all_electrode_sesoi_posterior.R`, and
`37_plot_2sd_bgam_diagnostics.R`.

### Diagnostics

Move fit checks and methodological investigations to `code/diagnostics/`,
including `02_diagnostics.R`, `19_practical_null_fpr_rebound.R`,
`28_delta_sesoi_probe.R`, `29_validate_heterogeneity_bgam.R`, and
`35_diagnose_paired_radius_bgam_pilot.R`.

### Sensitivity

Move prior, basis, threshold, alpha, SESOI, and Laplacian investigations to
`code/sensitivity/`. The existing `code/laplacian/` implementation becomes a
subdirectory of this section. These files remain available for scientific
audit but are not part of the formal reproduction command sequence.

### Pilot

Move reduced-scale paired-radius generation, validation, fitting, and
aggregation tools to `code/pilot/`. Pilot outputs and commands are labelled as
pilot results in documentation.

### Archive

Move the superseded single-radius implementation
`09_run_montecarlo_sem.R` and one-off exploratory calculations such as
`compute_trial_effective_sample_size.R` to `code/archive/`. Archive files are
retained as readable historical source but receive no compatibility wrapper
and no active Slurm launcher.

## Slurm cleanup

Only a small set of current launchers remains:

```text
slurm/
|-- formal/
|   |-- generate_paired_radius.sh
|   |-- fit_paired_radius_hpc2n.sh
|   |-- fit_paired_radius_dortmund.sh
|   `-- aggregate_paired_radius.sh
`-- pilot/
    `-- run_paired_radius_pilot.sh
```

The formal launchers call canonical pipeline modules, accept radius rule and
replication seed explicitly where applicable, and reproduce the central paper
configuration. Cluster-specific differences are limited to scheduler and
filesystem settings; they do not change model settings.

The following historical launcher families are deleted after their useful R
implementations have been classified:

- P=5 and P=10 Monte Carlo test launchers;
- old single-radius P=20 Monte Carlo launchers;
- the `k = c(12, 12)` test launcher;
- low-heterogeneity and `gaussian_hardlabel` generation launchers;
- obsolete basis-screening rescue and scaling launchers;
- Laplacian batch launchers;
- generic early `run_fit` and diagnostic launchers;
- temporary slide-generation and original-diagnostic launchers; and
- duplicated paired-radius `r1sd`, `fit_only`, and cluster-copy launchers.

No currently running Umea or Dortmund job is cancelled, altered, or
resubmitted. Remote cluster working directories are not changed as part of
this repository cleanup.

## Deletion safety

Before deleting a file, the implementation must establish all of the
following:

1. The file is not a canonical formal-pipeline module.
2. Any useful scientific implementation has been moved to an appropriate
   supporting directory or already exists elsewhere.
3. `rg` finds no remaining repository-owned call to the old path after caller
   updates.
4. The file is not needed by a currently documented formal command.
5. The resulting R and shell trees pass syntax and regression checks.

Data, fitted-model objects, result files, SSH configuration, logs, cluster
working copies, and thesis source are outside the deletion scope. Deletions
are recorded as normal Git changes and are therefore recoverable from history.

## Documentation

`README.md` becomes the short public entry point for formal reproduction.
`PROJECT_GUIDE.md` documents module responsibilities, data flow, outputs, and
the distinction between formal, pilot, diagnostic, sensitivity, reporting,
and archived code. Historical P=5, P=10, low-heterogeneity, single-radius, and
`k = c(12, 12)` settings are not described as current analysis settings.

The README and guide contain no commands pointing to deleted top-level or
Slurm paths.

## Verification

The cleanup is accepted only when all of the following pass:

- an inventory test asserts the canonical classification of every retained R
  script;
- every retained R file parses recursively;
- every retained Slurm file passes `bash -n`;
- `rg` finds no stale reference to deleted paths;
- all `tests/test_*.R` scripts pass with `PROJECT_ROOT` set;
- the thesis contract test passes against the local thesis source;
- a temporary paired generation produces 5, 9, and 21 active electrodes,
  exact spatial zeros, identical paired residuals, and explicit theoretical
  and realised truth;
- Git diff contains no data, fitted-model, results, SSH, log, or unrelated
  user-file change; and
- no remote push occurs without explicit user approval.

## Delivery

Implementation occurs on the isolated branch
`refactor/paper-reproducibility`. Commits are grouped by formal modules,
supporting-code moves, Slurm deletion, documentation, and final verification.
The final handoff includes a retained/moved/deleted file table and the exact
test evidence. The user reviews the branch before any push or merge.
