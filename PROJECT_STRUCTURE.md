# Project Structure

This project is organized around the active BGAM simulation and model-fitting pipeline.

## Main Working Files

```text
code/
```

Core R scripts for data generation, model fitting, diagnostics, and plotting.

Current active pipeline files:

```text
code/00_generate_data.R
code/00_generate_montecarlo.R
code/01_fit_bgam.R
code/09_run_montecarlo_sem.R
code/03_compare_methods.R
```

```text
slurm/
```

Cluster submission scripts.

Current hard-threshold SEM experiments:

```text
slurm/gen_montecarlo_n20_hard.sh
slurm/run_montecarlo_n20_hard_sem.sh
slurm/run_montecarlo_n20_hard_sem_k1212_test.sh
slurm/run_montecarlo_p5_hard_sem_test.sh
slurm/run_montecarlo_p10_hard_sem_test.sh
```

```text
EXPERIMENT_STATUS.md
```

Current SSH job IDs, data directories, model settings, and synchronization notes.

## Data And Results

```text
data/
```

Local data files. Large Monte Carlo datasets may exist on SSH only.

```text
results/
```

Model outputs, diagnostics, comparison CSVs, and generated figures.

```text
figures/
pred_vs_raw_compare/
```

Older and exploratory plots kept for reference.

## Documentation

```text
docs/reports/
```

Report-style `.tex` and `.pdf` files moved from the project root.

```text
notes/
```

Working notes, progress reports, and method rationale.

```text
references/
```

Reference books and large PDF references that were previously in the project root.

## Archive

```text
archive/
```

Old code, original imported model files, local backups, LaTeX build artifacts, and temporary PDF-edit images.

The latest cleanup moved build products and old snapshots into:

```text
archive/build_artifacts_20260614_111621
archive/pdf_edit_tmp_20260614_111621
archive/snapshots_20260614
```

Old duplicate `00_generate_data.R` files were renamed inside the snapshot
archive as:

```text
OLD_DO_NOT_USE_00_generate_data.R
```
