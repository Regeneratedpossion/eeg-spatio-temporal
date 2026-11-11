# Code for the thesis

This folder holds all code of the thesis. The scripts are numbered. Run them in the order of the numbers.

## Before you start

Run every script from the project root. The project root is the folder that holds the `.here` file.

```bash
Rscript code/06_tables/06_tables.R passed
```

## The steps

| Step | Folder | What it does | Where it runs |
|---|---|---|---|
| 00 | `00_config.R` | Sets the paths. Every other script loads it. | loaded automatically |
| 01 | `01_simulate/` | Simulates the datasets. | cluster or local |
| 02 | `02_fit/` | Fits the BGAM to one dataset. | cluster, about 20 h per fit |
| 03 | `03_conventional/` | Runs the conventional tests on one dataset. | cluster |
| 04 | `04_evaluate/` | Turns one BGAM fit into detection masks and measures. | cluster |
| 05 | `05_aggregate/` | Collects the results of all datasets. | cluster |
| 06 | `06_tables/` | Makes the tables of the thesis. | local |
| 07 | `07_figures/` | Makes the figures of the thesis. | local |

`R/` holds helper functions. The numbered scripts load them. You do not run them directly.

## The scripts of each step

**Step 01, simulate**

- `01_generate_dataset.R` simulates one dataset. It writes `sim_NNN.rds` (the EEG data) and `gt_NNN.rds` (the true effect) to `data/scenario*/`.
- `01_generate_data.R` holds the generator functions.
- `01_make_adjacency_61ch.py` and `01_adjacency_csv_to_rds.R` rebuild the electrode adjacency file in `resources/adjacency/`. You only need them if the electrode layout changes.

**Step 02, fit**

- `02_fit_bgam.R` fits the BGAM to one dataset. The model settings are in `R/fit_config.R`.

**Step 03, conventional methods**

- `03_evaluate_conventional.R` runs the uncorrected test, BH, BY, Holm, cluster mass and TFCE.
- `03_evaluate_conventional_no_space.R` runs the same tests without spatial adjacency. This is a control analysis.

**Step 04, evaluate**

- `04_evaluate_bgam.R` computes the detection masks and the measures of one BGAM fit.

**Step 05, aggregate**

- `05_aggregate_bgam_evaluation.R` summarises the evaluations of all seeds of one scenario.
- `05_make_passed_fits.R` writes `results/passed_fits.txt`. This file lists the 362 fits that pass the convergence rule (R-hat < 1.01, ESS > 400).

**Step 06, tables**

- `06_tables.R` makes the result tables of Chapter 3.

**Step 07, figures**

- `07_figure_*.R` make the figures.
- `07_figures/data/07_data_*.R` make the data that some figures need.
- `07_figures/README.md` lists which script makes which figure.

## Supporting scripts

`extra/` holds scripts for additional outputs. The main results do not need them. Some figures and tables do.

| Script | What it does |
|---|---|
| `04_extract_fixed_grid.R` | For one seed: saves the fitted values on a fixed grid. |
| `04_save_masks.R` | For one seed: saves the detection masks. |
| `04_per_chain_coverage.R` | For one fit: saves the onset and offset of each chain. |
| `05_aggregate_fixed_grid.R` | Collects the fixed grids of all seeds. |
| `05_agg_sesoi.R` | Compares the two SESOI rules. |
| `05_agg_threshold_sweep.R` | Summarises the performance across thresholds. |
| `05_agg_ring_fpr.R`, `05_agg_electrode_fpr.R` | Summarise the false positive rate per ring and per electrode. |
| `05_agg_no_space.R`, `05_agg_no_space_paired.R` | Summarise the control analysis without spatial adjacency. |

## Running on the cluster

`slurm/` holds the job scripts. Each job script has the same number as the step that it runs.

Two variables tell a job which dataset to use:

- `COND` is the scenario (1, 2, 3 or 4). You set it when you submit.
- `SLURM_ARRAY_TASK_ID` is the seed (1 to 100). Slurm sets it for each task of the array.

Set `R_LIBS_USER` and `CMDSTAN_PATH` in your shell before you submit. Then submit from the project root:

```bash
sbatch --export=ALL,COND=4 slurm/02_fit_array.sh
```

## Optional settings

You can change three folders with environment variables. The defaults work without any change.

| Variable | What it changes | Default |
|---|---|---|
| `PROJECT_ROOT` | the project folder | the folder with the `.here` file |
| `THESIS_DIR` | the folder of the results | `results/` |
| `FIG_OUT` | the output folder of figures and tables | a subfolder of `results/` |

## More information

- `inventory.csv` lists every table and figure with its script, input and output.
- The older code is on the branch `archive`.
