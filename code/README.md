# Code for the thesis

Run every script from the project root, the folder that holds the `.here` file.
Paths are found with `here::here()`; `PROJECT_ROOT` overrides the root, `THESIS_DIR` the folder that holds the
small input files and the outputs of tables and figures (default: `results/` in the project) and `FIG_OUT` the output folder of figures and tables.

| Step | Script | What it does |
|---|---|---|
| 00 | `00_config.R` | paths and the scenario map (1 = I, 2 = H, 3 = F, 4 = E) |
| 01 | `01_simulate/01_generate_dataset.R` (+ `01_generate_data.R`, the generator functions) | simulates one dataset; condition from `COND`, seed from `SLURM_ARRAY_TASK_ID` or the first argument; writes `data/single_node_cond*_seeds001_100/sim_NNN.rds` and `gt_NNN.rds` |
| 01 | `01_simulate/01_make_adjacency_61ch.py`, `01_simulate/01_adjacency_csv_to_rds.R` | rebuild `resources/adjacency/adjacency_xproj_yproj_61ch.rds` |
| 02 | `02_fit/02_fit_bgam.R` (+ `R/fit_config.R`) | fits the BGAM of one seed; condition from `COND`, seed from `SLURM_ARRAY_TASK_ID` |
| 03 | `03_conventional/03_evaluate_conventional.R`, `03_conventional/03_evaluate_conventional_no_space.R` | uncorrected, BH, BY, Holm, cluster mass, TFCE |
| 04 | `04_evaluate/04_evaluate_bgam.R` | detection masks and measures of the BGAM (also holds the shared truth and metric functions) |
| 05 | `05_aggregate/05_*.R` | per-seed outputs and aggregation (fixed grid, masks, per-chain onsets, thresholds, SESOI, rings, no-space control) |
| 05 | `05_aggregate/05_make_passed_fits.R` | writes `resources/passed_fits.txt` from the per-seed `evaluation_main.csv` (R-hat < 1.01, ESS > 400) |
| 06 | `06_tables/06_tables.R` | result tables of Chapter 3 |
| 07 | `07_figures/07_figure_*.R`, `07_figures/07_data_*.R`, `07_figures/07_ring3_analysis.R` | figures of the thesis and their data |

`slurm/` has one submit script per step, numbered like the code. Submit them from the project root,
for example `sbatch --export=ALL,COND=E --job-name=ds_cr_condE_fit --output=logs/ds_cr_condE_fit_%A_%a.log slurm/02_fit_array.sh`.
Set `R_LIBS_USER` and `CMDSTAN_PATH` in your shell before submitting; the scripts no longer set them.

`resources/passed_fits.txt` lists the 362 fits that passed the convergence criterion (condition, seed).
The older code is on the branch `archive`. The full list of tables, figures, scripts and checks
is in `inventory.csv`.
