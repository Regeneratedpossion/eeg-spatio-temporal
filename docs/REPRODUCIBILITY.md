# Reproducibility

## What can be reproduced?

The repository includes selected evaluation results needed for thesis tables and figures. Large fitted models (`fit.rds`) and some cluster-only outputs are not included.

- **Tables:** run `Rscript code/06_tables/06_tables.R passed` from the project root.
- **Figures:** use the scripts in `code/07_figures/` and check `inventory.csv` for required inputs.
- **Full simulation:** run steps 01–05 using the scripts in `slurm/` on LIDO or a comparable environment.

## Checks completed

- The passed-fit list matches the evaluated datasets: **362 of 400 fits**.
- `06_tables.R` produces identical outputs before and after cleanup in both `all` and `passed` modes.
- Small-data tests found no change in the cleaned simulation, conventional-evaluation and BGAM-evaluation calculations.
- On LIDO, scenario 1 with seed 1 and its real `fit.rds`: the cleaned BGAM evaluation, conventional evaluation, detection masks, per-chain coverage and fixed-grid scripts give the same outputs as the stored cluster results.
- On LIDO, the BGAM evaluation was run again on one fit of each scenario (scenario 1 seed 1, scenario 2 seed 2, scenario 3 seed 1, scenario 4 seed 1). `evaluation_main.csv` is identical to the stored file in all four. The job script `slurm/04_evaluate_bgam_array.sh` gives the same file for scenario 1, seed 1.
- The no-space summary (**Table F.1**) and the electrode-level FPR summary, recomputed from the stored results on LIDO, equal the stored files. The no-space inputs have no missing values.
- Older ground-truth files store `true_offset` = 0.342 s, the code uses 0.340 s. The evaluation reads 0.340 s from `evaluation_settings()`, and the 0.01 s time grid has no point between the two values, so the results do not change.
- The older `sesoi_all.tex` had scenarios 2 and 3 swapped. It is regenerated; all other tables are identical to the regenerated ones.
- All scripts in `code/07_figures/` run on the stored results.

These checks do not replace running the code with all original fitted models. The conventional, mask, per-chain coverage and fixed-grid comparisons cover scenario 1, seed 1 only. The BGAM evaluation comparison covers one fit per scenario.

## Model checks

`code/08_ppc/` checks the BGAM fits with posterior predictive checks and a residual diagnostic. Both run on the 362 passed fits (80, 86, 97 and 99 fits in scenarios 1 to 4). The 38 fits that did not pass the convergence rule are not checked. They are listed in `results/ppc_checks/not_passed_fits.csv` with the scenario, the seed and the failed criterion.

- The checks read the saved `fit.rds`. They do not refit the model. They need the cluster: `slurm/08_ppc_array.sh` and `slurm/08_resid_source_array.sh`.
- The posterior predictive checks use `posterior_predict()` with 200 posterior draws and the participant effects of the observed participants. The residual diagnostic uses `posterior_epred()` only.
- The tables in `results/ppc_checks/` give the median over datasets. The numbers in brackets are the 5% and 95% quantiles over datasets.
- Four datasets (one per scenario) were run again with the final scripts. The output is identical to the first run.
- The checks describe the residual structure of the fits. They do not show that the posterior odds are calibrated.

## Known differences

- Table F.1 (no-space control): two cells differ by 0.001 from the values recomputed from the stored results. Scenario 1, T: −0.01246 recomputed, −0.013 in the table. Scenario 4, T: −0.00649 recomputed, −0.007 in the table. The recomputed values equal the stored file. The table in the thesis is not changed by this repository.

Do not edit numerical results by hand. Use `inventory.csv` to trace each result to its source code and data.
