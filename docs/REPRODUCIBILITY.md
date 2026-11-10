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

These checks do not replace running the code with all original fitted models.

## Still to verify

- Compare the updated BGAM evaluation and mask outputs with real `fit.rds` files on LIDO.
- Check missing-value handling in the no-space analysis used for **Table F.1**.
- Verify electrode-level FPR outputs where the required files are available.
- Check older ground-truth metadata for the reported offset difference (**0.342 s versus 0.340 s**).
- The older `sesoi_all.tex` differs from regenerated values for scenarios 2 and 3; the thesis uses the verified `passed` tables.

Do not edit numerical results by hand. Use `inventory.csv` to trace each result to its source code and data.
