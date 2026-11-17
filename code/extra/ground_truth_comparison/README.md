# Ground truth comparison

This folder holds an exploratory figure. The two smooth ground truths in it are alternatives. The thesis does not use them.

`ground_truth_comparison.png` shows three rows at the same seven time points:

1. The thesis ground truth (Gaussian cut at the half-maximum radius, ERP-template waveform).
2. A Gaussian with a smooth taper.
3. A near-flat-top bump.

The red plane joins the four inner electrodes. The blue plane joins the four outer electrodes. Each label gives the mean effect of that plane in µV.

## Run

Run from the project root, the folder that holds the `.here` file:

```bash
Rscript code/extra/ground_truth_comparison/ground_truth_comparison.R
```

The script uses base R only. It reads these files:

- `code/R/superpopulation_truth.R`
- `resources/brms_meeg/code/erp_template.R`
- `resources/electrodes/electrode_layout_61ch.csv`
