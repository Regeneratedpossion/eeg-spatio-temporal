# P=20 Interim Methods Slides Design

## Scope

Replace the obsolete `P=5` methods-comparison frame and synchronize the two
following FP/FN frames. The three-frame block must describe only the current
`P=20`, `alpha=0.05` conventional-method analysis over the 18 seeds completed
for every heterogeneity scenario. The old pilot results must not remain in this
block.

## Data and naming

- Use the same 18 common seeds for low, medium, and high heterogeneity.
- Label the unadjusted method **Raw p-value**, never “Uncorrected” or “Raw p”.
- Include Raw p-value, BH, BY, Holm, cluster mass, and TFCE.
- State that cluster mass and TFCE use 1,024 permutations.
- Mark the results as interim (`18/30 MC`) and avoid presenting them as final.
- Report timing errors in milliseconds. Where a method has no valid interval,
  state the valid timing denominator (high BY: 16/18; high Holm: 15/18).

## Frame 1: Methods comparison

Use one readable summary frame rather than reproducing every cell of the full
thesis tables. Show mean F1 and IoU for all six methods in each heterogeneity
scenario. Add a compact timing summary for the strongest method in each
scenario, including onset and offset bias, MAE, RMSE, and SD. Identify Holm for
low and medium heterogeneity and BY for high heterogeneity. Footnotes define the
interim sample and permutation count.

## Frame 2: FP and FN definitions

Keep the cellwise confusion-matrix definitions. Change the detection threshold
to `alpha=0.05`. Define FPR as `FP/(FP+TN)`, FNR as `FN/(FN+TP)`, F1 from
precision and recall, and IoU as `TP/(TP+FP+FN)`. Clarify that these metrics use
the complete time-by-electrode mask and the current simulation truth mask.

## Frame 3: FP/FN diagnosis

Replace the old `alpha=0.10`, `P=5` diagnostic table with summaries computed
from the same 18 common seeds as Frame 1. Organize the table by heterogeneity
scenario and method, reporting mean FPR, FNR, F1, and IoU. The accompanying
interpretation should explain the observed precision-recall trade-off without
claiming that the interim ranking is final.

## Verification

- Compile `thesis_progress_slides.tex` with `latexmk` in halt-on-error mode.
- Check the log for undefined references and new overfull boxes.
- Extract PDF text to confirm that `P=5`, `alpha=0.10`, and the old method label
  are absent from the three-frame block.
- Render and visually inspect all three updated frames for clipping and legible
  type.
