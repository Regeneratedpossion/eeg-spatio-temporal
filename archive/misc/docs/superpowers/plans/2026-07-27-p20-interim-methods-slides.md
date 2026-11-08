# P=20 Interim Methods Slides Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the obsolete P=5 comparison block with four internally consistent P=20 interim frames based on 18 common Monte Carlo seeds.

**Architecture:** Treat the local per-seed CSV files as the numeric source of truth, derive one reproducible summary table, and then replace only the four relevant Beamer frames. Verification covers both numerical agreement and rendered-page legibility.

**Tech Stack:** R CSV summarization, LaTeX Beamer, `latexmk`, Poppler PDF tools.

## Global Constraints

- Use P=20 and the same 18 common seeds in low, medium, and high heterogeneity.
- Use alpha=0.05 and 1,024 permutations for cluster mass and TFCE.
- Label the unadjusted method `Raw p-value`.
- Report timing errors in milliseconds and identify valid timing denominators.
- Mark every result frame as interim 18/30 MC.
- Preserve unrelated changes in `slides/thesis_progress_slides.tex`.

---

### Task 1: Derive the slide summaries

**Files:**
- Read: `/private/tmp/conventional_n30_interim/sim_*.csv`
- Create: `/private/tmp/conventional_n30_interim/slide_summary_common18.csv`

**Interfaces:**
- Consumes: per-seed columns `scenario`, `seed`, `method`, timing errors, `tp`, `fp`, `fn`, `tn`, `f1`, and `iou`.
- Produces: one row per scenario and method with timing Bias/MAE/RMSE/SD/count plus mean FPR, FNR, F1, and IoU.

- [ ] **Step 1: Verify the common-seed invariant**

Run an R command that reads every `sim_*.csv`, intersects seed sets across all three scenarios, and asserts that the intersection has length 18.

- [ ] **Step 2: Compute the exact summary**

For onset and offset error vectors, compute `mean(x)`, `mean(abs(x))`, `sqrt(mean(x^2))`, `sd(x)`, and `sum(!is.na(x))`. Compute FPR as `fp/(fp+tn)`, FNR as `fn/(fn+tp)`, and take Monte Carlo means of FPR, FNR, F1, and IoU.

- [ ] **Step 3: Cross-check known values**

Assert low-Holm onset bias is `-7.2222` ms within tolerance, medium-Holm mean IoU is `0.686`, and high-BY valid timing count is 16.

### Task 2: Replace the four Beamer frames

**Files:**
- Modify: `slides/thesis_progress_slides.tex` in the block beginning `Methods comparison ($P=5$ pilot, 30 MC)` and ending before `\appendix`.

**Interfaces:**
- Consumes: `/private/tmp/conventional_n30_interim/slide_summary_common18.csv`.
- Produces: two comparison matrices, one definitions frame, and one FP/FN diagnosis matrix.

- [ ] **Step 1: Replace the comparison frame**

Use `apply_patch` to create separate onset and offset frames. Each contains low, medium, and high blocks with Raw p-value, BH, BY, Holm, Cluster mass, and TFCE in fixed order. The onset columns are Bias, MAE, RMSE, SD, and valid n; the offset frame adds mean F1 and IoU.

- [ ] **Step 2: Update the definitions frame**

Set alpha to 0.05 and display FPR, FNR, F1, and IoU formulas for the complete time-by-electrode mask.

- [ ] **Step 3: Replace the diagnosis frame**

Use the 18-row scenario-by-method matrix with mean FPR, FNR, F1, and IoU. Add a short interim interpretation: Raw p-value, cluster mass, and TFCE have high recall but broad masks; Holm performs best in low/medium, whereas BY is more balanced in high heterogeneity.

- [ ] **Step 4: Check forbidden stale copy**

Search the edited frame block for `P=5`, `0.10`, `Uncorrected`, and bare method label `Raw $p$`; expected result: no matches.

### Task 3: Compile and inspect the PDF

**Files:**
- Regenerate: `slides/thesis_progress_slides.pdf`
- Inspect: `slides/thesis_progress_slides.log`

**Interfaces:**
- Consumes: updated Beamer source.
- Produces: compiled, visually checked presentation PDF.

- [ ] **Step 1: Compile in halt-on-error mode**

Run `LC_ALL=C LANG=C latexmk -pdf -interaction=nonstopmode -halt-on-error thesis_progress_slides.tex` from `slides/`; expected exit status 0.

- [ ] **Step 2: Check structural warnings and PDF text**

Search the log for undefined references and overfull boxes, then extract the four frame titles with `pdftotext -layout` and confirm current labels and values appear.

- [ ] **Step 3: Render and inspect all four frames**

Locate the four page numbers from extracted text, render them at 150 dpi with `pdftoppm`, and visually inspect for clipping, collisions, and illegible type. If needed, adjust table spacing and repeat compilation.

- [ ] **Step 4: Review the final diff**

Run `git diff -- slides/thesis_progress_slides.tex` and confirm no content outside the intended frame block changed.
