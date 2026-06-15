# TODO

## Immediate (this week)

- [ ] Send `notes/progress_report.md` + Figure 1 PDF to advisor
- [ ] Decide target venue (workshop / conference / journal)
- [ ] Decide method-comparison scope (3 / 5 / 8 methods)

## Short term (1-2 weeks)

- [ ] Implement `code/03_compare_methods.R`
  - [ ] BGAM (already done)
  - [ ] Mass univariate t-test + FDR (BH) — pure R, ~30 lines
  - [ ] Cluster mass — needs MNE-Python via reticulate
  - [ ] (optional) TFCE — same call, different threshold spec
  - [ ] (optional) GAMM 1D × N — heavy, needs cluster resources
- [ ] Implement `code/05_visualise_methods.R`
  - [ ] Side-by-side topomap comparison
  - [ ] onset/offset bias bar chart
- [ ] Run full diagnostics on cluster (`sbatch slurm/run_diagnostics.sh`)

## Medium term (2-4 weeks)

- [ ] Run 100 simulations to estimate metric variance
  - [ ] Adapt `brms_meeg/code/errors_simulation_cluster.R` template
  - [ ] Heavy compute — needs ultralong partition
- [ ] Sensitivity analyses
  - [ ] k = c(8, 8) vs c(10, 12) vs c(15, 15)
  - [ ] sds prior std = 0.5, 1, 2.5
  - [ ] n_participants = 5, 10, 20
- [ ] Compare BGAM (3D joint) vs GAMM (1D × N) — direct test of
  whether 3D pooling beats N independent 1D fits

## Long term (1-3 months)

- [ ] Real M/EEG dataset application
- [ ] Manuscript draft
  - [ ] Background + Nalborczyk-style intro
  - [ ] Methods (model spec + simulation protocol)
  - [ ] Results (this pilot + Monte Carlo + real data)
  - [ ] Discussion (trade-offs: amplitude shrinkage vs pattern recovery)
- [ ] Code release on GitHub
- [ ] Submit

## Possible refinements (only if needed)

- [ ] `bs = "gp"` (Matérn) basis — may reduce spatial diffusion
- [ ] R-INLA + SPDE — for computational scaling
- [ ] GPP (Gaussian predictive processes) — for very large n
- [ ] Posterior odds visualisation per (electrode, time) — paper Figure 4 style

## Bug / clean-up

- [ ] `code/08_slide_topomap_n20.R` produces a contour warning when
      cond1 amplitude is near zero — cosmetic only
- [ ] PNG produced via type="cairo" — verify on cluster install
