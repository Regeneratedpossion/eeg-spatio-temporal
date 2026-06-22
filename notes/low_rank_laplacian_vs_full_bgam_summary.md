# Low-rank Laplacian vs full BGAM summary

Date: 2026-06-22

## Question

Can we reduce the spatial dimension first using a low-rank spatial basis, fit Bayesian temporal GAMMs to basis scores, and reconstruct the sensor x time posterior effect map?

## Key comparison to send

### Full spatiotemporal BGAM, k = (12, 12), P = 5

Model:

- Sensor-level full BGAM
- Tensor-product smooth over `xproj`, `yproj`, and `time`
- `k = c(12, 12)`
- `iter = 3000`, `warmup = 1000`

Results:

- `tp` basis:
  - runtime: ~58.2 h
  - divergent transitions: 0
  - max Rhat: 1.0056
  - temporal recovery: pass
  - top spatial channel: pass
  - peak recovery: ~82.4%
  - mean recovery: ~98.4%

- `ps` basis:
  - runtime: ~17.8 h
  - divergent transitions: 0
  - max Rhat: 1.0056
  - temporal recovery: pass
  - top spatial channel: pass
  - peak recovery: ~82.3%
  - mean recovery: ~98.5%

Interpretation:

The full BGAM recovers the signal well, but it is computationally expensive. The `ps` basis is much faster than `tp` with similar recovery.

### Low-rank Laplacian basis, K = 20 with constant/global component, P = 20

Model:

- Construct graph Laplacian basis from sensor adjacency matrix
- Include one constant/global spatial component
- Use K = 20 basis components in total
- Project sensor-level maps onto basis scores
- Fit Bayesian temporal GAMMs to basis scores
- Reconstruct sensor x time posterior effect map
- `iter = 3000`, `warmup = 1000`

Results:

- spatiotemporal correlation: 0.962
- effect-window spatiotemporal correlation: 0.964
- peak spatial correlation: 0.976
- peak recovery: 81.4%
- truth-weighted mean recovery: 104.5%
- top-25% truth-region mean recovery: 102.7%
- temporal recovery: pass
- top spatial channel: pass

Convergence note:

- The K = 20 + constant version still had a small convergence issue due to divergent transitions in high-order basis components.
- However, recovery quality was close to the full BGAM once the constant/global component was included.

Interpretation:

The low-rank Laplacian approach is much faster and becomes much more accurate after adding the constant/global spatial component. It looks defensible as a computational approximation to the full spatiotemporal BGAM, not as a cost-free replacement.

## Suggested email text

Hi Ladislas,

I wanted to ask whether the following could be a reasonable direction.

Instead of fitting the full spatiotemporal BGAM directly on the sensor-level data, I am considering a two-step low-rank spatial basis approach. For example, we could first construct a graph Laplacian basis from the sensor adjacency matrix, include a constant/global spatial component, project the sensor-level data onto the resulting low-rank basis, then fit Bayesian temporal GAMMs to the basis scores, and finally reconstruct the sensor x time posterior effect map.

The motivation is to keep the sensor topology through an adjacency-informed basis, while reducing the spatial dimension before Bayesian fitting.

Some preliminary results are:

Full BGAM, k = (12,12), P = 5:

- tp basis: runtime ~58.2 h, divergent = 0, max Rhat = 1.006, peak recovery ~82.4%, mean recovery ~98.4%
- ps basis: runtime ~17.8 h, divergent = 0, max Rhat = 1.006, peak recovery ~82.3%, mean recovery ~98.5%

Laplacian low-rank basis:

- K = 10 without constant/global component: runtime ~25 min, divergent = 0, max Rhat < 1.006, peak recovery = 44.8%, mean recovery = 52.2%
- K = 20 with constant/global component: spatiotemporal correlation = 0.962, peak recovery = 81.4%, truth-weighted mean recovery = 104.5%, top-25% truth-region mean recovery = 102.7%

So the low-rank approach seems much faster, and the amplitude recovery improves substantially once the constant/global spatial component is included. It may therefore be a useful approximation to the full spatiotemporal BGAM, rather than a replacement without cost.

Do you think this is a defensible direction to pursue?

Best,
Jiahui

## Local code and result files

Full BGAM k = (12,12):

- Code: `code/11_basis_screening.R`
  - `SCREEN_SAVE_FIT=true`
- Slurm: `slurm/run_basis_screening_full_ab.sh`
- Diagnostics:
  - `results/diagnostics/june_full/summary_fit_A_k1212_tp_iter3000_warmup1000.csv`
  - `results/diagnostics/june_full/summary_fit_B_k1212_ps_iter3000_warmup1000.csv`

Laplacian low-rank basis:

- Single workflow code: `code/12_laplacian_basis_brms_fixed.R`
  - `LAPLACE_ACTION=fit`
  - `LAPLACE_ACTION=diagnose_hard`
  - `LAPLACE_ACTION=diagnose_continuous`
- Main note: `notes/laplacian_k20_constant_result.md`
- Slurm:
  - `slurm/run_laplacian_p20_gaussian_hardlabel_k20_constant_long.sh`
  - `slurm/run_laplacian_p20_gaussian_hardlabel_k20_constant_recovery.sh`
  - `slurm/run_laplacian_p20_gaussian_hardlabel_k20_constant_continuous_diag.sh`
- Diagnostics:
  - `results/diagnostics/laplacian_basis_fixed/summary_laplacian_k20_seed001_iter3000_warmup1000_adapt0p999_constant.csv`
  - `results/diagnostics/laplacian_continuous/summary_continuous_laplacian_k20_seed001_iter3000_warmup1000_adapt0p999_constant.csv`
