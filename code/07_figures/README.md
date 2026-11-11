# Figures of the thesis

Run every script from the repository root. Figures go to `results/diagnostics/`
unless the table says otherwise (`FIG_OUT` overrides the folder).
The figure numbers are those of `thesis_script.tex`. The LaTeX label is the label used in the thesis; `inventory.csv` lists the inputs
of each figure and whether it was regenerated from a clean clone.

| Thesis figure | LaTeX label | Script | Output |
|---|---|---|---|
| Figure 2.2 | `fig:data_generation` | `07_figure_data_generation.R` | `fig_data_generation.png` |
| Figure 2.3 | `fig:electrode_adjacency` | `07_figure_electrode_adjacency.R` | `electrode_adjacency.pdf` |
| Figure 3.1 | `fig:bgam_fit` | `07_figure_bgam_fit_temporal.R` | `bgam_fit_temporal.png` |
| Figure 3.2 | `fig:bgam_spatial` | `07_figure_bgam_fit_spatial.R` | `bgam_fit_spatial.png` |
| Figure 3.3 | `fig:bgam_prob` | `07_figure_bgam_prob.R` | `bgam_prob_combined.png` |
| Figure 3.4 | `fig:detection_counts` | `07_figure_detection_combined.R` | `detection_combined.png`, `detection_counts.csv` |
| Figure 3.5 | `fig:tradeoff` | `07_figure_main_comparison.R` | `fig_tradeoff.pdf` |
| Figures C.5, C.6 | `fig:main_comparison`, `fig:main_comparison_appendix` | `07_figure_main_comparison.R` | `fig_main_comparison.pdf`, `fig_main_comparison_appendix.pdf` |
| Figure B.1 | `fig:lognormal` | `07_figure_lognormal.R` | `fig_lognormal.pdf`, `fig_lognormal.png` |
| Figures C.1 to C.4 | `fig:detection_topo_1..4` | `07_figure_detection_topo.R`, `07_figure_detection_topo_crop.py` | `results/heatmaps/topo_{1,2,3,4}.png` |
| Figure C.7 | `fig:bgam_spatial_scenarios` | `07_figure_bgam_spatial_scenarios.R` | `bgam_fit_spatial_scenarios.png` |
| Figure D.1 | `fig:divergences` | `07_figure_divergences.R` | `fig_divergences.pdf` |
| Figure E.1 | `fig:one_dataset` | `07_figure_one_dataset.R` | `fig_one_dataset.pdf`, `fig_one_dataset.png` |

`07_figure_main_comparison.R` makes several figures (and `figure_summary.csv`, `nvalid_rows.tex`).

## Data scripts

These write small input files for the figures to `results/fig_data/`.

| Script | Output | Used by |
|---|---|---|
| `07_data_truth_250.R` | `truth_250.csv` | `07_figure_bgam_fit_spatial.R`, `07_figure_bgam_fit_temporal.R`, `07_figure_bgam_spatial_scenarios.R` |
| `07_data_truth_250ms_s1s3.R` | `truth_250ms_s1s3.csv` | `07_ring3_analysis.R` |
| `07_data_div_vs_error.R` | `div_vs_error.csv` | `07_figure_divergences.R` |

`07_ring3_analysis.R` writes `results/fig_data/ring3_250ms_summary.csv` and `bias_map_250ms.pdf`.
