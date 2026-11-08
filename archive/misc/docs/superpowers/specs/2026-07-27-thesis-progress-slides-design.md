# Thesis Progress Slides Design

## Purpose and audience

The presentation is a 15-minute thesis progress update for the thesis supervisor. It should emphasise the current research question, the verified simulation and modelling setup, the results obtained so far, known limitations, and the next steps. The thesis PDF is the primary source of terminology and methodological definitions.

## Presentation structure

The main presentation will contain approximately 13 slides:

1. Title and date
2. Research question and contribution
3. Overall simulation design
4. Continuous spatio-temporal signal and group-level truth
5. Amplitude and noise heterogeneity
6. Noise structure
7. Comparison methods
8. BGAM and posterior decision rule
9. Evaluation metrics
10. Full BGAM recovery result
11. Spatio-temporal detection result
12. Current diagnostics and limitations
13. Current conclusions and next steps

Technical detail, secondary diagnostics, basis sensitivity, convergence information, and supplementary plots will be retained in an appendix.

## Editing workflow

Slides will be reviewed and changed one at a time. Before changing a slide, its proposed scientific content and wording will be confirmed. After modification, the LaTeX source will be compiled and the resulting PDF checked.

## Slide 1: title page

The title page will contain only:

- Title: `Spatio-Temporal Bayesian GAMs for Simulated EEG Effect Recovery`
- Author: `Jiahui Fan`
- Presentation date

It will not contain a subtitle such as “Thesis progress update.”

## Content constraints

- Do not present unfinished analyses as completed results.
- Preserve the distinction between the continuous injected effect and hard labels used for selected evaluation metrics.
- Represent participant-level onset, amplitude, and noise-scale heterogeneity consistently with the thesis.
- Move detail to the appendix rather than overcrowding the main 15-minute presentation.
