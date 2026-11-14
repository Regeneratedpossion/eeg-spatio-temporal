# Thesis Figure 3.6 (fig:ppc): observed value and 90% posterior predictive interval of every converged dataset.
# Uses window mean and peak amplitude of the difference wave. Datasets are ordered by the observed value.
# Usage: Rscript code/08_ppc/fig_ppc.R <results/ppc_checks dir> <output dir>
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly = TRUE)
x <- read.csv(file.path(args[[1]], "ppc_all_datasets.csv"))
x <- x[x$statistic %in% c("mean_160_340_uV", "peak_amplitude_uV"), ]
x$statistic <- factor(x$statistic, c("mean_160_340_uV", "peak_amplitude_uV"), c("Window mean (µV)", "Peak amplitude (µV)"))
x$scenario <- factor(x$scenario, 1:4, paste("Scenario", 1:4))
x$order <- ave(x$observed, x$statistic, x$scenario, FUN = function(v) rank(v, ties.method = "first"))
fig <- ggplot(x, aes(order)) +
  geom_linerange(aes(ymin = rep_q05, ymax = rep_q95), colour = "#7F77DD", alpha = 0.6, linewidth = 0.4) +
  geom_point(aes(y = observed), size = 0.5) +
  facet_grid(statistic ~ scenario, scales = "free", switch = "y") +
  labs(x = "Dataset (ordered by the observed value)", y = NULL) +
  theme_classic() + theme(strip.background = element_blank(), strip.placement = "outside")
ggsave(file.path(args[[2]], "fig_ppc.pdf"), fig, width = 8, height = 4.5)
ggsave(file.path(args[[2]], "fig_ppc.png"), fig, width = 8, height = 4.5, dpi = 110)
