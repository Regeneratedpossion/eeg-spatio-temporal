# Thesis Figure G.1 (fig:ppc_sd): observed residual SD of condition 2 and its 90% posterior predictive interval in every dataset.
# Datasets are ordered by the observed value. Each panel has its own axes.
# Usage: Rscript code/08_ppc/fig_ppc_sd.R <results/ppc_checks dir> <output dir>
suppressPackageStartupMessages(library(ggplot2))
args <- commandArgs(trailingOnly = TRUE)
x <- read.csv(file.path(args[[1]], "ppc_all_datasets.csv"))
x <- x[x$statistic == "cond2_sd", ]
x$scenario <- factor(x$scenario, 1:4, paste("Scenario", 1:4))
x$order <- ave(x$observed, x$scenario, FUN = function(v) rank(v, ties.method = "first"))
fig <- ggplot(x, aes(order)) +
  geom_linerange(aes(ymin = rep_q05, ymax = rep_q95), colour = "#7F77DD", linewidth = 0.5) +
  geom_point(aes(y = observed), size = 0.25) +
  facet_wrap(~scenario, nrow = 1, scales = "free") +
  labs(x = "Dataset (ordered by the observed value)", y = "Residual SD, condition 2 (\u00b5V)") +
  theme_classic() + theme(strip.background = element_blank())
ggsave(file.path(args[[2]], "fig_ppc_sd.pdf"), fig, width = 8, height = 2.6)
ggsave(file.path(args[[2]], "fig_ppc_sd.png"), fig, width = 8, height = 2.6, dpi = 110)
