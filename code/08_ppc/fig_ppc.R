# Thesis Figure 3.6 (fig:ppc): posterior predictive check of the difference waveform.
# Scenario 1 seed 1 and scenario 3 seed 36. Needs the draw files made by code/08_ppc/ppc_diff_draws.R.
# Usage: Rscript code/08_ppc/fig_ppc.R <dir with ppc_obs_*.csv and ppc_rep_*.csv> <output dir>
suppressPackageStartupMessages({library(ggplot2); library(patchwork)})
args <- commandArgs(trailingOnly = TRUE)
S <- args[[1]]; out <- args[[2]]
win <- c(160, 340)
load_w <- function(tag) {
  o <- read.csv(file.path(S, sprintf("ppc_obs_%s.csv", tag))); r <- read.csv(file.path(S, sprintf("ppc_rep_%s.csv", tag)))
  o$time <- o$time * 1000; r$time <- r$time * 1000; list(obs = o, rep = r)
}
panel_wave <- function(w, title) ggplot() +
  annotate("rect", xmin = win[1], xmax = win[2], ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.5) +
  geom_line(data = w$rep, aes(time, value, group = draw), colour = "#7F77DD", alpha = 0.08, linewidth = 0.3) +
  geom_line(data = w$obs, aes(time, obs), colour = "black", linewidth = 0.9) +
  labs(title = title, x = "Time (ms)", y = expression("Difference ("*mu*"V)")) + theme_classic()
panel_stat <- function(w) {
  in_w <- w$rep$time >= win[1] & w$rep$time <= win[2]
  t_rep <- tapply(w$rep$value[in_w], w$rep$draw[in_w], mean)
  t_obs <- mean(w$obs$obs[w$obs$time >= win[1] & w$obs$time <= win[2]])
  ggplot(data.frame(t = as.numeric(t_rep)), aes(t)) +
    geom_histogram(bins = 30, fill = "#7F77DD", colour = "white", alpha = 0.7) +
    geom_vline(xintercept = quantile(t_rep, c(0.05, 0.95)), linetype = "dashed", colour = "grey40", linewidth = 0.5) +
    geom_vline(xintercept = t_obs, colour = "black", linewidth = 0.9) +
    labs(x = expression(T[mean]~"("*mu*"V)"), y = "Count") + theme_classic()
}
w1 <- load_w("scenario1_seed001"); w3 <- load_w("scenario3_seed036")
fig <- (panel_wave(w1, "Weak effect") | panel_wave(w3, "Strong effect")) / (panel_stat(w1) | panel_stat(w3)) + plot_annotation(tag_levels = "A")
ggsave(file.path(out, "fig_ppc.pdf"), fig, width = 8, height = 6)
ggsave(file.path(out, "fig_ppc.png"), fig, width = 8, height = 6, dpi = 110)
