source(here::here("code/00_config.R")); setwd(thesis_dir)   # reads and writes relative to the thesis folder
# Figure: BGAM analysis of ONE simulated dataset (scenario 3, seed 36)
library(ggplot2); library(dplyr); library(tidyr); library(patchwork)
cells <- read.csv("fig_data/cells_scenario3_seed036.csv")
thr <- 10 / 11
seed_label <- "Scenario 3, seed 36 (converged dataset whose F1 equals the scenario median)"

th <- theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank(),
                                           plot.title = element_text(face = "bold"))
win <- annotate("rect", xmin = 160, xmax = 340, ymin = -Inf, ymax = Inf, fill = "grey90")
two <- filter(cells, electrode %in% c("CZ", "PZ")) |>
  mutate(electrode = factor(electrode, levels = c("CZ", "PZ")))

pA <- ggplot(two, aes(time * 1000)) + win +
  geom_ribbon(aes(ymin = g_lo, ymax = g_hi), fill = "#2166AC", alpha = 0.25) +
  geom_line(aes(y = g_mean), colour = "#2166AC", linewidth = 0.8) +
  geom_line(aes(y = truth), linetype = "dashed") +
  facet_wrap(~electrode, ncol = 1, scales = "free_y", strip.position = "left") +
  labs(title = "A  Estimate", x = "Time (ms)", y = expression("Condition difference ("*mu*V*")")) + th +
  theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, face = "bold"))

pB <- ggplot(two, aes(time * 1000, prob)) + win +
  geom_hline(yintercept = thr, linetype = "dashed", colour = "grey40") +
  geom_line(colour = "#2166AC", linewidth = 0.8) +
  geom_point(data = filter(two, det_bgam), colour = "#1B7837", size = 1.2) +
  facet_wrap(~electrode, ncol = 1) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  labs(title = expression(bold("B  Posterior probability ")*P(g>delta)),
       x = "Time (ms)", y = "Probability") + th + theme(strip.text = element_blank())

rk <- cells |> distinct(electrode, d) |> arrange(d, electrode) |> mutate(row = row_number())
grid <- cells |> left_join(rk[, c("electrode", "row")], by = "electrode") |>
  pivot_longer(c(det_bgam, det_cm), names_to = "method", values_to = "det") |>
  mutate(method = recode(method, det_bgam = "BGAM", det_cm = "Cluster mass"))
n_act <- n_distinct(cells$electrode[cells$active])
lab <- grid |> distinct(electrode, row) |> filter(electrode %in% c("CZ", "PZ", "FP1"))
pC <- ggplot(grid, aes(time * 1000, row)) +
  geom_tile(aes(fill = det), width = 10, height = 1) +
  annotate("rect", xmin = 155, xmax = 345, ymin = 0.5, ymax = n_act + 0.5,
           fill = NA, colour = "black", linewidth = 0.5) +
  scale_fill_manual(values = c(`FALSE` = "grey95", `TRUE` = "#1B7837"),
                    labels = c("not detected", "detected"), name = NULL) +
  scale_y_reverse(expand = c(0, 0), breaks = lab$row, labels = lab$electrode) + scale_x_continuous(expand = c(0, 0)) +
  facet_wrap(~method, ncol = 1, strip.position = "top") +
  labs(title = "C  Detections", x = "Time (ms)", y = "Electrode") +
  th + theme(panel.grid = element_blank(), legend.position = "bottom", strip.text = element_text(face = "bold", hjust = 0))

fig <- (pA | pB | pC) + plot_layout(widths = c(1, 1, 1.2))
ggsave(file.path(fig_dir(thesis_path("diagnostics")), "fig_one_dataset.pdf"), fig, width = 10, height = 5)
ggsave(file.path(fig_dir(thesis_path("diagnostics")), "fig_one_dataset.png"), fig, width = 10, height = 5, dpi = 150)
