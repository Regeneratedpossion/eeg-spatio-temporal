source(here::here("code/00_config.R")); setwd(thesis_dir)   # reads and writes relative to the thesis folder
# Why the participant amplitude A_p is log-normal rather than normal.
library(ggplot2); library(dplyr)
omegas <- c(0.26, 0.70)
d <- expand.grid(x = seq(-1.5, 3.5, by = 0.005), omega = omegas) |>
  mutate(sA = sqrt(exp(omega^2) - 1),
         Normal = dnorm(x, 1, sA),
         `Log-normal` = ifelse(x > 0, dlnorm(x, -omega^2 / 2, omega), 0)) |>
  tidyr::pivot_longer(c(Normal, `Log-normal`), names_to = "dist", values_to = "dens") |>
  mutate(panel = sprintf("omega[A] == %.2f", omega))
neg <- d |> filter(dist == "Normal", x < 0)
lab <- d |> distinct(omega, sA, panel) |>
  mutate(p = pnorm(0, 1, sA), txt = sprintf("P(A[p] < 0) == '%.1f%%'", 100 * p))
p <- ggplot(d, aes(x, dens, colour = dist, linetype = dist)) +
  geom_area(data = neg, aes(x, dens), fill = "#B2182B", alpha = 0.3,
            colour = NA, inherit.aes = FALSE) +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_line(linewidth = 0.8) +
  geom_text(data = lab, aes(x = -1.4, y = 1.45, label = txt), parse = TRUE,
            hjust = 0, size = 3, colour = "#B2182B", inherit.aes = FALSE) +
  facet_wrap(~panel, labeller = label_parsed) +
  scale_colour_manual(values = c("Normal" = "grey40", "Log-normal" = "#2166AC"), name = NULL) +
  scale_linetype_manual(values = c("Normal" = "dashed", "Log-normal" = "solid"), name = NULL) +
  labs(x = expression(A[p] / A), y = "Density") +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(file.path(fig_dir(thesis_path("diagnostics")), "fig_lognormal.pdf"), p, width = 6.5, height = 3)
ggsave(file.path(fig_dir(thesis_path("diagnostics")), "fig_lognormal.png"), p, width = 6.5, height = 3, dpi = 150)
print(lab[, c("omega", "sA", "p")])
