# Thesis Figure D.1 (fig:divergences).
source(here::here("code/00_config.R")); setwd(thesis_dir)   # reads and writes relative to the thesis folder
library(ggplot2); library(dplyr); library(tidyr)
set.seed(1)   # jitter
d <- read.csv("fig_data/div_vs_error.csv") |>
  pivot_longer(c(onset_error_ms, offset_error_ms), names_to = "measure", values_to = "err") |>
  mutate(measure = recode(measure, onset_error_ms = "Onset error (ms)", offset_error_ms = "Offset error (ms)"))
rho <- d |> group_by(scenario, measure) |> summarise(r = cor(divergences, err, method = "spearman", use = "complete.obs"), .groups = "drop") |>
  mutate(lab = sprintf("r = %.2f", r))
p <- ggplot(d, aes(divergences, err)) +
  geom_hline(yintercept = 0, colour = "grey80") +
  geom_jitter(height = 0.4, width = 0, alpha = 0.5, size = 1, colour = "#4575b4") +
  geom_text(data = rho, aes(x = Inf, y = Inf, label = lab), hjust = 1.1, vjust = 1.5, size = 2.8, inherit.aes = FALSE) +
  facet_grid(measure ~ scenario, labeller = labeller(scenario = function(x) paste("Scenario", x)), scales = "free_x") +
  labs(x = "Number of divergent transitions", y = "Error (ms)") +
  theme_minimal(base_size = 9) + theme(panel.grid.minor = element_blank())
ggsave(file.path(fig_dir(thesis_path("diagnostics")), "fig_divergences.pdf"), p, width = 8, height = 4)
