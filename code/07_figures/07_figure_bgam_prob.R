# Thesis Figure 3.3 (fig:bgam_prob).
source(here::here("code/00_config.R")); source(here::here("code/R/passed_fits.R"))   # run from the project root
out_dir <- fig_dir(thesis_path("diagnostics"))
setwd(project_dir)
suppressMessages({library(dplyr); library(ggplot2); library(patchwork)})
sc <- c("1" = 1, "2" = 2, "3" = 3, "4" = 4)
chs <- c("CZ", "C2", "PZ", "FP1")
d <- bind_rows(lapply(names(sc), function(c) bind_rows(lapply(pgs(sprintf("results/scenario%s/seed_*/fixed_grid/six_channel.csv", c)), read.csv)) |> mutate(scenario = sc[[c]]))) |>
  filter(channel %in% chs) |>
  mutate(A = ifelse(scenario %in% c(1, 2), 0.5, 1.5), omega = ifelse(scenario %in% c(1, 3), "0.26", "0.70"))
s <- d |> group_by(A, channel, time) |> summarise(m = median(model_p_sig), lo = quantile(model_p_sig, 0.1), hi = quantile(model_p_sig, 0.9), .groups = "drop")
s$channel <- factor(s$channel, chs)
thr <- 10 / 11
panel <- function(a, strips = TRUE) {
  ggplot(filter(s, A == a), aes(time * 1000)) +
    geom_vline(xintercept = c(160, 340), linetype = "dashed", colour = "grey60", linewidth = 0.3) +
    geom_hline(yintercept = thr, linetype = "dashed", colour = "grey30", linewidth = 0.4) +
    geom_text(data = data.frame(channel = factor("CZ", chs), time = 0.5, m = thr - 0.04)[a == 1.5, ], aes(y = m), label = "odds = 10", hjust = 1, vjust = 1, size = 2.6, colour = "grey30") +
    geom_ribbon(aes(ymin = lo, ymax = hi), fill = "#4575b4", alpha = 0.3) +
    geom_line(aes(y = m), colour = "#4575b4", linewidth = 0.8) +
    facet_wrap(~channel, ncol = 1, strip.position = "right") +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    labs(title = bquote(A == .(a)), x = "Time (ms)", y = "Posterior probability") +
    theme_minimal(base_size = 9) + theme(strip.text.y.right = if (strips) element_text(angle = 0, face = "bold", size = 10) else element_blank(), panel.grid.minor = element_blank())
}
p <- panel(0.5, strips = FALSE) | panel(1.5)
ggsave(file.path(out_dir, "bgam_prob_combined.png"), p, width = 7.5, height = 7, dpi = 200, bg = "white")
print(s |> filter(channel %in% c("PZ","FP1"), A == 1.5, time %in% c(0.19,0.25,0.32)) |> as.data.frame())
print(s |> filter(channel == "PZ", A == 1.5, m > thr) |> summarise(from = min(time), to = max(time)))
