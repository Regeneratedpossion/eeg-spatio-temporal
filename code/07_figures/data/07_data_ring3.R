source(here::here("code/00_config.R")); setwd(thesis_dir)   # reads and writes relative to the thesis folder
# Ring-3 check at 250 ms: detection frequency and mean bias (posterior mean minus truth) per electrode.
suppressPackageStartupMessages({library(dplyr); library(ggplot2)})
truth <- read.csv("fig_data/truth_250ms_s1s3.csv", colClasses = c(scenario = "character"))
rd <- function(s) {
  f <- Sys.glob(sprintf("fig_data/topo_%s/seed_*/fixed_grid/topomap_3t.csv", s))
  x <- do.call(rbind, lapply(f, read.csv)); x <- x[abs(x$time - 0.25) < 1e-9, ]
  x$scenario <- s; x
}
d <- rbind(rd("1"), rd("3"))
d$detected <- d$model_p_sig > 10 / 11   # posterior odds above 10
e <- d |> group_by(scenario, channel, xproj, yproj) |>
  summarise(mean_post = mean(model_mean_diff), det_rate = mean(detected), .groups = "drop") |>
  left_join(truth, by = c("scenario", "channel"))
e$bias <- e$mean_post - e$effect
e$scenario_lab <- factor(ifelse(e$scenario == "1", "Scenario 1 (weak)", "Scenario 3 (strong)"),
                         levels = c("Scenario 1 (weak)", "Scenario 3 (strong)"))
write.csv(e, file.path(fig_dir(thesis_path("fig_data")), "ring3_250ms_summary.csv"), row.names = FALSE)

cat("\nMean bias by ring (posterior mean minus truth, 250 ms):\n")
print(as.data.frame(e |> group_by(scenario, ring) |> summarise(n_el = n(), mean_bias = mean(bias), max_bias = max(bias), mean_det = mean(det_rate), .groups = "drop")), digits = 3)
cat("\nRing-3 electrodes with detection rate > 0 at 250 ms:\n")
print(as.data.frame(e |> filter(ring == 3, det_rate > 0) |> arrange(scenario, desc(det_rate)) |>
        select(scenario, channel, xproj, yproj, det_rate, bias)), digits = 3)

lim <- max(abs(e$bias))
p <- ggplot(e, aes(xproj, yproj, fill = bias)) +
  geom_point(shape = 21, size = 4.6, colour = ifelse(e$ring == 3, "black", "grey70"), stroke = ifelse(e$ring == 3, 0.9, 0.3)) +
  geom_point(data = filter(e, ring == 0), shape = 21, size = 1.2, fill = "black", colour = "black") +
  coord_equal() + facet_wrap(~ scenario_lab) +
  scale_fill_gradient2(low = "#0072B2", mid = "white", high = "#D55E00", midpoint = 0, limits = c(-lim, lim),
                       name = "Bias") +
  theme_void(base_size = 9) +
  theme(strip.text = element_text(face = "bold", margin = margin(b = 4)), legend.position = "bottom")
ggsave(file.path(fig_dir(thesis_path("diagnostics")), "bias_map_250ms.pdf"), p, width = 6.6, height = 3.6)
