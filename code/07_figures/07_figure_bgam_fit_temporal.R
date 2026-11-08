source(here::here("code/00_config.R")); source(here::here("code/R/passed_fits.R"))   # run from the project root
out_dir <- fig_dir(thesis_path("diagnostics")); truth_file <- thesis_path("fig_data/truth_250.csv")
setwd(project_dir)
suppressMessages({library(dplyr); library(ggplot2); library(patchwork)})
mk <- function(c, scen) { ps <- .man$seed[.man$cond == c]
  ce <- do.call(rbind, lapply(ps, function(s) readRDS(sprintf("data/single_node_cond%s_seeds001_100/gt_%03d.rds", c, s))$continuous_effect[, c("channel", "time", "effect")]))
  ce <- aggregate(effect ~ channel + time, ce, mean); ce$scenario <- scen; ce }
truth4 <- rbind(mk("I", 1), mk("H", 2), mk("F", 3), mk("E", 4))
cat("check CZ 250ms:", truth4$effect[truth4$channel == "CZ" & abs(truth4$time - 0.25) < 1e-9 & truth4$scenario %in% c(1, 3)], "\n")
sc <- c(I = 1, H = 2, F = 3, E = 4)
chs <- c("CZ", "C2", "PZ", "FP1")
d <- bind_rows(lapply(names(sc), function(c) bind_rows(lapply(pgs(sprintf("results/single_node_cond%s_ds_cr_seeds001_100/seed_*/fixed_grid/six_channel.csv", c)), read.csv)) |> mutate(scenario = sc[[c]]))) |>
  filter(channel %in% chs) |>
  mutate(A = ifelse(scenario %in% c(1, 2), 0.5, 1.5), omega = ifelse(scenario %in% c(1, 3), "0.26", "0.70"))
s <- d |> group_by(A, omega, channel, time) |> summarise(m = mean(model_mean_diff), lo = quantile(model_mean_diff, 0.025), hi = quantile(model_mean_diff, 0.975), .groups = "drop")
tr <- truth4 |> filter(channel %in% chs) |> mutate(A = ifelse(scenario %in% c(1, 2), 0.5, 1.5)) |> filter(scenario %in% c(1, 3))
s$channel <- factor(s$channel, chs); tr$channel <- factor(tr$channel, chs)
cols <- c("0.26" = "#5E3C99", "0.70" = "#B2ABD2")
boxes <- data.frame(channel = factor("CZ", chs), xmin = c(120, 300), xmax = c(200, 380))
panel <- function(a, strips = TRUE) {
  ggplot() +
    geom_rect(data = boxes, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), fill = "grey70", alpha = 0.25) +
    geom_ribbon(data = filter(s, A == a, omega == "0.70"), aes(time * 1000, ymin = lo, ymax = hi, fill = omega), alpha = 0.40) +
    geom_ribbon(data = filter(s, A == a, omega == "0.26"), aes(time * 1000, ymin = lo, ymax = hi, fill = omega), alpha = 0.40) +
    geom_line(data = filter(s, A == a, omega == "0.70"), aes(time * 1000, m, colour = omega), linewidth = 1.6) +
    geom_line(data = filter(s, A == a, omega == "0.26"), aes(time * 1000, m, colour = omega), linewidth = 0.7) +
    geom_line(data = filter(tr, A == a), aes(time * 1000, effect, linetype = "Truth"), colour = "black", linewidth = 0.6) +
    scale_linetype_manual(values = c(Truth = "dashed"), name = NULL) +
    facet_wrap(~channel, ncol = 1, strip.position = "right") +
    scale_colour_manual(values = cols, name = expression(omega[A])) + scale_fill_manual(values = cols, name = expression(omega[A])) +
    labs(title = bquote(A == .(a)), x = "Time (ms)", y = "Condition difference (\u00b5V)") +
    theme_minimal(base_size = 11) + theme(strip.text.y.right = if (strips) element_text(angle = 0, face = "bold", size = 11) else element_blank(), panel.grid.minor = element_blank())
}
pz <- local({ source(here::here("code/R/bgam_zoom_panel.R"), local = TRUE); pz })
pA <- (panel(0.5, strips = FALSE) | panel(1.5)) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
p <- wrap_elements(pA) / wrap_elements(pz) + plot_layout(heights = c(7, 3.8)) + plot_annotation(tag_levels = "A")
ggsave(file.path(out_dir, "bgam_fit_temporal.png"), p, width = 7.5, height = 10.8, dpi = 200, bg = "white", device = "png", type = "cairo")
