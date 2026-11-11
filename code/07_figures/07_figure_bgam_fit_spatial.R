source(here::here("code/00_config.R")); source(here::here("code/R/passed_fits.R"))   # run from the project root
out_dir <- fig_dir(thesis_path("diagnostics")); truth_file <- thesis_path("fig_data/truth_250.csv")
setwd(project_dir)
suppressMessages({library(dplyr); library(ggplot2); library(patchwork)})
source("code/R/utils.R")
sc <- c("1" = 1, "2" = 2, "3" = 3, "4" = 4)
est <- bind_rows(lapply(names(sc), function(c) bind_rows(lapply(pgs(sprintf("results/scenario%s/seed_*/fixed_grid/topomap_3t.csv", c)), read.csv)) |> mutate(scenario = sc[[c]]))) |>
  filter(abs(time - 0.25) < 1e-9) |>
  group_by(scenario, channel, xproj, yproj) |> summarise(est = mean(model_mean_diff), .groups = "drop")
d <- left_join(est, read.csv(truth_file)[, c("scenario", "channel", "effect")], by = c("scenario", "channel"))
peak <- d |> group_by(scenario) |> summarise(pk = max(effect)); d <- left_join(d, peak, by = "scenario")
d <- d |> mutate(Truth = effect / pk, Estimate = est / pk, Diff = (est - effect) / pk)
active <- readRDS("data/scenario1/gt_001.rds")$active_electrodes
d$active <- d$channel %in% active
sens <- distinct(d, channel, xproj, yproj)
hr <- max(sqrt(sens$xproj^2 + sens$yproj^2)) * 1.05
dec <- scalp_decorations(hr)
grid <- expand.grid(xproj = seq(-hr, hr, length.out = 110), yproj = seq(-hr, hr, length.out = 110)) |> filter(sqrt(xproj^2 + yproj^2) <= hr)
d$A <- ifelse(d$scenario %in% c(1, 2), 0.5, 1.5)
act <- d |> filter(active) |> distinct(channel, xproj, yproj)
dm <- d |> group_by(A, channel, xproj, yproj) |> summarise(Estimate = mean(Estimate), Diff = mean(Diff), .groups = "drop")
cz <- sens[sens$channel == "CZ", ]
sigma_space <- 0.5; r_act <- sigma_space * sqrt(2 * log(2))
dm$dist <- sqrt((dm$xproj - cz$xproj)^2 + (dm$yproj - cz$yproj)^2)
dm$effect <- factor(dm$A, c(0.5, 1.5), c("A = 0.5", "A = 1.5"))
tr <- tibble::tibble(dist = seq(0, max(dm$dist), length.out = 500)) |> mutate(rel = exp(-dist^2 / (2 * sigma_space^2)) * (dist <= r_act))
dA <- c(0.045 / 0.486, 0.051 / 1.457)
cols <- c("A = 0.5" = "#D55E00", "A = 1.5" = "#0072B2")
pp <- ggplot() + geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = r_act, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_hline(yintercept = dA, colour = cols, linetype = "dotted", linewidth = 0.9) +
  annotate("text", x = 1.78, y = dA[1], label = "\u03b4 (A = 0.5)", colour = cols[1], hjust = 0, vjust = -0.2, size = 3.6) +
  annotate("text", x = 1.78, y = dA[2], label = "\u03b4 (A = 1.5)", colour = cols[2], hjust = 0, vjust = 1.2, size = 3.6) + scale_x_continuous(limits = c(-0.05, 2.25), breaks = seq(0, 1.5, 0.5)) +
  geom_line(data = tr, aes(dist, rel, linetype = "Truth"), linewidth = 0.6) +
  geom_point(data = dm, aes(dist, Estimate, colour = effect, shape = effect), size = 2.4, stroke = 0.7, alpha = 0.9, position = position_dodge(width = 0.04)) +
  scale_colour_manual(values = cols, name = "BGAM estimate") + scale_shape_manual(values = c(1, 16), name = "BGAM estimate") + scale_linetype_manual(values = "solid", name = NULL) +
  labs(x = "Distance from CZ", y = "Effect at 250 ms\n(fraction of the true peak)") +
  theme_minimal(base_size = 11) + theme(legend.position = "right")
dl <- ceiling(max(abs(dm$Diff)) * 20) / 20
hd <- function(a) {
  e <- filter(dm, A == a); g <- grid; g$v <- idw_interp(g, e[, c("xproj", "yproj")], e$Diff, power = 2, k = 6)
  ggplot(g, aes(xproj, yproj)) + geom_raster(aes(fill = v), interpolate = TRUE) +
    geom_path(data = dec$head_circle, aes(x, y), inherit.aes = FALSE, linewidth = 0.3) +
    geom_path(data = dec$nose, aes(x, y), inherit.aes = FALSE, linewidth = 0.3) +
    geom_path(data = dec$ear_left, aes(x, y), inherit.aes = FALSE, linewidth = 0.3) +
    geom_path(data = dec$ear_right, aes(x, y), inherit.aes = FALSE, linewidth = 0.3) +
    geom_point(data = sens, aes(xproj, yproj), inherit.aes = FALSE, size = 0.25) +
    geom_point(data = act, aes(xproj, yproj), inherit.aes = FALSE, shape = 1, size = 1.8, stroke = 0.4) +
    coord_equal(xlim = c(-hr * 1.12, hr * 1.12), ylim = c(-hr * 1.12, hr * 1.18), expand = FALSE) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, limits = c(-dl, dl),
                         name = "Estimate \u2212 truth (fraction of the true peak)", guide = guide_colourbar(barwidth = 12, barheight = 0.6, title.position = "top")) +
    ggtitle(sprintf("A = %s", a)) + theme_void(base_size = 11) + theme(plot.title = element_text(hjust = 0.5, size = 11), legend.position = "bottom")
}
bc <- (hd(0.5) | hd(1.5)) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
p <- (pp / bc) + plot_layout(heights = c(1, 1.1)) + plot_annotation(tag_levels = "A")
ggsave(file.path(out_dir, "bgam_fit_spatial.png"), p, width = 6.5, height = 7.0, dpi = 200, bg = "white", device = "png", type = "cairo")
cat("diff limit:", dl, "\n")
