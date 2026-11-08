source(here::here("code/00_config.R")); source(here::here("code/R/passed_fits.R"))   # run from the project root
out_dir <- fig_dir(thesis_path("diagnostics")); truth_file <- thesis_path("fig_data/truth_250.csv")
setwd(project_dir)
suppressMessages({library(dplyr); library(ggplot2); library(patchwork)})
source("code/R/utils.R")
sc <- c(I = 1, H = 2, F = 3, E = 4)
est <- bind_rows(lapply(names(sc), function(c) bind_rows(lapply(pgs(sprintf("results/single_node_cond%s_ds_cr_seeds001_100/seed_*/fixed_grid/topomap_3t.csv", c)), read.csv)) |> mutate(scenario = sc[[c]]))) |>
  filter(abs(time - 0.25) < 1e-9) |>
  group_by(scenario, channel, xproj, yproj) |> summarise(est = mean(model_mean_diff), .groups = "drop")
d <- left_join(est, read.csv(truth_file)[, c("scenario", "channel", "effect")], by = c("scenario", "channel"))
peak <- d |> group_by(scenario) |> summarise(pk = max(effect)); d <- left_join(d, peak, by = "scenario")
d <- d |> mutate(Truth = effect / pk, Estimate = est / pk, Diff = (est - effect) / pk)
active <- readRDS("data/single_node_condI_seeds001_100/gt_001.rds")$active_electrodes
d$active <- d$channel %in% active
sens <- distinct(d, channel, xproj, yproj)
hr <- max(sqrt(sens$xproj^2 + sens$yproj^2)) * 1.05
dec <- scalp_decorations(hr)
grid <- expand.grid(xproj = seq(-hr, hr, length.out = 110), yproj = seq(-hr, hr, length.out = 110)) |> filter(sqrt(xproj^2 + yproj^2) <= hr)
order_sc <- 1:4
act <- d |> filter(active) |> distinct(channel, xproj, yproj)
lab <- c("1" = "Scenario~1~(A == 0.5*','~omega[A] == \"0.26\")", "2" = "Scenario~2~(A == 0.5*','~omega[A] == \"0.70\")",
         "3" = "Scenario~3~(A == 1.5*','~omega[A] == \"0.26\")", "4" = "Scenario~4~(A == 1.5*','~omega[A] == \"0.70\")")
mk <- function(col, name) bind_rows(lapply(order_sc, function(sn) {
  e <- filter(d, scenario == sn); g <- grid
  g$v <- idw_interp(g, e[, c("xproj", "yproj")], e[[col]], power = 2, k = 6)
  g$scen <- factor(lab[[as.character(sn)]], lab[as.character(order_sc)]); g$panel <- name; g }))
g1 <- bind_rows(mk("Truth", "Truth"), mk("Estimate", "BGAM estimate")); g1$panel <- factor(g1$panel, c("Truth", "BGAM estimate"))
g2 <- mk("Diff", "Estimate \u2212 truth")
base <- function(g, lim, nm, strip_left) {
  ggplot(g, aes(xproj, yproj)) + geom_raster(aes(fill = v), interpolate = TRUE) +
    geom_path(data = dec$head_circle, aes(x, y), inherit.aes = FALSE, linewidth = 0.3) +
    geom_path(data = dec$nose, aes(x, y), inherit.aes = FALSE, linewidth = 0.3) +
    geom_path(data = dec$ear_left, aes(x, y), inherit.aes = FALSE, linewidth = 0.3) +
    geom_path(data = dec$ear_right, aes(x, y), inherit.aes = FALSE, linewidth = 0.3) +
    geom_point(data = sens, aes(xproj, yproj), inherit.aes = FALSE, size = 0.25, colour = "black") +
    geom_point(data = act, aes(xproj, yproj), inherit.aes = FALSE, shape = 1, size = 1.8, stroke = 0.4, colour = "black") +
    facet_grid(scen ~ panel, switch = "y", labeller = labeller(scen = label_parsed)) +
    coord_equal(xlim = c(-hr * 1.12, hr * 1.12), ylim = c(-hr * 1.12, hr * 1.18), expand = FALSE) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b", midpoint = 0, limits = c(-lim, lim), name = nm,
                         guide = guide_colourbar(barwidth = 8, barheight = 0.5, title.position = "top")) +
    theme_void(base_size = 9) +
    theme(legend.position = "bottom", strip.text.x = element_text(face = "bold", size = 9, margin = margin(b = 3)),
          strip.text.y.left = if (strip_left) element_text(angle = 90, size = 8) else element_blank())
}
dl <- ceiling(max(abs(g2$v), na.rm = TRUE) * 20) / 20
p <- base(g1, 1, "Truth and estimate (fraction of the true peak)", TRUE) | base(g2, dl, "Difference (fraction of the true peak)", FALSE)
p <- p + plot_layout(widths = c(2.25, 1))
ggsave(file.path(out_dir, "bgam_fit_spatial_scenarios.png"), p, width = 7.5, height = 9.5, dpi = 200, bg = "white")
