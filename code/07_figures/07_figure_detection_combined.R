# Figure 6: detected electrodes over time (A, B) and detection rate by distance from CZ (C), converged fits only.
source(here::here("code/00_config.R")); source(here::here("code/R/passed_fits.R"))   # run from the project root
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(patchwork)})
thesis <- thesis_dir; out_dir <- fig_dir(thesis_path("diagnostics")); csv_dir <- fig_dir(thesis_path("fig_data"))
xy <- read.csv(file.path(thesis, "fig_data/adjacency_edges.txt"), nrows = 61)
setwd(project_path("results"))   # the masks are in results/<batch>/seed_*/evaluation/

sc <- c(I = 1, H = 2, F = 3, E = 4)
res <- list(); rate <- list()
for (s in names(sc)) for (d in sort(pgs(sprintf("single_node_cond%s_ds_cr_seeds001_100/seed_*", s)))) for (f in c("bgam", "conv")) {
  x <- readRDS(file.path(d, "evaluation", sprintf("masks_%s.rds", f)))
  win <- x$times >= 0.16 - 1e-9 & x$times <= 0.34 + 1e-9
  for (m in names(x$masks)) {
    mk <- x$masks[[m]]
    res[[length(res) + 1]] <- data.frame(A = sc[[s]], method = m, time = x$times,
      active = rowSums(mk[, x$active, drop = FALSE]), inactive = rowSums(mk[, !x$active, drop = FALSE]))
    rate[[length(rate) + 1]] <- data.frame(A = sc[[s]], method = m, channel = x$channels, rate = colMeans(mk[win, , drop = FALSE]))
  }
  if (f == "bgam") tr <- data.frame(time = x$times, active = rowSums(x$truth[, x$active, drop = FALSE]))
}
m <- bind_rows(res) |> group_by(A, method, time) |> summarise(active = mean(active), inactive = mean(inactive), .groups = "drop")
tr$inactive <- 0

cz <- xy[xy$channel == "CZ", ]
xy$dist <- sqrt((xy$xproj - cz$xproj)^2 + (xy$yproj - cz$yproj)^2)
r_act <- 0.5 * sqrt(2 * log(2))
dr <- bind_rows(rate) |> left_join(xy[, c("channel", "dist")], by = "channel") |>
  mutate(dist = round(dist, 2)) |>
  group_by(A, method, dist) |> summarise(rate = mean(rate), .groups = "drop")

cols <- c(BGAM = "#D55E00", cluster_mass = "#0072B2", TFCE = "#56B4E9", uncorrected = "#7F7F7F",
          BH = "#009E73", BY = "#CC79A7", Holm = "#E69F00")
lab <- c(BGAM = "BGAM", cluster_mass = "Cluster mass", TFCE = "TFCE", uncorrected = "Uncorrected",
         BH = "BH", BY = "BY", Holm = "Holm")
ttl <- c("1" = "Scenario 1\nA = 0.5, ωA = 0.26", "2" = "Scenario 2\nA = 0.5, ωA = 0.70",
         "3" = "Scenario 3\nA = 1.5, ωA = 0.26", "4" = "Scenario 4\nA = 1.5, ωA = 0.70")
ord <- 1:4
m <- m |> filter(method %in% names(cols)); m$method <- factor(m$method, names(cols))
dr <- dr |> filter(method %in% names(cols)); dr$method <- factor(dr$method, names(cols))

th <- theme_minimal(base_size = 9.5) +
  theme(panel.grid.minor = element_blank(), plot.title = element_text(size = 9.5, face = "bold"))
cnt <- function(a, col, ylab, ymax, first, top) {
  ggplot() + geom_vline(xintercept = c(160, 340), linetype = "dashed", colour = "grey40", linewidth = 0.3) +
    geom_line(data = filter(m, A == a, method != "BGAM"), aes(time * 1000, .data[[col]], colour = method), linewidth = 0.4, alpha = 0.7) +
    geom_line(data = filter(m, A == a, method == "BGAM"), aes(time * 1000, .data[[col]], colour = method), linewidth = 1.1) +
    geom_step(data = tr, direction = "mid", aes(time * 1000, .data[[col]], linetype = "Truth"), colour = "black", linewidth = 0.45) +
    scale_colour_manual(values = cols, labels = lab, name = NULL) +
    scale_linetype_manual(values = "dashed", name = NULL) +
    scale_x_continuous(breaks = c(0, 250, 500)) + coord_cartesian(ylim = c(0, ymax)) +
    labs(title = if (top) ttl[[as.character(a)]] else NULL, x = "Time (ms)", y = if (first) ylab else NULL) + th
}
dist_row <- function(a, first) {
  d <- filter(dr, A == a)
  hi <- c("BGAM", "cluster_mass")
  ggplot(d, aes(dist, rate, colour = method)) +
    geom_vline(xintercept = r_act, linetype = "dashed", colour = "grey40", linewidth = 0.3) +
    geom_point(data = filter(d, !method %in% hi), colour = "grey80", size = 1.3, show.legend = FALSE) +
    geom_point(data = filter(d, method %in% hi), size = 1.9, alpha = 0.9, show.legend = FALSE) +
    annotate("text", x = r_act, y = 0.97, label = "active radius", hjust = -0.1, size = 2.6, colour = "grey30") +
    scale_colour_manual(values = cols, labels = lab, name = NULL) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Distance from CZ", y = if (first) "Detection rate" else NULL) + th
}
ym <- max(m$inactive) * 1.05
rowA <- lapply(seq_along(ord), function(i) cnt(ord[i], "active", "Detected active electrodes", 9, i == 1, TRUE))
rowB <- lapply(seq_along(ord), function(i) cnt(ord[i], "inactive", "Detected inactive electrodes", ym, i == 1, FALSE))
rowC <- lapply(seq_along(ord), function(i) dist_row(ord[i], i == 1))
p <- wrap_plots(c(rowA, rowB, rowC), nrow = 3) + plot_layout(guides = "collect") +
  plot_annotation(tag_levels = list(c("A", rep("", 3), "B", rep("", 3), "C", rep("", 3)))) &
  theme(legend.position = "bottom")
ggsave(file.path(out_dir, "detection_combined.png"), p, width = 7.5, height = 7.6, dpi = 200, bg = "white")
write.csv(m, file.path(csv_dir, "detection_counts.csv"), row.names = FALSE)
cat("n dist groups", length(unique(dr$dist)), "\n")
