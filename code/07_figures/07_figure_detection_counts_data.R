source(here::here("code/00_config.R")); source(here::here("code/R/passed_fits.R"))   # run from the project root
out_dir <- fig_dir(thesis_path("diagnostics")); csv_dir <- fig_dir(thesis_path("fig_data"))
suppressMessages({library(dplyr); library(ggplot2); library(patchwork)})
setwd(project_path("results"))
sc <- c("1" = 1, "2" = 2, "3" = 3, "4" = 4)
res <- list()
for (s in names(sc)) for (d in sort(pgs(sprintf("single_node_cond%s_ds_cr_seeds001_100/seed_*", s)))) for (f in c("bgam", "conv")) {
  x <- readRDS(file.path(d, "evaluation", sprintf("masks_%s.rds", f)))
  for (m in names(x$masks)) {
    mk <- x$masks[[m]]
    res[[length(res) + 1]] <- data.frame(A = sc[[s]], seed = basename(d), cond = s, method = m, time = x$times,
      active = rowSums(mk[, x$active, drop = FALSE]), inactive = rowSums(mk[, !x$active, drop = FALSE]))
  }
  if (f == "bgam") tr <- data.frame(time = x$times, active = rowSums(x$truth[, x$active, drop = FALSE]))
}
r <- bind_rows(res)
m <- r |> group_by(A, method, time) |> summarise(active = mean(active), inactive = mean(inactive), .groups = "drop")
write.csv(m, file.path(csv_dir, "detection_counts.csv"), row.names = FALSE)
cols <- c(BGAM = "#D55E00", cluster_mass = "#0072B2", TFCE = "#56B4E9", uncorrected = "#7F7F7F", BH = "#009E73", BY = "#CC79A7", Holm = "#E69F00")
lab <- c(BGAM = "BGAM", cluster_mass = "Cluster mass", TFCE = "TFCE", uncorrected = "Uncorrected", BH = "BH", BY = "BY", Holm = "Holm")
m <- m |> filter(method %in% names(cols)); m$method <- factor(m$method, names(cols))
tr$inactive <- 0
ttl <- c("1" = "Scenario 1\nA = 0.5, \u03c9A = 0.26", "2" = "Scenario 2\nA = 0.5, \u03c9A = 0.70", "3" = "Scenario 3\nA = 1.5, \u03c9A = 0.26", "4" = "Scenario 4\nA = 1.5, \u03c9A = 0.70")
mk <- function(a, col, ylab, ymax, first, top) {
  ggplot() + annotate("rect", xmin = 155, xmax = 345, ymin = -Inf, ymax = Inf, fill = "grey85", alpha = 0.5) +
    geom_line(data = filter(m, A == a, method != "BGAM"), aes(time * 1000, .data[[col]], colour = method), linewidth = 0.4, alpha = 0.7) +
    geom_line(data = filter(m, A == a, method == "BGAM"), aes(time * 1000, .data[[col]], colour = method), linewidth = 1.1) +
    geom_step(data = transform(tr, A = a), direction = "mid", aes(time * 1000, .data[[col]], linetype = "Truth"), colour = "black", linewidth = 0.45) +
    scale_colour_manual(values = cols, labels = lab, name = NULL) + scale_linetype_manual(values = "dashed", name = NULL) +
    scale_x_continuous(breaks = c(0, 250, 500)) +
    coord_cartesian(ylim = c(0, ymax)) + labs(title = if (top) ttl[[as.character(a)]] else NULL, x = "Time (ms)", y = if (first) ylab else NULL) +
    theme_minimal(base_size = 8) + theme(panel.grid.minor = element_blank(), plot.title = element_text(size = 8, face = "bold"))
}
ord <- 1:4; ym <- max(m$inactive) * 1.05
top <- lapply(seq_along(ord), function(i) mk(ord[i], "active", "Detected active electrodes", 9, i == 1, TRUE))
bot <- lapply(seq_along(ord), function(i) mk(ord[i], "inactive", "Detected inactive electrodes", ym, i == 1, FALSE))
p <- wrap_plots(c(top, bot), nrow = 2) + plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave(file.path(out_dir, "detection_counts.png"), p, width = 7.5, height = 5, dpi = 200, bg = "white")
