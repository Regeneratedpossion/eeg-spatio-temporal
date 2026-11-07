# Run from the project root. Reads the masks written by code/05_save_masks.R (results/<batch>/seed_*/evaluation).
source(here::here("code/00_config.R")); source(here::here("code/R/passed_fits.R"))
out_dir <- fig_dir(here::here("results/heatmaps"))
setwd(here::here("results"))
suppressMessages({library(ggplot2); library(dplyr)})
source(here::here("code/R/utils.R"))
sens <- read.csv(here::here("resources/electrodes/electrode_layout_61ch.csv"))
hr <- max(sqrt(sens$xproj^2 + sens$yproj^2)) * 1.05
dec <- scalp_decorations(hr)
ord <- c("truth","BGAM","uncorrected","BH","BY","Holm","cluster_mass","TFCE")
lab <- c(truth="Active set (truth)", BGAM="BGAM", uncorrected="Uncorrected", BH="BH", BY="BY", Holm="Holm", cluster_mass="Cluster mass", TFCE="TFCE")
tms <- c(0.16, 0.25, 0.34)
for (sc in c("I","F","H","E")) {
  acc <- list(); n <- 0
  for (d in sort(pgs(sprintf("single_node_cond%s_ds_cr_seeds001_100/seed_*", sc)))) for (f in c("bgam","conv")) {
    x <- readRDS(file.path(d, "evaluation", sprintf("masks_%s.rds", f)))
    for (m in names(x$masks)) acc[[m]] <- if (is.null(acc[[m]])) x$masks[[m]] * 1 else acc[[m]] + x$masks[[m]]
    if (f == "bgam") { tr <- x$truth; times <- x$times; chs <- x$channels; n <- n + 1 }
  }
  acc$truth <- tr * n
  ti <- sapply(tms, \(t) which.min(abs(times - t)))
  df <- do.call(rbind, lapply(ord, \(m) do.call(rbind, lapply(seq_along(ti), \(k)
    data.frame(method = m, time = sprintf("%.0f ms", tms[k] * 1000), channel = chs, freq = acc[[m]][ti[k], ] / n)))))
  df <- left_join(df, sens, by = "channel")
  df$method <- factor(lab[df$method], lab[ord])
  p <- ggplot(df, aes(xproj, yproj)) +
    geom_path(data = dec$head_circle, aes(x, y), inherit.aes = FALSE, linewidth = 0.3) +
    geom_path(data = dec$nose, aes(x, y), inherit.aes = FALSE, linewidth = 0.3) +
    geom_path(data = dec$ear_left, aes(x, y), inherit.aes = FALSE, linewidth = 0.3) +
    geom_path(data = dec$ear_right, aes(x, y), inherit.aes = FALSE, linewidth = 0.3) +
    geom_point(aes(fill = freq), shape = 21, size = 1.6, stroke = 0.15, color = "grey30") +
    facet_grid(method ~ time, switch = "y") +
    coord_equal(xlim = c(-hr * 1.12, hr * 1.12), ylim = c(-hr * 1.12, hr * 1.18), expand = FALSE) +
    scale_fill_viridis_c(limits = c(0, 1), name = "Detection\nfrequency") +
    labs(title = sprintf("Scenario %d", match(sc, c("I","H","F","E")))) +
    theme_void(base_size = 8) +
    theme(strip.text.x = element_text(face = "bold"), strip.text.y.left = element_text(angle = 0, hjust = 1, face = "bold"),
          plot.background = element_rect(fill = "white", color = NA), plot.title = element_text(size = 8))
  ggsave(file.path(out_dir, sprintf("topo_%s.png", sc)), p, width = 14, height = 18, units = "cm", dpi = 300)
}
