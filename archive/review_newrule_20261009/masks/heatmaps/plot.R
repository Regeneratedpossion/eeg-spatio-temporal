library(ggplot2)
ord <- c("truth","BGAM","uncorrected","BH","BY","Holm","cluster_mass","TFCE")
for (sc in c("I","F","H","E")) {
  acc <- list(); n <- 0
  for (d in sort(Sys.glob(sprintf("single_node_cond%s_ds_cr_seeds001_100/seed_*", sc)))) {
    for (f in c("bgam","conv")) {
      x <- readRDS(file.path(d, "evaluation", sprintf("masks_%s.rds", f)))
      for (m in names(x$masks)) acc[[m]] <- if (is.null(acc[[m]])) x$masks[[m]] * 1 else acc[[m]] + x$masks[[m]]
      if (f == "bgam") { tr <- x$truth; times <- x$times; chs <- x$channels; n <- n + 1 }
    }
  }
  acc$truth <- tr * n
  chorder <- chs[order(-colSums(tr), chs)]
  df <- do.call(rbind, lapply(ord, function(m) {
    a <- acc[[m]]
    data.frame(method = m, time = rep(times * 1000, times = length(chs)),
               channel = rep(chs, each = length(times)), freq = as.vector(a) / n)
  }))
  df$method <- factor(df$method, ord)
  df$channel <- factor(df$channel, rev(chorder))
  p <- ggplot(df, aes(time, channel, fill = freq)) + geom_raster() +
    facet_wrap(~method, nrow = 2) +
    scale_fill_viridis_c(limits = c(0, 1), name = "Detection\nfrequency") +
    labs(title = sprintf("Scenario %s: detection frequency over %d seeds", sc, n), x = "Time (ms)", y = "Electrode (true-signal electrodes at top)") +
    theme_minimal(base_size = 8) + theme(axis.text.y = element_blank(), panel.grid = element_blank())
  ggsave(sprintf("heatmaps/heatmap_%s.png", sc), p, width = 11, height = 6, dpi = 150)
}
