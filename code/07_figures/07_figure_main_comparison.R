source(here::here("code/00_config.R")); source(here::here("code/R/passed_fits.R"))   # run from the project root
# Builds the main comparison figure and the trade-off figure from the per-seed evaluation CSVs.
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})

data_dir <- thesis_path("fig_data")
out_dir <- fig_dir(thesis_path("diagnostics"))
scen_dir <- c("1" = "I", "2" = "H", "3" = "F", "4" = "E")
scen_lab <- c("1" = "'Scenario 1:'~A == '0.5'*','~omega[A] == '0.26'", "2" = "'Scenario 2:'~A == '0.5'*','~omega[A] == '0.70'",
              "3" = "'Scenario 3:'~A == '1.5'*','~omega[A] == '0.26'", "4" = "'Scenario 4:'~A == '1.5'*','~omega[A] == '0.70'")

methods <- c("BGAM", "cluster_mass", "TFCE", "uncorrected", "BH", "BY", "Holm")
method_lab <- c(BGAM = "BGAM", cluster_mass = "Cluster mass", TFCE = "TFCE",
                uncorrected = "Uncorrected", BH = "BH", BY = "BY", Holm = "Holm")
cols <- c(BGAM = "#D55E00", cluster_mass = "#0072B2", TFCE = "#56B4E9", uncorrected = "#7F7F7F",
          BH = "#009E73", BY = "#CC79A7", Holm = "#E69F00")
shps <- c(BGAM = 16, cluster_mass = 15, TFCE = 17, uncorrected = 18, BH = 3, BY = 4, Holm = 8)

keep <- c("rule", "baseline_end", "odds", "method", "alpha", "F1", "FPR", "onset_error_ms", "offset_error_ms")
read_one <- function(f) {x <- read.csv(f); x[, intersect(names(x), keep)]}

load_scenario <- function(s) {
  d <- file.path(data_dir, sprintf("single_node_cond%s_ds_cr_seeds001_100", scen_dir[[s]]))
  b <- do.call(rbind, lapply(pgs(file.path(d, "seed_*/evaluation/evaluation_sensitivity.csv")), read_one))
  b <- b[b$rule == "global" & abs(b$baseline_end - 0.10) < 1e-9, ]
  b$method <- "BGAM"; b$thr <- b$odds
  g <- do.call(rbind, lapply(pgs(file.path(d, "seed_*/evaluation/evaluation_conventional.csv")), read_one))
  g$thr <- g$alpha
  cols_keep <- c("method", "thr", "F1", "FPR", "onset_error_ms", "offset_error_ms")
  out <- rbind(b[, cols_keep], g[, cols_keep])
  out$scenario <- s
  out
}

all_data <- do.call(rbind, lapply(names(scen_dir), load_scenario))

rmse_se <- function(e) {
  e <- e[is.finite(e)]; n <- length(e)
  r <- sqrt(mean(e^2))
  c(est = r, se = sd(e^2) / (2 * r * sqrt(n)), n = n)
}
mean_se <- function(x) c(est = mean(x), se = sd(x) / sqrt(length(x)), n = length(x))

summarise_cell <- function(d) {
  on <- rmse_se(d$onset_error_ms); off <- rmse_se(d$offset_error_ms)
  f1 <- mean_se(d$F1); fpr <- mean_se(d$FPR)
  data.frame(
    metric = c("Onset RMSE (ms)", "Offset RMSE (ms)", "F1", "FPR"),
    est = c(on["est"], off["est"], f1["est"], fpr["est"]),
    se = c(on["se"], off["se"], f1["se"], fpr["se"]),
    n_valid = c(on["n"], off["n"], NA, NA), n = nrow(d))
}

summ <- all_data |>
  group_by(scenario, method, thr) |>
  group_modify(~ summarise_cell(.x)) |>
  ungroup()

# Main setting: BGAM at odds 10, other methods at alpha = 0.05.
main <- summ |> filter((method == "BGAM" & thr == 10) | (method != "BGAM" & abs(thr - 0.05) < 1e-9))
main$method <- factor(main$method, levels = methods)
main$metric <- factor(main$metric, levels = c("Onset RMSE (ms)", "Offset RMSE (ms)", "F1", "FPR"))

fig_a <- function(scens, file, height = 7.6) {
  d <- main |> filter(scenario %in% scens)
  d$scenario <- factor(scen_lab[d$scenario], levels = scen_lab[scens])
  d$lab <- ifelse(!is.na(d$n_valid) & d$n_valid < d$n, paste0("n=", d$n_valid), NA)
  p <- ggplot(d, aes(method, est, colour = method, shape = method)) +
    geom_errorbar(aes(ymin = est - 1.96 * se, ymax = est + 1.96 * se), width = 0.0, linewidth = 0.5) +
    geom_point(size = 2.4, fill = "white") +
    geom_text(aes(label = lab), size = 2.4, vjust = -1.1, colour = "grey25", na.rm = TRUE) +
    facet_grid(metric ~ scenario, scales = "free_y", switch = "y", labeller = labeller(scenario = label_parsed, metric = label_value)) +
    scale_colour_manual(values = cols, labels = method_lab) +
    scale_shape_manual(values = shps, labels = method_lab) +
    scale_x_discrete(labels = method_lab) +
    scale_y_continuous(expand = expansion(mult = c(0.08, 0.18))) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 9) +
    theme(legend.position = "none", strip.placement = "outside", strip.background = element_blank(),
          strip.text = element_text(face = "bold"), panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(), axis.text.x = element_text(angle = 40, hjust = 1))
  ggsave(file.path(out_dir, file), p, width = 6.6, height = height)
}

# Range of F1, FPR and timing MAE over all thresholds, one row per scenario.
fig_c <- function(scens, file, height = 7.4) {
  mae <- function(d) mean(abs(c(d$onset_error_ms, d$offset_error_ms)), na.rm = TRUE)
  d <- all_data |> filter(scenario %in% scens) |>
    group_by(scenario, method, thr) |>
    summarise(F1 = mean(F1), FPR = mean(FPR), MAE = mae(pick(everything())), .groups = "drop") |>
    pivot_longer(c(F1, FPR, MAE), names_to = "metric", values_to = "est")
  allpts <- d
  d <- d |> group_by(scenario, method, metric) |>
    summarise(lo = min(est), hi = max(est),
              main = est[(method == "BGAM" & thr == 10) | (method != "BGAM" & abs(thr - 0.05) < 1e-9)],
              strict = est[if (first(method) == "BGAM") thr == max(thr) else thr == min(thr)],
              .groups = "drop")
  d$method <- factor(d$method, levels = rev(methods))
  m_lab <- c(F1 = "'F1'", FPR = "'False positive rate'", MAE = "'Timing MAE (ms)'")
  d$metric <- factor(m_lab[d$metric], levels = m_lab)
  r_lab <- c("1" = "A == '0.5'*','~omega[A] == '0.26'", "3" = "A == '0.5'*','~omega[A] == '0.70'",
             "2" = "A == '1.5'*','~omega[A] == '0.26'", "4" = "A == '1.5'*','~omega[A] == '0.70'")
  d$row <- factor(r_lab[as.character(d$scenario)], levels = r_lab[as.character(scens)])
  allpts$method <- factor(allpts$method, levels = rev(methods))
  allpts$metric <- factor(m_lab[allpts$metric], levels = m_lab)
  allpts$row <- factor(r_lab[as.character(allpts$scenario)], levels = r_lab[as.character(scens)])
  r_lab <- c("1" = "atop('Scenario 1',A == '0.5'*','~omega[A] == '0.26')", "2" = "atop('Scenario 2',A == '0.5'*','~omega[A] == '0.70')",
             "3" = "atop('Scenario 3',A == '1.5'*','~omega[A] == '0.26')", "4" = "atop('Scenario 4',A == '1.5'*','~omega[A] == '0.70')")
  d$row <- factor(r_lab[as.character(d$scenario)], levels = r_lab[as.character(scens)])
  allpts$method <- factor(allpts$method, levels = rev(methods))
  allpts$metric <- factor(m_lab[allpts$metric], levels = m_lab)
  allpts$row <- factor(r_lab[as.character(allpts$scenario)], levels = r_lab[as.character(scens)])
  main_plot <- ggplot(d, aes(y = method, colour = method)) +
    geom_segment(aes(x = lo, xend = hi, yend = method), linewidth = 1.2, lineend = "round") +
    geom_point(data = allpts, aes(x = est, y = method), size = 0.9, shape = 16, colour = "black", inherit.aes = FALSE) +
    geom_point(aes(x = strict), size = 2.2, shape = 1, stroke = 0.7) +
    geom_point(aes(x = main), size = 3, shape = 16) +
    facet_grid(row ~ metric, scales = "free_x", switch = "y", labeller = label_parsed) +
    scale_colour_manual(values = cols) +
    scale_y_discrete(labels = method_lab) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none", strip.placement = "outside",
          strip.text.y.left = element_text(angle = 0, hjust = 0, size = 10),
          strip.text.x = element_text(face = "bold"), panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank())
  key <- ggplot(data.frame(x = c(1, 4, 7), y = 1, s = c(16, 16, 1), z = c(0.9, 3, 2.2), lab = c("threshold tested", "main setting", "strictest threshold")),
                aes(x, y)) +
    geom_point(aes(size = z, shape = s), colour = c("black", "grey30", "grey30"), stroke = 0.7) +
    geom_text(aes(label = lab, x = x + 0.3), hjust = 0, size = 3.4) +
    scale_size_identity() + scale_shape_identity() + xlim(0.7, 10) + theme_void()
  p <- patchwork::wrap_plots(main_plot, key, ncol = 1, heights = c(1, 0.06))
  ggsave(file.path(out_dir, file), p, width = 7.8, height = height)
}

# N_valid table for the appendix (rows: method, columns: scenario x boundary).
nvalid_table <- function(file) {
  t <- main |> filter(metric %in% c("Onset RMSE (ms)", "Offset RMSE (ms)")) |>
    mutate(b = ifelse(grepl("Onset", metric), "on", "off")) |>
    select(scenario, method, b, n_valid) |>
    pivot_wider(names_from = c(scenario, b), values_from = n_valid, names_sep = "_")
  t <- t[order(t$method), ]
  ord <- c("1_on", "1_off", "2_on", "2_off", "3_on", "3_off", "4_on", "4_off")
  lines <- apply(t, 1, function(r) paste0(method_lab[[as.character(r[["method"]])]], " & ",
                                          paste(trimws(r[ord]), collapse = " & "), " \\\\"))
  writeLines(lines, file.path(out_dir, file))
}

fig_a(c("1", "2"), "fig_main_comparison.pdf")
fig_a(c("3", "4"), "fig_main_comparison_appendix.pdf")
fig_c(c("1", "2", "3", "4"), "fig_tradeoff.pdf")
nvalid_table("nvalid_rows.tex")
write.csv(summ, file.path(fig_dir(data_dir), "figure_summary.csv"), row.names = FALSE)
cat("done\n")

# ---- False positive rate by ring (2 x 2 panels) ----
ring_cols <- c("FPR_ring1", "FPR_ring2", "FPR_ring3")
read_ring <- function(f, keep_extra) {x <- read.csv(f); x[, intersect(names(x), c(keep_extra, ring_cols))]}
ring_data <- do.call(rbind, lapply(names(scen_dir), function(s) {
  d <- file.path(data_dir, sprintf("single_node_cond%s_ds_cr_seeds001_100", scen_dir[[s]]))
  b <- do.call(rbind, lapply(pgs(file.path(d, "seed_*/evaluation/evaluation_sensitivity.csv")),
                             read_ring, keep_extra = c("rule", "baseline_end", "odds")))
  b <- b[b$rule == "global" & abs(b$baseline_end - 0.10) < 1e-9 & b$odds == 10, ]
  b$method <- "BGAM"
  g <- do.call(rbind, lapply(pgs(file.path(d, "seed_*/evaluation/evaluation_conventional.csv")),
                             read_ring, keep_extra = c("method", "alpha")))
  g <- g[abs(g$alpha - 0.05) < 1e-9, ]
  x <- rbind(b[, c("method", ring_cols)], g[, c("method", ring_cols)])
  x$scenario <- s
  x
}))
ring_long <- ring_data |>
  pivot_longer(all_of(ring_cols), names_to = "ring", values_to = "fpr") |>
  filter(!is.na(fpr)) |>
  group_by(scenario, method, ring) |>
  summarise(fpr = mean(fpr), .groups = "drop")
ring_long$ring <- factor(sub("FPR_ring", "Ring ", ring_long$ring), levels = paste("Ring", 1:3))
ring_long$method <- factor(ring_long$method, levels = methods)
col_lab <- c("A == '0.5'~(weak)", "A == '1.5'~(strong)")
row_lab <- c("omega[A] == '0.26'~(low~heterogeneity)", "omega[A] == '0.70'~(high~heterogeneity)")
ring_long$col <- factor(ifelse(ring_long$scenario %in% c("1", "2"), col_lab[1], col_lab[2]), levels = col_lab)
ring_long$row <- factor(ifelse(ring_long$scenario %in% c("1", "3"), row_lab[1], row_lab[2]), levels = row_lab)
p_ring <- ggplot(ring_long, aes(ring, fpr, colour = method, shape = method, group = method)) +
  geom_line(linewidth = 0.5) + geom_point(size = 2) +
  facet_grid(row ~ col, labeller = label_parsed) +
  scale_colour_manual(values = cols, labels = method_lab, name = NULL) +
  scale_shape_manual(values = shps, labels = method_lab, name = NULL) +
  guides(colour = guide_legend(nrow = 1), shape = guide_legend(nrow = 1)) +
  labs(x = NULL, y = "False positive rate") +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom", legend.key.size = unit(0.4, "cm"), strip.background = element_blank(),
        strip.text = element_text(face = "bold"), panel.grid.minor = element_blank())
ggsave(file.path(out_dir, "ring_fpr.pdf"), p_ring, width = 6.6, height = 5.6)
