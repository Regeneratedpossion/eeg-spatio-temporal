# Thesis Figure 2.2 (fig:data_generation).
# Figure: how the simulated data are generated (Section 2.2). Uses the real generator of the simulation.
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(patchwork)})
source(here::here("code/00_config.R"))   # run from the project root
coords <- read.csv(thesis_path("fig_data/adjacency_edges.txt"), nrows = 61)
out_dir <- fig_dir(thesis_path("diagnostics"))
setwd(project_dir); Sys.setenv(PROJECT_ROOT = project_dir)
source(here::here("code/01_simulate/01_generate_data.R"))
set.seed(3)

# One participant of scenario 3 (A = 1.5, omega_A = 0.26)
params <- simulation_parameters(n_ppt = 20, scenario = "cond_3", signal_mode = "hard_signal")
sensors <- compute_spatial_weight_base(coords, params)
sensors <- apply_active_radius(sensors, params$sigma_space * sqrt(2 * log(2)), "half-max")
L <- compute_noise_covariance(sensors, params)$L
base <- generate_dataset_base(params, sensors, L)
pp <- base$participant_parameters
d <- base$base_data
d$w <- sensors$spatial_weight[match(d$channel, sensors$channel)]
tt <- sort(unique(d$time))

ord <- sensors |> arrange(spatial_distance) |> mutate(row = row_number())
n_act <- sum(ord$spatial_weight > 0)
cz <- "CZ"
ms <- function(df) df |> mutate(time = time * 1000)
dz <- d |> filter(channel == cz) |> ms()
p1 <- pp$participant[1]
noise_only <- dz |> filter(participant == p1, condition == "cond1", trial <= 5) |> transmute(trial, time, value = eeg_noise)
tr2 <- dz |> filter(participant == p1, condition == "cond2") |> transmute(trial, time, value = w * signal_unit + eeg_noise)
m2 <- tr2 |> group_by(time) |> summarise(value = mean(value), .groups = "drop")
m1 <- dz |> filter(participant == p1, condition == "cond1") |> group_by(time) |> summarise(value = mean(eeg_noise), .groups = "drop")
clean <- dz |> filter(condition == "cond2", trial == 1) |> transmute(participant, time, value = w * signal_unit)
p1d <- d |> filter(participant == p1) |> mutate(y = w * signal_unit + eeg_noise) |> group_by(channel, condition, time) |> summarise(y = mean(y), .groups = "drop") |>
  tidyr::pivot_wider(names_from = condition, values_from = y) |> transmute(channel, time = time * 1000, value = cond2 - cond1)
p1d$row <- ord$row[match(p1d$channel, ord$channel)]
gt <- build_ground_truth(params, sensors, pp)$continuous_effect |> filter(channel == cz) |> transmute(time = time * 1000, value = effect)
gt <- gt[gt$time %in% unique(dz$time), ]

th <- theme_minimal(base_size = 9.5) + theme(plot.title = element_text(face = "bold", size = 9.5), panel.grid.minor = element_blank())
lim <- 1.1 * params$bump_amplitude
heat <- function(df, title) ggplot(df, aes(time, row, fill = value)) + geom_raster() +
  geom_hline(yintercept = n_act + 0.5, linetype = "dashed", linewidth = 0.3) +
  geom_vline(xintercept = c(160, 340), linetype = "dashed", colour = "grey40", linewidth = 0.3) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-lim, lim), oob = scales::squish, name = expression(mu * V)) +
  scale_y_reverse(expand = c(0, 0)) + scale_x_continuous(expand = c(0, 0)) +
  labs(title = title, x = "Time (ms)", y = "Electrode (by distance from CZ)") +
  th + theme(axis.text.y = element_blank(), panel.grid = element_blank())
win <- geom_vline(xintercept = c(160, 340), linetype = "dashed", colour = "grey40", linewidth = 0.3)
ylab <- expression("Amplitude (" * mu * V * ")")

pA <- ggplot(sensors, aes(xproj, yproj)) +
  geom_point(aes(fill = spatial_weight), shape = 21, size = 2.4, colour = "grey40") +
  geom_point(data = filter(sensors, spatial_weight > 0), shape = 21, size = 3.6, colour = "black", stroke = 0.5) +
  scale_fill_gradient(low = "white", high = "#B2182B", limits = c(0, 1), name = expression(w[c])) +
  coord_equal() + labs(title = expression(bold("A  Spatial weight ") * w[c])) +
  th + theme(axis.text = element_blank(), axis.title = element_blank(), panel.grid = element_blank())

wf <- do.call(rbind, lapply(params$ronset, function(o) data.frame(time = Xf, tau = o, f = temporal_signal_at_onset(o))))
pB <- ggplot(wf, aes(time, f)) + geom_line(aes(group = tau), colour = "grey50", linewidth = 0.12) +
  stat_summary(fun = mean, geom = "line", colour = "#2166AC", linewidth = 0.8) + win +
  labs(title = expression(bold("B  Waveform ") * f(t * ";" ~ tau[p])), x = "Time (ms)", y = "f") + th

om <- c(0.26, 0.70)
amp <- expand.grid(x = seq(0.01, 3, by = 0.01), omega = om) |> mutate(dens = dlnorm(x, -omega^2 / 2, omega), omega = factor(sprintf("%.2f", omega)))
pC <- ggplot(amp, aes(x, dens, colour = omega)) + geom_line(linewidth = 0.7) + geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  scale_colour_manual(values = c("#5E3C99", "#B2ABD2"), name = expression(omega[A])) +
  labs(title = expression(bold("C  Amplitude ") * A[p]), x = expression(A[p] / A), y = "Density") +
  th + theme(legend.position = c(0.8, 0.75))

pD <- ggplot(tr2, aes(time, value, group = trial)) + win + geom_line(colour = "grey75", linewidth = 0.2) +
  geom_line(data = m1, aes(group = NULL), colour = "#2166AC", linewidth = 0.6) +
  geom_line(data = m2, aes(group = NULL), colour = "black", linewidth = 0.8) +
  labs(title = "D  One participant (CZ, 20 trials)", x = "Time (ms)", y = ylab) + th
pE <- heat(p1d, "E  Same participant, all electrodes")
pF <- ggplot(clean, aes(time, value, group = participant)) + win + geom_line(colour = "grey60", linewidth = 0.25) +
  geom_line(data = gt, aes(group = NULL), colour = "black", linewidth = 1) +
  labs(title = expression(bold("F  All participants and truth ") * bar(g)[c * "," * t]), x = "Time (ms)", y = ylab) + th
bar <- guides(fill = guide_colourbar(barwidth = 5, barheight = 0.5, title.position = "left", title.vjust = 0.8))
pA <- pA + bar + theme(legend.position = "bottom")
pE <- pE + bar + theme(legend.position = "bottom")
fig <- (pA | pB | pC) / (pD | pE | pF)
ggsave(file.path(out_dir, "fig_data_generation.png"), fig, width = 8, height = 5.8, dpi = 200, bg = "white")
cat("participant:", unlist(pp[1, ]), "\n")
