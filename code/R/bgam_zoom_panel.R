if (!exists("thesis_path")) source(here::here("code/00_config.R"))   # run from the project root
if (!exists("pgs")) source(here::here("code/R/passed_fits.R"))
out_dir <- fig_dir(thesis_path("diagnostics")); truth_file <- thesis_path("fig_data/truth_250.csv")
setwd(project_dir)
mk <- function(c, scen) { ps <- .man$seed[.man$cond == c]
  ce <- do.call(rbind, lapply(ps, function(s) readRDS(sprintf("data/scenario%s/gt_%03d.rds", c, s))$continuous_effect[, c("channel", "time", "effect")]))
  ce <- aggregate(effect ~ channel + time, ce, mean); ce$scenario <- scen; ce }
truth4 <- rbind(mk("1", 1), mk("2", 2), mk("3", 3), mk("4", 4))
cat("check CZ 250ms:", truth4$effect[truth4$channel == "CZ" & abs(truth4$time - 0.25) < 1e-9 & truth4$scenario %in% c(1, 3)], "\n")
suppressMessages({library(dplyr); library(ggplot2)})
sc <- c("1" = 1, "2" = 2, "3" = 3, "4" = 4)
d <- bind_rows(lapply(names(sc), function(c) bind_rows(lapply(pgs(sprintf("results/scenario%s/seed_*/fixed_grid/six_channel.csv", c)), read.csv)) |> mutate(scenario = sc[[c]]))) |>
  filter(channel == "CZ") |> mutate(A = ifelse(scenario %in% c(1, 2), 0.5, 1.5))
est <- d |> group_by(A, time) |> summarise(v = mean(model_mean_diff), .groups = "drop")
tr <- truth4 |> filter(channel == "CZ", scenario %in% c(1, 3)) |> mutate(A = ifelse(scenario == 1, 0.5, 1.5)) |> select(A, time, v = effect)
dl <- data.frame(A = c(0.5, 1.5), v = c(0.045, 0.051))
win <- function(x) x |> mutate(w = factor(ifelse(time <= 0.2, "Onset", ifelse(time >= 0.3, "Offset", NA)), c("Onset", "Offset"))) |> filter(!is.na(w), time >= 0.12, time <= 0.38) |> mutate(row = factor(sprintf("A = %s", A), c("A = 0.5", "A = 1.5")))
est <- win(est); tr <- win(tr)
dl <- merge(dl, data.frame(w = factor(c("Onset", "Offset"), c("Onset", "Offset")))) |> mutate(row = factor(sprintf("A = %s", A), c("A = 0.5", "A = 1.5")))
vl <- data.frame(w = factor(c("Onset", "Offset"), c("Onset", "Offset")), x = c(150, 350))
pz <- ggplot() +
  geom_vline(data = vl, aes(xintercept = x), colour = "grey75", linewidth = 0.6) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_hline(data = dl, aes(yintercept = v, linetype = "δ (mean SESOI)"), colour = "grey30", linewidth = 0.4) +
  geom_line(data = tr, aes(time * 1000, v, linetype = "Truth"), colour = "black", linewidth = 0.5) +
  geom_line(data = est, aes(time * 1000, v, linetype = "BGAM estimate"), colour = "#4575b4", linewidth = 0.8) +
  geom_point(data = est, aes(time * 1000, v), colour = "#4575b4", size = 1.6) +
  scale_linetype_manual(values = c("BGAM estimate" = "solid", "Truth" = "dashed", "δ (mean SESOI)" = "dotted"), name = NULL) +
  guides(linetype = guide_legend(override.aes = list(colour = c("#4575b4", "black", "grey30")))) +
  facet_grid(row ~ w, scales = "free_x") + coord_cartesian(ylim = c(-0.05, 0.25)) +
  labs(x = "Time (ms)", y = "Condition difference at CZ (µV)") +
  theme_minimal(base_size = 11) + theme(legend.position = "bottom", strip.text = element_text(face = "bold"), panel.grid.minor = element_blank())
if (sys.nframe() == 0) ggsave(file.path(out_dir, "zoom.png"), pz, width = 6.5, height = 4.2, dpi = 200, bg = "white")
