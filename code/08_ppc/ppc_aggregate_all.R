# Collects the PPC and residual-source results of all passed fits and lists the fits that did not pass.
# Run from the project root after all array tasks of slurm/08_ppc_array.sh are done.
root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = getwd()))
d <- Sys.getenv("PPC_OUT_ROOT", unset = file.path(root, "results", "ppc_checks"))
sel <- read.table(file.path(root, "results", "passed_fits.txt"), col.names = c("scenario", "seed"))

ppc <- do.call(rbind, lapply(seq_len(nrow(sel)), function(i) {
  f <- file.path(d, sprintf("scenario%d_seed_%03d", sel$scenario[i], sel$seed[i]), "ppc_summary.csv")
  stopifnot(file.exists(f))
  x <- read.csv(f); x$scenario <- sel$scenario[i]; x$seed <- sel$seed[i]; x
}))
ppc$outside90 <- ppc$below_90 | ppc$above_90
write.csv(ppc, file.path(d, "ppc_all_datasets.csv"), row.names = FALSE)

# Median [5%, 95% quantile over datasets]
q <- function(v) sprintf("%.3f [%.3f, %.3f]", median(v), quantile(v, .05), quantile(v, .95))
ppc_tab <- do.call(rbind, lapply(unique(ppc$statistic), function(st) do.call(rbind, lapply(1:4, function(sc) {
  z <- ppc[ppc$statistic == st & ppc$scenario == sc, ]
  data.frame(statistic = st, scenario = sc, n = nrow(z), observed = q(z$observed), rep_median = q(z$rep_median),
             outside90 = sum(z$outside90), share_outside90 = round(mean(z$outside90), 3))
}))))
write.csv(ppc_tab, file.path(d, "ppc_scenario_summary.csv"), row.names = FALSE)

rs <- do.call(rbind, lapply(seq_len(nrow(sel)), function(i) {
  f <- file.path(d, "residual_source", sprintf("scenario%d_seed_%03d.csv", sel$scenario[i], sel$seed[i]))
  stopifnot(file.exists(f)); read.csv(f)
}))
stopifnot(nrow(rs) == 2L * nrow(sel))
write.csv(rs, file.path(d, "residual_source_all_datasets.csv"), row.names = FALSE)
cols <- c("gt_noise_lag1", "bgam_lag1", "lag1_difference", "gt_noise_lag2", "bgam_lag2", "lag2_difference", "mean_rmse")
rs_tab <- do.call(rbind, lapply(1:4, function(sc) do.call(rbind, lapply(c("cond1", "cond2"), function(cn) {
  z <- rs[rs$scenario == sc & rs$condition == cn, ]
  data.frame(scenario = sc, condition = cn, n = nrow(z), setNames(lapply(cols, function(c) q(z[[c]])), cols))
}))))
write.csv(rs_tab, file.path(d, "residual_source_summary.csv"), row.names = FALSE)

# Fits that did not pass the convergence rule (R-hat < 1.01, bulk and tail ESS > 400)
files <- Sys.glob(file.path(root, "results", "scenario*", "seed_*", "evaluation", "evaluation_main.csv"))
ev <- do.call(rbind, lapply(files, function(f) {
  x <- read.csv(f)[1, c("max_rhat_g", "min_ess_bulk_g", "min_ess_tail_g")]
  data.frame(scenario = as.integer(sub(".*scenario([0-9]).*", "\\1", f)), seed = as.integer(sub(".*seed_([0-9]+).*", "\\1", f)), x)
}))
bad <- ev[!paste(ev$scenario, ev$seed) %in% paste(sel$scenario, sel$seed), ]
bad$reason <- apply(bad[c("max_rhat_g", "min_ess_bulk_g", "min_ess_tail_g")], 1, function(r)
  paste(c(if (r[1] >= 1.01) "R-hat >= 1.01", if (r[2] <= 400) "bulk ESS <= 400", if (r[3] <= 400) "tail ESS <= 400"), collapse = "; "))
bad <- bad[order(bad$scenario, bad$seed), ]
write.csv(bad, file.path(d, "not_passed_fits.csv"), row.names = FALSE)
cat("evaluated:", nrow(ev), "| passed:", nrow(sel), "| not passed:", nrow(bad), "\n")
print(table(bad$scenario)); print(ppc_tab, row.names = FALSE); print(rs_tab, row.names = FALSE)
