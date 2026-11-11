source(here::here("code/R/aggregate_utils.R"))
man <- read_passed_fits()
sc <- c("1" = 1, "2" = 2, "3" = 3, "4" = 4)
for (c in names(sc)) {
  d <- scenario_results_dir(here::here("results"), c); s <- man$seed[man$cond == c]
  g <- function(ev, m, v) sapply(s, function(k) { y <- read.csv(sprintf("%s/seed_%03d/%s/evaluation_conventional.csv", d, k, ev)); y[y$main & y$method == m, v] })
  for (v in c("F1", "onset_error_ms")) for (ev in c("evaluation", "evaluation_no_space")) {
    x <- g(ev, "TFCE", v) - g(ev, "cluster_mass", v); ok <- !is.na(x)
    cat(sprintf("scen %d %-15s %-20s n=%d diff=%.4f mcse=%.4f\n", sc[[c]], v, ev, sum(ok), mean(x[ok]), sd(x[ok]) / sqrt(sum(ok))))
  }
  # NA counts in no-space onset
  cat(sprintf("   NA onset (cm/tfce, no_space): %d/%d\n", sum(is.na(g("evaluation_no_space","cluster_mass","onset_error_ms"))), sum(is.na(g("evaluation_no_space","TFCE","onset_error_ms")))))
}
