source(here::here("code/R/aggregate_utils.R"))
man <- read_passed_fits()
sc <- c("1" = 1, "2" = 2, "3" = 3, "4" = 4)
out <- list()
for (c in names(sc)) for (ev in c("evaluation", "evaluation_no_space")) {
  d <- scenario_results_dir(here::here("results"), c)
  f <- sprintf("%s/seed_%03d/%s/evaluation_conventional.csv", d, man$seed[man$cond == c], ev)
  x <- do.call(rbind, lapply(f, function(p) { y <- read.csv(p); y[y$main & y$method %in% c("cluster_mass", "TFCE"), ] }))
  a <- aggregate(cbind(F1, precision, recall, FPR, onset_error_ms, offset_error_ms) ~ method, x, mean)
  a$n <- as.vector(table(x$method)[a$method]); a$scenario <- sc[[c]]; a$adjacency <- ifelse(ev == "evaluation", "space+time", "time only")
  out[[length(out) + 1]] <- a
}
r <- do.call(rbind, out); r[3:8] <- round(r[3:8], 3)
write.csv(r, "results/nospace_summary.csv", row.names = FALSE)
print(r[order(r$scenario, r$method, r$adjacency), c("scenario","method","adjacency","n","F1","precision","recall","FPR","onset_error_ms","offset_error_ms")], row.names = FALSE)
