rd <- function(f) tryCatch({x <- read.csv(f); x[, intersect(names(x), c("rule","baseline_end","odds","method","alpha","F1","FPR","onset_error_ms","offset_error_ms"))]}, error = function(e) NULL)
mae <- function(e) mean(abs(e[is.finite(e)])); nv <- function(e) sum(is.finite(e))
cell <- function(c, m, t, x) data.frame(scenario = c, method = m, thr = t, n = nrow(x), F1 = mean(x$F1), FPR = mean(x$FPR),
  on_MAE = mae(x$onset_error_ms), off_MAE = mae(x$offset_error_ms), N_valid = nv(x$onset_error_ms))
out <- list()
for (c in c("I", "F", "H", "E")) {
  d <- sprintf("single_node_cond%s_ds_cr_seeds001_100", c)
  b <- do.call(rbind, lapply(Sys.glob(file.path(d, "seed_*/evaluation/evaluation_sensitivity.csv")), rd))
  b <- b[b$rule == "global" & abs(b$baseline_end - 0.10) < 1e-9, ]
  for (o in sort(unique(b$odds))) out[[length(out) + 1]] <- cell(c, "BGAM", o, b[b$odds == o, ])
  g <- do.call(rbind, lapply(Sys.glob(file.path(d, "seed_*/evaluation/evaluation_conventional.csv")), rd))
  for (m in unique(g$method)) for (a in sort(unique(g$alpha))) out[[length(out) + 1]] <- cell(c, m, a, g[g$method == m & g$alpha == a, ])
}
r <- do.call(rbind, out)
write.csv(r, "threshold_sweep_20261009.csv", row.names = FALSE)
cat(nrow(r), "rows; seeds per cell:", range(r$n), "\n")
