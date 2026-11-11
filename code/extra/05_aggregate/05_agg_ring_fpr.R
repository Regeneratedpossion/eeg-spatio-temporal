# Run from the project root. Reads results/single_node_cond*/seed_*/evaluation; writes to results/fig_data (AGG_OUT overrides).
source(here::here("code/R/aggregate_utils.R"))
res_dir <- here::here("results"); out_dir <- aggregate_out_dir()
rd <- function(f) tryCatch({x <- read.csv(f); if (!("FPR_ring1" %in% names(x))) return(NULL); x[, c(intersect(names(x), c("method","main")), "FPR_ring1","FPR_ring2","FPR_ring3")]}, error = function(e) NULL)
out <- list()
for (c in c("1", "3", "2", "4")) {
  d <- scenario_results_dir(res_dir, c)
  b <- do.call(rbind, lapply(Sys.glob(file.path(d, "seed_*/evaluation/evaluation_main.csv")), rd)); b$method <- "BGAM"
  g <- do.call(rbind, lapply(Sys.glob(file.path(d, "seed_*/evaluation/evaluation_conventional.csv")), rd)); g <- g[g$main, ]
  x <- rbind(b[, c("method","FPR_ring1","FPR_ring2","FPR_ring3")], g[, c("method","FPR_ring1","FPR_ring2","FPR_ring3")])
  x <- x[!is.na(x$FPR_ring1), ]
  for (m in unique(x$method)) { y <- x[x$method == m, ]
    out[[length(out)+1]] <- data.frame(scenario = c, method = m, n = nrow(y), ring1 = mean(y$FPR_ring1), ring2 = mean(y$FPR_ring2), ring3 = mean(y$FPR_ring3)) }
}
r <- do.call(rbind, out); write.csv(r, file.path(out_dir, "ring_fpr_20261009.csv"), row.names = FALSE); print(range(r$n))
