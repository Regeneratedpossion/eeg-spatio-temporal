# Run from the project root. Reads results/single_node_cond*/seed_*/evaluation; writes to results/fig_data (AGG_OUT overrides).
res_dir <- here::here("results"); out_dir <- Sys.getenv("AGG_OUT", unset = here::here("results", "fig_data")); dir.create(file.path(out_dir, "tables"), showWarnings = FALSE, recursive = TRUE)
rd <- function(f) tryCatch({x <- read.csv(f); if (!("post_offset" %in% names(x))) return(NULL); x[x$odds == 10, c("rule","baseline_end","F1","FPR","recall","post_offset")]}, error = function(e) NULL)
out <- list()
for (c in c("I","F","H","E")) {
  d <- file.path(res_dir, sprintf("single_node_cond%s_ds_cr_seeds001_100", c))
  b <- do.call(rbind, lapply(Sys.glob(file.path(d, "seed_*/evaluation/evaluation_sensitivity.csv")), rd))
  for (r in c("zero","global")) { y <- b[b$rule == r & (r == "zero" | abs(b$baseline_end - 0.10) < 1e-9), ]
    out[[length(out)+1]] <- data.frame(scenario = c, rule = r, n = nrow(y), F1 = mean(y$F1), FPR = mean(y$FPR), post_offset = mean(y$post_offset), recall = mean(y$recall)) }
}
r <- do.call(rbind, out); write.csv(r, file.path(out_dir, "tables", "sesoi_rule_20261009.csv"), row.names = FALSE); print(r)
