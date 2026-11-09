out <- list()
for (c in c("1", "3", "2", "4")) {
  d <- sprintf("single_node_cond%s_ds_cr_seeds001_100", c)
  b <- do.call(rbind, lapply(Sys.glob(file.path(d, "seed_*/evaluation/detection_by_electrode.csv")), read.csv)); b$method <- "BGAM"
  g <- do.call(rbind, lapply(Sys.glob(file.path(d, "seed_*/evaluation/detection_by_electrode_conventional.csv")), read.csv))
  x <- rbind(b[, c("method","channel","rate")], g[, c("method","channel","rate")])
  a <- aggregate(rate ~ method + channel, x, mean); a$n <- as.integer(aggregate(rate ~ method + channel, x, length)$rate); a$scenario <- c
  out[[c]] <- a
}
r <- do.call(rbind, out); write.csv(r, "electrode_fpr_20261009.csv", row.names = FALSE); cat("rows", nrow(r), "n range", range(r$n), "\n")
