# Original: a scratch script (chk.R). Usage: PROJECT_ROOT=<project> Rscript make_div_vs_error.R passed_fits.txt div_vs_error.csv
suppressMessages(library(dplyr))
root <- Sys.getenv("PROJECT_ROOT"); stopifnot(nzchar(root))   # the simulation project (on the cluster or the local copy)
args <- commandArgs(TRUE); passed_file <- args[1]; out_file <- args[2]
passed <- read.table(passed_file, col.names = c("cond","seed"))
sc <- c("1" = 1, "2" = 2, "3" = 3, "4" = 4)
ex <- list(); dv <- list()
for (cd in names(sc)) for (s in 1:100) {
  d <- sprintf("%s/results/scenario%s/seed_%03d", root, cd, s)
  isp <- any(passed$cond == cd & passed$seed == s)
  if (!file.exists(file.path(d, "convergence.csv"))) next
  if (!isp) {
    f <- file.path(d, "evaluation/convergence_effect.csv")
    cv <- read.csv(file.path(d, "convergence.csv"))
    if (!file.exists(f)) { ex[[length(ex)+1]] <- data.frame(cond=cd, seed=s, have=FALSE, max_rhat=cv$max_rhat, min_tail=cv$min_ess_tail); next }
    x <- read.csv(f)
    gt <- readRDS(sprintf("%s/data/scenario%s/gt_%03d.rds", root, cd, s))
    act <- x$channel %in% gt$active_electrodes; win <- x$time >= 0.16-1e-9 & x$time <= 0.34+1e-9
    bad <- x$ess_tail < 400
    ex[[length(ex)+1]] <- data.frame(cond=cd, seed=s, have=TRUE, max_rhat=max(x$rhat), min_tail=min(x$ess_tail), n_bad=sum(bad), nbad_rhat=sum(x$rhat>1.01),
       f_act=mean(act[bad]), f_win=mean(win[bad]), n=nrow(x))
  } else {
    nu <- read.csv(file.path(d, "nuts_summary.csv"))
    m <- read.csv(file.path(d, "evaluation/evaluation_main.csv"))
    dv[[length(dv)+1]] <- data.frame(scenario = sc[[cd]], seed = s, divergences = nu$divergences, onset_error_ms = m$onset_error_ms[1], offset_error_ms = m$offset_error_ms[1])
  }
}
ex <- bind_rows(ex); dv <- bind_rows(dv)
cat("excluded fits:", nrow(ex), " with convergence_effect:", sum(ex$have), "\n")
e <- filter(ex, have)
cat("all tail<400 in >=1 cell:", all(e$n_bad >= 1), "; rhat>1.01:", sum(e$max_rhat > 1.01), "\n")
cat("median n_bad:", median(e$n_bad), " median min_tail:", median(e$min_tail), "\n")
cat("mean share on active electrodes:", mean(e$f_act), " in window:", mean(e$f_win), "\n")
cat("excluded per scenario:", table(sc[ex$cond]), "\n")
write.csv(dv, out_file, row.names = FALSE)
cat("passed fits:", nrow(dv), "with divergences:", sum(dv$divergences > 0), "\n")
for (s in 1:4) { d <- filter(dv, scenario == s); cat("scen", s, "n", nrow(d), "spearman onset", round(cor(d$divergences, d$onset_error_ms, method="spearman", use="complete.obs"),3), "offset", round(cor(d$divergences, d$offset_error_ms, method="spearman", use="complete.obs"),3), "\n") }
