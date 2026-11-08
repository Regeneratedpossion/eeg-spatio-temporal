# Builds the result tables. Usage: Rscript make_tables.R [all|passed]
suppressMessages(library(dplyr))
mode <- commandArgs(TRUE)[1]; if (is.na(mode)) mode <- "passed"
source(here::here("code/00_config.R"))   # run from the project root: Rscript fig_data/tables/make_tables.R passed
tdir <- thesis_path("fig_data/tables"); odir <- fig_dir(tdir)
root0 <- thesis_path("fig_data")
sc <- c(I = 1, H = 2, F = 3, E = 4)
man <- read.table(here::here("resources/passed_fits.txt"), col.names = c("cond", "seed"))
rd <- function(c, f, root = root0) bind_rows(lapply(sort(Sys.glob(sprintf("%s/single_node_cond%s_ds_cr_seeds001_100/seed_*/evaluation/%s", root, c, f))), function(p) {
  x <- read.csv(p); x$seed <- as.integer(sub("seed_", "", basename(dirname(dirname(p))))); x$scenario <- sc[[c]]; x$cond <- c; x }))
b <- bind_rows(lapply(names(sc), rd, f = "evaluation_main.csv"))
g <- bind_rows(lapply(names(sc), rd, f = "evaluation_conventional.csv")) |> filter(main == TRUE)
bg <- b |> filter(rule == "global", odds == 10) |> mutate(method = "BGAM")
z <- bind_rows(lapply(names(sc), rd, f = "evaluation_sensitivity.csv", root = root0)) |> filter(rule == "zero", odds == 10)
if (mode == "passed") { key <- paste(man$cond, man$seed); f <- function(d) d[paste(d$cond, d$seed) %in% key, ]; bg <- f(bg); g <- f(g); z <- f(z) }
cat("fits per scenario:", table(bg$scenario), "\n")
all <- bind_rows(bg, g)
M <- c("BGAM", "cluster_mass", "TFCE", "uncorrected", "BH", "BY", "Holm")
lab <- c("Model-based & BGAM", "Cluster-based & Cluster mass", " & TFCE", "Cell-wise & Uncorrected", " & BH", " & BY", " & Holm")
ord <- 1:4
fmt <- function(x, d) { s <- formatC(x, format = "f", digits = d); sub("^-", "$-$", s) }
panel <- function(vals, d, best) { # vals: matrix methods x cols
  r <- round(vals, d); out <- character(length(M))
  for (i in seq_along(M)) out[i] <- paste(sapply(seq_len(ncol(vals)), function(j) {
    s <- fmt(vals[i, j], d); if (!is.na(best) && !is.na(r[i, j]) && r[i, j] == (if (best == "min") min(r[, j], na.rm = TRUE) else max(r[, j], na.rm = TRUE))) paste0("\\textbf{", s, "}") else s }), collapse = " & ")
  paste0("    ", lab, " & ", out, " \\\\")
}
tim <- function(stat, best) {
  v <- sapply(M, function(m) unlist(lapply(ord, function(s) { d <- all |> filter(method == m, scenario == s); sapply(c("onset_error_ms", "offset_error_ms"), function(k) { e <- d[[k]][!is.na(d$onset)]; stat(e) }) })))
  panel(t(v), 1, best)
}
st <- list(Bias = mean, MAE = function(e) mean(abs(e)), RMSE = function(e) sqrt(mean(e^2)), SD = sd)
tt <- c("\\begin{tabular}{llcccccccc}", "    \\toprule", NULL)
body <- c()
for (nm in names(st)) body <- c(body, "    \\midrule", sprintf("    \\multicolumn{10}{l}{\\textit{%s. %s}} \\\\", LETTERS[match(nm, names(st))], nm), tim(st[[nm]], if (nm == "Bias") NA else "min"))
writeLines(body, file.path(odir, sprintf("timing_%s.tex", mode)))
ov <- function(col, d, best, scale = 1, valid = FALSE) {
  v <- sapply(M, function(m) sapply(ord, function(s) { x <- all |> filter(method == m, scenario == s); y <- x[[col]]; if (valid) y <- y[(x$TP + x$FP) > 0]; mean(y, na.rm = TRUE) * scale }))
  panel(t(v), d, best)
}
body <- c("    \\midrule", "    \\multicolumn{6}{l}{\\textit{A. $F_1$}} \\\\", ov("F1", 3, "max"),
  "    \\midrule", "    \\multicolumn{6}{l}{\\textit{B. Precision}} \\\\", ov("precision", 3, "max", valid = TRUE),
  "    \\midrule", "    \\multicolumn{6}{l}{\\textit{C. Recall}} \\\\", ov("recall", 3, "max"),
  "    \\midrule", "    \\multicolumn{6}{l}{\\textit{D. False positive rate (per 1,000 cells)}} \\\\", ov("FPR", 2, "min", scale = 1000))
writeLines(body, file.path(odir, sprintf("overlap_%s.tex", mode)))
t3 <- sapply(c("map_cor", "recovery", "leakage"), function(k) sapply(ord, function(s) mean(bg[[k]][bg$scenario == s])))
writeLines(c(paste("    Map correlation &", paste(fmt(t3[, 1], 3), collapse = " & "), "\\\\"), paste("    Recovery &", paste(fmt(t3[, 2], 3), collapse = " & "), "\\\\"), paste("    Leakage &", paste(fmt(t3[, 3], 3), collapse = " & "), "\\\\")), file.path(odir, sprintf("bgam_%s.tex", mode)))
sesoi <- c(); sr <- read.csv(file.path(tdir, if (mode == "passed") "sesoi_rule_passed.csv" else "sesoi_rule_20261009.csv")); sr$s <- match(sr$scenario, c("I", "H", "F", "E"))
for (s in 1:4) for (r in c("zero", "global")) { d <- sr[sr$s == s & sr$rule == r, ]
  sesoi <- c(sesoi, sprintf("    %d & %s & %s & %s & %s & %s \\\\", s, if (r == "zero") "$\\delta = 0$" else "SESOI", fmt(d$FPR, 3), fmt(d$post_offset, 3), fmt(d$recall, 3), fmt(d$F1, 3))) }
writeLines(sesoi, file.path(odir, sprintf("sesoi_%s.tex", mode)))
