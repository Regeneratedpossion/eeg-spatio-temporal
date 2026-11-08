# Keeps only the paths of fits that passed the convergence criterion (Section 2.2.5).
if (!exists("thesis_path")) source(here::here("code/00_config.R"))
.man <- read.table(here::here("resources/passed_fits.txt"), col.names = c("cond", "seed"))
.key <- paste(.man$cond, .man$seed)
pg <- function(p) p[paste(sub(".*single_node_cond([A-Z]).*", "\\1", p), as.integer(sub(".*seed_([0-9]+).*", "\\1", p))) %in% .key]
pgs <- function(...) pg(Sys.glob(...))
