# Control: cluster mass and TFCE with no spatial neighbours (time-only adjacency).
# Same code as 03_evaluate_conventional.R; only the MNE adjacency is set to zero.
suppressPackageStartupMessages(library(reticulate))
project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = here::here()))
source(file.path(project_root, "code", "03_evaluate_conventional.R"))
mne_cell_pvalues_full <- mne_cell_pvalues
mne_cell_pvalues <- function(arr, adj, s, tfce) mne_cell_pvalues_full(arr, adj * 0, s, tfce)
a <- commandArgs(trailingOnly = TRUE)
out <- evaluate_conventional(a[[1]], a[[2]], a[[3]])
print(out[out$main, c("method", "F1", "onset_error_ms", "offset_error_ms")], row.names = FALSE, digits = 3)
