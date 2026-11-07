###########################################################
# Laplacian workflow configuration.
#
# Fixed model used here:
#   P20 Gaussian hardlabel data
#   K=20 graph-Laplacian basis, including the constant component
#   Bayesian temporal GAMM on each spatial basis score
###########################################################

library(brms)
library(dplyr)
library(tidyr)
library(posterior)

k_basis <- 20L
task_ids <- seq_len(k_basis)

data_subdir <- "montecarlo_p20_gaussian_hardlabel"
out_subdir <- "laplacian_basis_brms_p20_gaussian_hardlabel_constant_long/k20"
result_tag <- "laplacian_k20_seed001_iter3000_warmup1000_adapt0p999"

seed <- 1L
iter <- 3000L
warmup <- 1000L
time_k <- 10L
adapt_delta <- 0.999
max_treedepth <- 13L
cz_neighbours <- c("CZ", "FCZ", "CPZ", "C1", "C2", "FC1", "FC2",
                   "CP1", "CP2", "Cz", "FCz", "CPz")

model_formula <- bf(score | se(score_sem) ~
                        condition +
                        s(time, by = condition, k = time_k) +
                        (1 + condition | participant))
