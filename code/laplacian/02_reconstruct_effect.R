###########################################################
# Reconstruct sensor x time posterior effects.
#
# Each fitted model gives the condition effect for one basis score.
# Multiplying that effect by its spatial basis vector and summing over
# bases reconstructs the posterior effect map.
###########################################################

load_laplace_reconstruction <- function() {
    # Fitted basis models.
    fit_dir <- file.path(project_root(), "results", out_subdir)
    summary_file <- file.path(fit_dir, "summary_all.csv")
    conv_summary <- read.csv(summary_file)
    fit_files <- conv_summary$fit_file[order(conv_summary$basis_task)]

    # Ground truth and time grid from the simulation data.
    ground_truth <- readRDS(data_path(data_subdir, "gt_001.rds"))
    sim <- readRDS(data_path(data_subdir, "sim_001.rds"))
    sim_model <- aggregate_trials(sim)
    times <- sort(unique(sim_model$time))
    condition_levels <- c("cond1", "cond2")

    condition_score_draws <- function(obj, times) {
        participant_id <- c(obj$participant_id, obj$score_data$participant[1])[1]

        base <- data.frame(
            time = times,
            score = 0,
            score_sem = 1,
            participant = participant_id
        )
        nd1 <- base
        nd2 <- base
        nd1$condition <- factor("cond1", levels = condition_levels)
        nd2$condition <- factor("cond2", levels = condition_levels)

        posterior_epred(obj$fit, newdata = nd2, re_formula = NA) -
            posterior_epred(obj$fit, newdata = nd1, re_formula = NA)
    }

    first_obj <- readRDS(fit_files[1])
    channels <- first_obj$channels
    first_score_diff <- condition_score_draws(first_obj, times)
    effect_array <- array(
        0,
        dim = c(nrow(first_score_diff), length(times), length(channels)),
        dimnames = list(NULL, as.character(times), channels)
    )

    add_basis_effect <- function(obj, score_diff) {
        sweep(array(score_diff, dim = dim(effect_array)),
              3, obj$basis_vector[channels], "*")
    }

    effect_array <- effect_array + add_basis_effect(first_obj, first_score_diff)
    for (file in fit_files[-1]) {
        obj <- readRDS(file)
        effect_array <- effect_array +
            add_basis_effect(obj, condition_score_draws(obj, times))
    }

    list(
        tag = result_tag,
        gt = ground_truth,
        sim_model = sim_model,
        conv_summary = conv_summary,
        effect_array = effect_array,
        times = times,
        channels = channels
    )
}
