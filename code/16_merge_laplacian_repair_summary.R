###########################################################
# 16_merge_laplacian_repair_summary.R
#
# Merge a small set of repaired Laplacian basis fits into a full
# 1:K summary so downstream reconstruction diagnostics can use
# the original good basis fits plus repaired problematic basis fits.
###########################################################

source(file.path(Sys.getenv("PROJECT_ROOT", unset = getwd()), "code/utils.R"))

base_subdir <- Sys.getenv("LAPLACE_BASE_SUBDIR")
repair_subdir <- Sys.getenv("LAPLACE_REPAIR_SUBDIR")
out_subdir <- Sys.getenv("LAPLACE_MERGED_SUBDIR")
if (!nzchar(base_subdir) || !nzchar(repair_subdir) || !nzchar(out_subdir)) {
    stop("Set LAPLACE_BASE_SUBDIR, LAPLACE_REPAIR_SUBDIR, and LAPLACE_MERGED_SUBDIR.")
}

iter <- as.integer(Sys.getenv("LAPLACE_ITER", "3000"))
warmup <- as.integer(Sys.getenv("LAPLACE_WARMUP", "1000"))
base_adapt_tag <- Sys.getenv("LAPLACE_BASE_ADAPT_TAG", "0p999")
repair_adapt_tag <- Sys.getenv("LAPLACE_REPAIR_ADAPT_TAG", "0p9999")
merged_adapt_tag <- Sys.getenv("LAPLACE_MERGED_ADAPT_TAG", repair_adapt_tag)
repair_tasks <- as.integer(strsplit(Sys.getenv("LAPLACE_REPAIR_TASKS", "11,19"),
                                    ",")[[1]])

summary_name <- function(adapt_tag) {
    sprintf("summary_all_iter%d_warmup%d_adapt%s.csv", iter, warmup, adapt_tag)
}

base_file <- file.path(project_root(), "results", base_subdir,
                       summary_name(base_adapt_tag))
repair_file <- file.path(project_root(), "results", repair_subdir,
                         summary_name(repair_adapt_tag))
out_dir <- file.path(project_root(), "results", out_subdir)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(out_dir, summary_name(merged_adapt_tag))

base <- read.csv(base_file, stringsAsFactors = FALSE)
repair <- read.csv(repair_file, stringsAsFactors = FALSE)

missing_repair <- setdiff(repair_tasks, repair$basis_task)
if (length(missing_repair)) {
    stop(sprintf("Repair summary missing basis tasks: %s",
                 paste(missing_repair, collapse = ",")))
}

merged <- base[!(base$basis_task %in% repair_tasks), , drop = FALSE]
merged <- rbind(merged, repair[repair$basis_task %in% repair_tasks, ,
                               drop = FALSE])
merged <- merged[order(merged$basis_task), , drop = FALSE]

if (!identical(merged$basis_task, seq_len(max(merged$basis_task)))) {
    stop("Merged summary does not contain a complete ordered 1:K basis set.")
}
if (!all(file.exists(merged$fit_file))) {
    stop("Some merged fit_file paths do not exist.")
}

write.csv(merged, out_file, row.names = FALSE)
cat("Saved:", out_file, "\n")
print(merged[, c("basis_task", "basis_label", "n_divergent",
                 "n_treedepth", "max_rhat", "min_ess_bulk",
                 "min_ess_tail", "fit_file")])
