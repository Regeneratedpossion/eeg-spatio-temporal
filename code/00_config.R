# Settings shared by the thesis scripts (tables, figures and figure data).
# Run every script from the project root, the folder that holds the .here file.
#   PROJECT_ROOT  overrides the project folder (default: found with here::here())
#   THESIS_DIR    folder with the small input files and the outputs of tables and figures (default: results/)
#   FIG_OUT       if set, figures and csv files are written there instead of their usual folder
if (!exists("thesis_dir")) {
  project_dir <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = here::here()))
  thesis_dir  <- normalizePath(Sys.getenv("THESIS_DIR", unset = file.path(project_dir, "results")))
  fig_out     <- Sys.getenv("FIG_OUT", unset = "")
  thesis_path  <- function(...) file.path(thesis_dir, ...)
  project_path <- function(...) file.path(project_dir, ...)
  fig_dir <- function(default) {
    d <- if (nzchar(fig_out)) fig_out else default
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
    d
  }
  # condition letter of each scenario: 1 = I, 2 = H, 3 = F, 4 = E
  scenario_cond <- c("1" = "I", "2" = "H", "3" = "F", "4" = "E")
}
