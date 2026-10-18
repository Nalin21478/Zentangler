#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

parse_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!grepl("^--", key)) {
      i <- i + 1L
      next
    }
    nm <- sub("^--", "", key)
    val <- if (i + 1L <= length(args) && !grepl("^--", args[[i + 1L]])) args[[i + 1L]] else TRUE
    out[[nm]] <- val
    i <- i + if (isTRUE(val)) 1L else 2L
  }
  out
}

row_value <- function(row, name, default = NULL) {
  if (!(name %in% colnames(row))) return(default)
  val <- row[[name]][1]
  if (length(val) == 0L || is.na(val) || !nzchar(as.character(val))) return(default)
  val
}

safe_tag <- function(x) {
  x <- ifelse(is.na(x) | !nzchar(as.character(x)), "NA", as.character(x))
  gsub("[^A-Za-z0-9]+", "", x)
}

mae_path_from_grid <- function(mae_name, mae_dir) {
  file.path(
    mae_dir,
    switch(
      as.character(mae_name),
      pathways = "Blueberry_MAE_pathways.RData",
      ecs = "Blueberry_MAE_ecs.RData",
      all_omics = "Blueberry_MAE_all_omics.RData",
      stop("Unknown mae value: ", mae_name, call. = FALSE)
    )
  )
}

load_blueberry_mae <- function(mae_name, mae_dir) {
  mae_file <- mae_path_from_grid(mae_name, mae_dir)
  if (!file.exists(mae_file)) {
    stop("MAE file not found: ", mae_file, call. = FALSE)
  }
  env <- new.env(parent = emptyenv())
  load(mae_file, envir = env)
  if (!exists("mae", envir = env, inherits = FALSE)) {
    stop("Expected object 'mae' in ", mae_file, call. = FALSE)
  }
  get("mae", envir = env, inherits = FALSE)
}

write_if_nonempty <- function(x, file) {
  if (is.null(x)) return(invisible(FALSE))
  if (is.data.frame(x) && nrow(x) == 0L) return(invisible(FALSE))
  utils::write.csv(x, file, row.names = FALSE)
  invisible(TRUE)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (is.null(args$grid)) stop("Missing --grid", call. = FALSE)
if (is.null(args$`task-id`)) stop("Missing --task-id", call. = FALSE)

out_dir <- if (!is.null(args$`out-dir`)) args$`out-dir` else "blueberry_negative_control_results"
mae_dir <- if (!is.null(args$`mae-dir`)) args$`mae-dir` else file.path("data")
task_offset <- if (!is.null(args$`task-offset`)) as.integer(args$`task-offset`) else 0L
n_perm <- if (!is.null(args$`n-perm`)) as.integer(args$`n-perm`) else 500L
q_threshold <- if (!is.null(args$`q-threshold`)) as.numeric(args$`q-threshold`) else 0.25
seed <- if (!is.null(args$seed)) as.integer(args$seed) else 1L

grid <- utils::read.csv(args$grid, stringsAsFactors = FALSE, check.names = FALSE)
array_task_id <- as.integer(args$`task-id`)
task_id <- array_task_id + task_offset
if (!is.finite(task_id) || task_id < 1L || task_id > nrow(grid)) {
  stop(
    "effective task-id must be between 1 and ", nrow(grid),
    ". Got array task-id ", array_task_id, " + offset ", task_offset,
    " = ", task_id,
    call. = FALSE
  )
}

row <- grid[task_id, , drop = FALSE]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Map a fusion label to the zentangler outcome-model arguments:
#   "multiview" (default; cooperative learning, rho in [0, 1] tuned by CV)
#   "early" / "late" (supplementary IntegratedLearner-style fusion)
# The legacy grid value "intermediate" is read as "multiview".
fusion_args <- function(label) {
  label <- as.character(label)[1]
  if (label %in% c("multiview", "intermediate", "cooperative")) {
    return(list(fusion = "multiview", integrated_learner_fusion = "late"))
  }
  if (label %in% c("early", "late")) {
    return(list(fusion = "integrated_learner", integrated_learner_fusion = label))
  }
  stop("Unknown fusion label: ", label, call. = FALSE)
}

suppressPackageStartupMessages({
  library(Zentangler)
  library(MultiAssayExperiment)
})

mae_name <- as.character(row_value(row, "mae"))
# The grid keeps its fusion_mode column (early / late / multiview) so existing
# CSVs still work; it is translated to fusion / integrated_learner_fusion here.
fusion_mode <- as.character(row_value(row, "fusion_mode"))
fusion_settings <- fusion_args(fusion_mode)
sis_n <- as.integer(row_value(row, "sis_n"))
glmnet_alpha <- as.numeric(row_value(row, "glmnet_alpha"))
outcome <- as.character(row_value(row, "outcome"))
n_active_mediators <- as.integer(row_value(row, "n_active_mediators", NA))

covars <- c(
  "Age",
  "sex",
  "bmi",
  "bp_med",
  "statin",
  "metabolicsyndromecriteria"
)

tag <- sprintf(
  "job%03d_%s_%s_sis%s_alpha%s_%s_perm%s",
  as.integer(row_value(row, "job_id", task_id)),
  safe_tag(mae_name),
  safe_tag(fusion_mode),
  safe_tag(sis_n),
  gsub("\\.", "p", as.character(glmnet_alpha)),
  safe_tag(outcome),
  safe_tag(n_perm)
)

out_rds <- file.path(out_dir, paste0(tag, ".rds"))
out_summary <- file.path(out_dir, paste0(tag, "_empirical_summary.csv"))
out_null <- file.path(out_dir, paste0(tag, "_null_summaries.csv"))
out_observed <- file.path(out_dir, paste0(tag, "_observed_summary.csv"))
out_mediators <- file.path(out_dir, paste0(tag, "_observed_mediators_active.csv"))
out_empirical_mediators <- file.path(out_dir, paste0(tag, "_empirical_mediator_summary.csv"))
out_grid_row <- file.path(out_dir, paste0(tag, "_grid_row.csv"))
out_session <- file.path(out_dir, paste0(tag, "_sessionInfo.txt"))

cat("Running Blueberry negative-control task\n")
cat("Task:", task_id, "of", nrow(grid), "\n")
print(row)

mae_obj <- load_blueberry_mae(mae_name = mae_name, mae_dir = mae_dir)
utils::write.csv(row, out_grid_row, row.names = FALSE)
writeLines(capture.output(sessionInfo()), out_session)

fit <- run_negative_controls(
  mae = mae_obj,
  x_var = "treatment_any",
  y_var = outcome,
  n_perm = n_perm,
  q_threshold = q_threshold,
  seed = seed,
  y_family = "gaussian",
  covariates = covars,
  screen_method = "sis",
  sis_n = sis_n,
  glmnet_alpha = glmnet_alpha,
  lambda_choice = "lambda.min",
  fusion = fusion_settings$fusion,
  integrated_learner_fusion = fusion_settings$integrated_learner_fusion,
  b_inference = "debiased",
  fdr_method = "BH",
  fdr_scope = "global"
)

saveRDS(fit, out_rds)
write_if_nonempty(fit$empirical_summary, out_summary)
write_if_nonempty(fit$null_summaries, out_null)
write_if_nonempty(fit$observed_summary, out_observed)
write_if_nonempty(zentangler_mediators(fit$observed_fit, q_threshold = q_threshold), out_mediators)
write_if_nonempty(fit$empirical_units, out_empirical_mediators)
# Per-permutation null mediator tables are no longer returned; per-mediator null
# frequencies and p-values are in fit$empirical_units (written above).

cat("Saved RDS:", out_rds, "\n")
cat("Observed active mediators:", n_active_mediators, "\n")
cat("Done.\n")
