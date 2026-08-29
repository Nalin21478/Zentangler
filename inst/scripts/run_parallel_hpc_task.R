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
  if (length(val) == 0L || is.na(val) || (!is.numeric(val) && !nzchar(as.character(val)))) return(default)
  val
}

split_arg <- function(x) {
  if (is.null(x) || length(x) == 0L || all(is.na(x)) || !nzchar(paste(x, collapse = ""))) return(character(0))
  out <- trimws(unlist(strsplit(as.character(x), ",", fixed = TRUE), use.names = FALSE))
  out[nzchar(out)]
}

safe_tag <- function(x) {
  x <- ifelse(is.na(x) | !nzchar(as.character(x)), "NA", as.character(x))
  gsub("[^A-Za-z0-9]+", "", x)
}

load_mae <- function(mae_file) {
  if (!file.exists(mae_file)) stop("MAE file not found: ", mae_file, call. = FALSE)
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

out_dir <- if (!is.null(args$`out-dir`)) args$`out-dir` else "parallel_hpc_results"
task_offset <- if (!is.null(args$`task-offset`)) as.integer(args$`task-offset`) else 0L
seed <- if (!is.null(args$seed)) as.integer(args$seed) else 1L

grid <- utils::read.csv(args$grid, stringsAsFactors = FALSE, check.names = FALSE)
array_task_id <- as.integer(args$`task-id`)
task_id <- array_task_id + task_offset
if (!is.finite(task_id) || task_id < 1L || task_id > nrow(grid)) {
  stop("effective task-id must be between 1 and ", nrow(grid), ". Got ", task_id, call. = FALSE)
}

row <- grid[task_id, , drop = FALSE]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(Zentangler)
  library(MultiAssayExperiment)
})

mae_name <- as.character(row_value(row, "mae"))
mae_file <- as.character(row_value(row, "mae_file"))
fusion_mode <- as.character(row_value(row, "fusion_mode"))
sis_n <- as.integer(row_value(row, "sis_n"))
glmnet_alpha <- as.numeric(row_value(row, "glmnet_alpha"))
outcome <- as.character(row_value(row, "outcome"))
x_var <- as.character(row_value(row, "x_var", "treatment_any"))
y_family <- as.character(row_value(row, "y_family", "gaussian"))
covars <- split_arg(row_value(row, "covariates", ""))
method_preset <- as.character(row_value(row, "method_preset", "custom"))
screen_method <- as.character(row_value(row, "screen_method", "sis"))
lambda_choice <- as.character(row_value(row, "lambda_choice", "lambda.min"))
b_inference <- as.character(row_value(row, "b_inference", "debiased_lasso"))
fdr_method <- as.character(row_value(row, "fdr_method", "BH"))
fdr_scope <- as.character(row_value(row, "fdr_scope", "global"))

tag <- sprintf(
  "job%04d_%s_%s_sis%s_alpha%s_%s",
  as.integer(row_value(row, "job_id", task_id)),
  safe_tag(mae_name),
  safe_tag(fusion_mode),
  safe_tag(sis_n),
  gsub("\\.", "p", as.character(glmnet_alpha)),
  safe_tag(outcome)
)

out_rds <- file.path(out_dir, paste0(tag, ".rds"))
out_summary <- file.path(out_dir, paste0(tag, "_observed_summary.csv"))
out_active <- file.path(out_dir, paste0(tag, "_observed_mediators_active.csv"))
out_grid_row <- file.path(out_dir, paste0(tag, "_grid_row.csv"))
out_session <- file.path(out_dir, paste0(tag, "_sessionInfo.txt"))

cat("Running parallel HPC task\n")
cat("Task:", task_id, "of", nrow(grid), "\n")
print(row)

mae_obj <- load_mae(mae_file)
utils::write.csv(row, out_grid_row, row.names = FALSE)
writeLines(capture.output(sessionInfo()), out_session)

fit <- fit_multiview_parallel_zentangler(
  mae = mae_obj,
  x_var = x_var,
  y_var = outcome,
  y_family = y_family,
  covariates = covars,
  method_preset = method_preset,
  screen_method = screen_method,
  sis_n = sis_n,
  glmnet_alpha = glmnet_alpha,
  lambda_choice = lambda_choice,
  fusion_mode = fusion_mode,
  b_inference = b_inference,
  fdr_method = fdr_method,
  fdr_scope = fdr_scope,
  seed = seed
)

active_tab <- fit[["mediators_active"]]
observed_summary <- data.frame(
  mae = mae_name,
  mae_file = mae_file,
  fusion_mode = fusion_mode,
  sis_n = sis_n,
  glmnet_alpha = glmnet_alpha,
  outcome = outcome,
  x_var = x_var,
  y_family = y_family,
  covariates = paste(covars, collapse = ","),
  method_preset = method_preset,
  screen_method = screen_method,
  lambda_choice = lambda_choice,
  b_inference = b_inference,
  fdr_method = fdr_method,
  fdr_scope = fdr_scope,
  n_active_mediators = if (is.null(active_tab)) 0L else nrow(active_tab),
  stringsAsFactors = FALSE
)

saveRDS(fit, out_rds)
write_if_nonempty(observed_summary, out_summary)
write_if_nonempty(active_tab, out_active)

cat("Saved RDS:", out_rds, "\n")
cat("Done.\n")
