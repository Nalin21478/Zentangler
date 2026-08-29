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

split_arg <- function(x, default = character(0)) {
  if (is.null(x) || length(x) == 0L || all(is.na(x)) || !nzchar(paste(x, collapse = ""))) return(default)
  out <- trimws(unlist(strsplit(as.character(x), ",", fixed = TRUE), use.names = FALSE))
  out[nzchar(out)]
}

as_num_vec <- function(x, default = numeric(0)) {
  vals <- suppressWarnings(as.numeric(split_arg(x)))
  vals <- vals[is.finite(vals)]
  if (length(vals) == 0L) default else vals
}

manifest <- function(manifest_file) {
  tab <- utils::read.csv(manifest_file, stringsAsFactors = FALSE, check.names = FALSE)
  needed <- c("mae", "mae_file")
  missing_needed <- setdiff(needed, colnames(tab))
  if (length(missing_needed) > 0L) {
    stop("Manifest is missing columns: ", paste(missing_needed, collapse = ", "), call. = FALSE)
  }
  tab$mae <- as.character(tab$mae)
  tab$mae_file <- as.character(tab$mae_file)
  tab <- tab[nzchar(tab$mae) & nzchar(tab$mae_file), , drop = FALSE]
  if (nrow(tab) == 0L) stop("Manifest contains no usable MAE rows.", call. = FALSE)
  rownames(tab) <- NULL
  tab
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (is.null(args$manifest)) stop("Missing --manifest", call. = FALSE)
if (is.null(args$out)) stop("Missing --out", call. = FALSE)

mae_tab <- manifest(args$manifest)
outcomes <- split_arg(args$outcomes)
fusion_modes <- split_arg(args$fusion, default = c("early", "late"))
sis_grid <- as_num_vec(args$sis, default = c(10, 25, 50, 100))
alpha_grid <- as_num_vec(args$alpha, default = c(0, 0.5, 1))
x_var <- if (!is.null(args$`x-var`)) as.character(args$`x-var`) else "treatment_any"
y_family <- if (!is.null(args$`y-family`)) as.character(args$`y-family`) else "gaussian"
covariates <- paste(split_arg(args$covariates), collapse = ",")
screen_method <- if (!is.null(args$`screen-method`)) as.character(args$`screen-method`) else "sis"
method_preset <- if (!is.null(args$`method-preset`)) as.character(args$`method-preset`) else "custom"
lambda_choice <- if (!is.null(args$`lambda-choice`)) as.character(args$`lambda-choice`) else "lambda.min"
b_inference <- if (!is.null(args$`b-inference`)) as.character(args$`b-inference`) else "debiased_lasso"
fdr_method <- if (!is.null(args$`fdr-method`)) as.character(args$`fdr-method`) else "BH"
fdr_scope <- if (!is.null(args$`fdr-scope`)) as.character(args$`fdr-scope`) else "global"

if (length(outcomes) == 0L) stop("Provide at least one outcome via --outcomes.", call. = FALSE)

grid <- expand.grid(
  mae = mae_tab$mae,
  fusion_mode = fusion_modes,
  sis_n = sis_grid,
  glmnet_alpha = alpha_grid,
  outcome = outcomes,
  stringsAsFactors = FALSE,
  KEEP.OUT.ATTRS = FALSE
)
grid <- merge(grid, mae_tab, by = "mae", all.x = TRUE, sort = FALSE)
grid$x_var <- x_var
grid$y_family <- y_family
grid$covariates <- covariates
grid$screen_method <- screen_method
grid$method_preset <- method_preset
grid$lambda_choice <- lambda_choice
grid$b_inference <- b_inference
grid$fdr_method <- fdr_method
grid$fdr_scope <- fdr_scope
grid$job_id <- seq_len(nrow(grid))
grid <- grid[, c(
  "job_id", "mae", "mae_file", "x_var", "outcome", "y_family",
  "covariates", "method_preset", "screen_method", "sis_n",
  "glmnet_alpha", "lambda_choice", "fusion_mode", "b_inference",
  "fdr_method", "fdr_scope"
), drop = FALSE]

utils::write.csv(grid, args$out, row.names = FALSE)
message("Wrote parallel HPC grid: ", args$out)
message("Rows: ", nrow(grid))
