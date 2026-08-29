#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) {
  stop("Usage: Rscript make_parallel_nc_grid.R <best_csv> <out_csv>", call. = FALSE)
}

best_file <- args[[1L]]
out_file <- args[[2L]]
grid <- utils::read.csv(best_file, stringsAsFactors = FALSE, check.names = FALSE)
needed <- c(
  "mae", "mae_file", "fusion_mode", "sis_n", "glmnet_alpha", "outcome",
  "n_active_mediators", "x_var", "y_family", "covariates",
  "method_preset", "screen_method", "lambda_choice", "b_inference",
  "fdr_method", "fdr_scope"
)
missing_needed <- setdiff(needed, colnames(grid))
if (length(missing_needed) > 0L) {
  stop("Best-results table is missing columns: ", paste(missing_needed, collapse = ", "), call. = FALSE)
}
grid <- grid[, needed, drop = FALSE]
grid$job_id <- seq_len(nrow(grid))
utils::write.csv(grid, out_file, row.names = FALSE)
message("Wrote parallel NC grid: ", out_file)
message("Rows: ", nrow(grid))
